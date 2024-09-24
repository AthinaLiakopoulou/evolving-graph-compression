package gr.uoa.di.networkanalysis;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

import it.unimi.dsi.bits.Fast;
import it.unimi.dsi.fastutil.io.BinIO;
import it.unimi.dsi.fastutil.io.FastMultiByteArrayInputStream;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.io.InputBitStream;
import it.unimi.dsi.io.OutputBitStream;
import it.unimi.dsi.sux4j.util.EliasFanoMonotoneLongBigList;
import it.unimi.dsi.webgraph.LazyIntIterator;
import me.lemire.integercompression.differential.IntegratedIntCompressor;

public class EvolvingMultiGraph {

    protected String graphFile;
    protected boolean headers;
    protected String basename;
    protected double aggregationFactor;

    protected BVMultiGraph graph;
    protected EliasFanoMonotoneLongBigList efindex;
    protected byte[] timestamps;
    protected FastMultiByteArrayInputStream timestampsStream;
    protected long minTimestamp;

    LongArrayList offsetsIndex = null;
    private long currentOffset;
    private final IntegratedIntCompressor compressor;
    DebugLogger encode_debugLogger = new DebugLogger("debug_encode.txt");
    DebugLogger decode_debugLogger = new DebugLogger("debug_decode.txt");
    DebugLogger percentages_debugLogger = new DebugLogger("percentages.txt");
    private int smallListCount = 0;    // Λίστες μικρότερες των 128
    private int compressedCount = 0;   // Συμπιεσμένες λίστες
    private int uncompressedCount = 0; // Μη συμπιεσμένες λίστες


    public EvolvingMultiGraph(String graphFile, boolean headers, String basename, double aggregationFactor) {
        super();
        this.graphFile = graphFile;
        this.headers = headers;
        this.basename = basename;
        this.aggregationFactor = aggregationFactor;
        this.compressor = new IntegratedIntCompressor();
    }

    protected long findMinimumTimestamp() {
        try (InputStream fileStream = new FileInputStream(graphFile);
                InputStream gzipStream = new GZIPInputStream(fileStream);
                Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
                BufferedReader buffered = new BufferedReader(decoder);
                ) {
            Long minTimestamp = null;
            minTimestamp = buffered.lines().mapToLong(line -> Long.parseLong(line.split("\t")[3])).min().orElse(0);
            return minTimestamp;
        } catch (IOException e) {
        }
        return 0;
    }

    protected long writeTimestampsToFile(List<Long> currentNeighborsTimestamps, OutputBitStream obs, long minTimestamp) throws IOException {
        // Sort timestamps in ascending order
        currentNeighborsTimestamps.sort(Long::compareTo);

        // Returns the number of bits appended to the file
        long ret = 0;
        long previousNeighborTimestamp = minTimestamp;
        ArrayList<Integer> periodsBetweenList = new ArrayList<>();

        // Calculate the periods between the current and previous timestamps
        for(Long seconds : currentNeighborsTimestamps) {
            double periodsBetween = TimestampComparerAggregator.timestampsDifference(previousNeighborTimestamp, seconds, aggregationFactor);
            periodsBetween = Fast.int2nat(Math.round(periodsBetween));
            periodsBetweenList.add((int) periodsBetween);  // Store the period
            previousNeighborTimestamp = seconds;
        }

        // Convert list to int array
        int[] periodsArray = periodsBetweenList.stream().mapToInt(i -> i).toArray();

        // If the size of the list is smaller than 128, increase the counter
        if (periodsArray.length < 128) {
            smallListCount++;
        }

        // Check if compression is beneficial
        int[] compressedData = compressor.compress(periodsArray);
        int uncompressedSizeBits = periodsArray.length * 32; // Total bits if uncompressed
        int compressedSizeBits = compressedData.length * 32; // Assuming each int is 32 bits

        // If compression is not effective or array is too small, store uncompressed
        if (periodsArray.length < 128 || compressedSizeBits >= uncompressedSizeBits) {
            obs.writeBit(0); // 0 indicates uncompressed
            obs.writeInt(periodsArray.length, 32); // Write the number of periods (to ensure correct reading)
            for (int period : periodsArray) {
                obs.writeInt(period, 32); // Write each period as an int (uncompressed)
            }
            ret += 1 + 32 + uncompressedSizeBits; // 1 bit for flag + 32 bits for length + uncompressed size

            // Increase the counter for uncompressed data
            uncompressedCount++;

            encode_debugLogger.log("Stored uncompressed data for array length: " + periodsArray.length);
            return ret;
        }

        // If compression is beneficial, store compressed data
        obs.writeBit(1); // 1 indicates compressed
        obs.writeInt(compressedSizeBits, 32); // Write compressed size in bits
        for (int value : compressedData) {
            obs.writeInt(value, 32); // Write each int as 32 bits
        }

        ret += 1 + 32 + compressedSizeBits; // 1 bit for flag + 32 bits for size + compressed size

        // Increase the counter for compressed data
        compressedCount++;

        encode_debugLogger.log("Stored compressed data, original length: " + periodsArray.length + ", compressed length: " + compressedData.length);
        return ret;
    }

    protected int[] readTimestampsFromFile(InputBitStream ibs) throws IOException {
        // Read the compression flag
        boolean isCompressed = ibs.readBit() == 1;
        int[] decompressedData;

        if (isCompressed) {
            // Read the size of the compressed data
            int compressedSize = ibs.readInt(32); // Size in bits
            int numCompressedInts = (compressedSize + 31) / 32; // Number of 32-bit ints
            int[] compressedData = new int[numCompressedInts];

            // Read the compressed data
            for (int i = 0; i < numCompressedInts; i++) {
                compressedData[i] = ibs.readInt(32); // Each int is 32 bits
            }

            // Decompress the data
            decompressedData = compressor.uncompress(compressedData);
        } else {
            // Read uncompressed data directly
            List<Integer> uncompressedList = new ArrayList<>();
            int length = ibs.readInt(32); // Number of uncompressed timestamps

            for (int i = 0; i < length; i++) {
                if (!ibs.hasNext()) {
                    throw new EOFException("Unexpected end of file while reading uncompressed timestamps.");
                }
                uncompressedList.add(ibs.readInt(32)); // Read each int
            }

            decompressedData = uncompressedList.stream().mapToInt(i -> i).toArray();
        }

        return decompressedData;
    }


    public void store() throws IOException, InterruptedException {
        ExecutorService executorService = Executors.newFixedThreadPool(2);
        executorService.execute(()-> {storeBVMultiGraph();});
        executorService.execute(()-> {
            storeTimestampsAndIndex();
            percentages_debugLogger.log("Sum of lists where lenght < 128: " + smallListCount);
            percentages_debugLogger.log("Compressed lists: " + compressedCount);
            percentages_debugLogger.log("Uncompressed lists: " + uncompressedCount);

            // Υπολογισμός ποσοστών
            int totalLists = compressedCount + uncompressedCount;
            if (totalLists > 0) {
                double compressedPercentage = (compressedCount * 100.0) / totalLists;
                double uncompressedPercentage = (uncompressedCount * 100.0) / totalLists;
                percentages_debugLogger.log("Percentage of compressed data: " + compressedPercentage + "%");
                percentages_debugLogger.log("Percentage of uncompressed data: " + uncompressedPercentage + "%");
            }
            percentages_debugLogger.close();
        });
        executorService.shutdown();
        executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
    }

    public void storeBVMultiGraph() {
        try (InputStream fileStream = new FileInputStream(graphFile);
                InputStream gzipStream = new GZIPInputStream(fileStream);
                ) {
            ArcListASCIIEvolvingGraph inputGraph = new ArcListASCIIEvolvingGraph(gzipStream, 0);
            BVMultiGraph.store(inputGraph, basename);
        } catch (FileNotFoundException e) {
        } catch (IOException e) {
        }

    }

    public void storeTimestampsAndIndex() {

        // Aggregate the minimum timestamp of the file
        minTimestamp = TimestampComparerAggregator.aggregateMinTimestamp(findMinimumTimestamp(), aggregationFactor);

        try (InputStream fileStream = new FileInputStream(graphFile);
                GZIPInputStream gzipStream = new GZIPInputStream(fileStream);
                InputStreamReader decoder = new InputStreamReader(gzipStream, "UTF-8");
                BufferedReader buffered = new BufferedReader(decoder);) {
            if(headers) {
                buffered.readLine();
            }
            // The file we will write the results to
            final OutputBitStream obs = new OutputBitStream(new FileOutputStream(basename+".timestamps"), 1024 * 1024);
            offsetsIndex= new LongArrayList();
            currentOffset = obs.writeLong(minTimestamp, 64);

            int currentNode = 0;
            ArrayList<Long> currentNeighborsTimestamps = new ArrayList<>();
            String line;

            ExecutorService executorService = Executors.newSingleThreadExecutor();

            while ((line = buffered.readLine()) != null) {
                //String[] tokens = line.split("\\s+");
                String[] tokens = line.split("\t");
                int node = Integer.parseInt(tokens[0]);
                long timestamp = Long.parseLong(tokens[3]);

                int previous = currentNode;

                // If you find a new currentNode in the file, write the results you have so far about the current node.
                if(node != currentNode) {
                    executorService.submit(new TimeStampsWriter(currentNeighborsTimestamps, obs, node, previous));

                    // Prepare the variables for the next currentNode
                    currentNode = node;
                    currentNeighborsTimestamps = new ArrayList<Long>();
                    currentNeighborsTimestamps.add(timestamp);

                }
                else {
                    currentNeighborsTimestamps.add(timestamp);
                }
            }
            // Write the last node. It was not written because no change in node != currentNode was detected
            executorService.submit(new TimeStampsWriter(currentNeighborsTimestamps, obs, 0, 0));

            executorService.shutdown();
            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);

            obs.close();
            buffered.close();

            // Perform compression of the index using EliasFano
            EliasFanoMonotoneLongBigList efmlbl = new EliasFanoMonotoneLongBigList(offsetsIndex);
            FileOutputStream fos = new FileOutputStream(basename+".efindex");
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(efmlbl);
            oos.close();
            fos.close();

            offsetsIndex = null;
        } catch (IOException e) {
        } catch (InterruptedException e) {
        }

    }

    public void load() throws Exception {
        // Graph
        graph = BVMultiGraph.load(basename);
        // EliasFano index
        FileInputStream fis = new FileInputStream(basename+".efindex");
        ObjectInputStream ois = new ObjectInputStream(fis);
        efindex = (EliasFanoMonotoneLongBigList) ois.readObject();
        ois.close();
        // Timestamps
        fis = new FileInputStream(basename+".timestamps");
        if(fis.getChannel().size() <= Integer.MAX_VALUE) {
            timestamps = new byte[(int) fis.getChannel().size()];
            BinIO.loadBytes(fis, timestamps);
            fis.close();
            InputBitStream ibs = new InputBitStream(timestamps);
            minTimestamp = ibs.readLong(64);
            ibs.close();
        } else {
            timestampsStream = new FastMultiByteArrayInputStream(fis, fis.getChannel().size());
            InputBitStream ibs = new InputBitStream(timestampsStream);
            minTimestamp = ibs.readLong(64);
            ibs.close();
        }

    }

    public boolean isNeighbor(int node, int neighbor) {
        LazyIntIterator it = graph.successors(node);
        int n = -1;
        while((n = it.nextInt()) != -1) {
			if(n == neighbor) return true;
			else if(n > neighbor) return false;
        }

        return false;
    }

    public boolean isNeighbor(int node, int neighbor, long t1, long t2) throws Exception {
        LazyIntIterator it = graph.successors(node);
        int n = -1, from = -1, to = -1, pos = 0;
        long t;
        // Find  the starting position
        while((n = it.nextInt()) != -1) {
            if(n == neighbor) {
                from = pos++;
                break;
            }
        }
        // Return false if it was not found at least once
        if(from == -1) {
            return false;
        }
        // Find the ending position
        while((n = it.nextInt()) == neighbor) {
            pos++;
        }
        to = pos-1;

        InputBitStream ibs;
        if (timestamps != null) {
            ibs = new InputBitStream(timestamps);
        } else {
            ibs = new InputBitStream(new FastMultiByteArrayInputStream(timestampsStream));
        }

        ibs.position(efindex.getLong(node));
        int[] decompressedData = readTimestampsFromFile(ibs);
        decode_debugLogger.log("Decompressed Data Length: " + decompressedData.length);

        // Close the InputBitStream
        ibs.close();

        // Initialize index for accessing decompressed data
        int index = 0;

        if (from < 0 || to >= decompressedData.length || from > to) {
            return false;
        }

        // Skip everything up to 'from'
        long previous = minTimestamp;
        for (int i = 0; i < from; i++) {
            t = TimestampComparerAggregator.reverse(previous, decompressedData[index++], aggregationFactor);
            previous = t;
        }

        // Scan all the timestamps in the range 'from' to 'to'
        for (int i = from; i <= to; i++) {
            t = TimestampComparerAggregator.reverse(previous, decompressedData[index++], aggregationFactor);
            if (t1 <= t && t <= t2) {
                return true;
            }
            previous = t;
        }

        return false;
    }

    public SuccessorIterator successors(int node) throws Exception {
        return new SuccessorIterator(node);
    }

    public class SuccessorIterator implements Iterator<Successor> {

        LazyIntIterator neighborsIterator;
        InputBitStream ibs;
        long previous;
        int[] decompressedData; // To store decompressed timestamps
        int index = 0;

        public SuccessorIterator(int node) throws Exception {
            neighborsIterator = graph.successors(node);
            if (timestamps != null) {
                ibs = new InputBitStream(timestamps);
            } else {
                ibs = new InputBitStream(new FastMultiByteArrayInputStream(timestampsStream));
            }
            long position = efindex.getLong(node);
            ibs.position(position);

            decompressedData = readTimestampsFromFile(ibs);
            decode_debugLogger.log("Decompressed Data Length: " + decompressedData.length);

            previous = minTimestamp;
            ibs.close(); // Close the stream after reading the timestamps
        }

        @Override
        public boolean hasNext() throws UnsupportedOperationException {
            throw new UnsupportedOperationException();
        }

        @Override
        public Successor next() throws NoSuchElementException {
            int neighbor = neighborsIterator.nextInt();
            if(neighbor == -1) {
                throw new NoSuchElementException();
            }
            long t;
            if (index < decompressedData.length) {
                t = Fast.nat2int(decompressedData[index]);
                t = TimestampComparerAggregator.reverse(previous, t, aggregationFactor);
                previous = t;
                index++;
            } else {
                throw new NoSuchElementException("No more timestamps available");
            }
            return new Successor(neighbor, t);
        }
    }

    public String getGraphFile() {
        return graphFile;
    }

    public boolean isHeaders() {
        return headers;
    }

    public String getBasename() {
        return basename;
    }

    public BVMultiGraph getGraph() {
        return graph;
    }

    public EliasFanoMonotoneLongBigList getEfindex() {
        return efindex;
    }

    public byte[] getTimestamps() {
        return timestamps;
    }

    public long getMinTimestamp() {
        return minTimestamp;
    }

    class TimeStampsWriter implements Runnable {


        private List<Long> currentNeighborsTimestamps;
        private OutputBitStream obs;
        private int node;
        private int previous;

        public TimeStampsWriter(List<Long> currentNeighborsTimestamps, OutputBitStream obs, int node, int previous) {
            this.currentNeighborsTimestamps = currentNeighborsTimestamps;
            this.obs = obs;
            this.node = node;
            this.previous = previous;
        }

        @Override
        public void run() {
            offsetsIndex.add(currentOffset);
            try {
                currentOffset += writeTimestampsToFile(currentNeighborsTimestamps, obs, minTimestamp);
            } catch (IOException e) {
            }
            // If at least one node was skipped, add the current offset to the index for each node in-between
            if(node > previous + 1) {
                for(int i = 0; i < node - previous - 1; i++) {
                    offsetsIndex.add(currentOffset);
                }
            }
        }

    }


}
