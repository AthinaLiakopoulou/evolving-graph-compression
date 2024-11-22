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
import me.lemire.integercompression.*;
import me.lemire.integercompression.differential.IntegratedIntCompressor;

public class EvolvingMultiGraph {

    protected String graphFile;
    protected boolean headers;
    protected String basename;
    protected long aggregationFactor;

    protected BVMultiGraph graph;
    protected EliasFanoMonotoneLongBigList efindex;
    protected byte[] timestamps;
    protected FastMultiByteArrayInputStream timestampsStream;
    protected long minTimestamp;

    LongArrayList offsetsIndex = null;
    private long currentOffset;

    private static final int DELTA_THRESHOLD = 10;

    public class CombinedCompressor {
        private final IntegratedIntCompressor integratedCompressor;
        private final FastPFOR fastPForCompressor;
        private final VariableByte variableByte;

        public CombinedCompressor() {
            this.integratedCompressor = new IntegratedIntCompressor();
            this.fastPForCompressor = new FastPFOR();
            this.variableByte = new VariableByte();
        }

        //compress returns compressor flag, compressed data
        public int[] compress(int[] data) {
            // Check if the data is suitable for delta compression
            if (isSmallDeltas(data)) {
                // Use IntegratedIntCompressor for delta compression
                int[] compressedData = integratedCompressor.compress(data);
                return prependCompressedSizeUncompressedSizeAndIndicator(compressedData, 0,data.length); // Prepend 0 for IntegratedIntCompressor
            } else {
                // Use FastPFOR for standard compression
                int[] compressedData = compressWithFastPFor(data);
                return prependCompressedSizeUncompressedSizeAndIndicator(compressedData, 1,data.length); // Prepend 1 for FastPFOR
            }
        }

        private int[] prependCompressedSizeUncompressedSizeAndIndicator(int[] compressedData, int indicator, int uncompressedSize) {
            // Create a new array with space for size, indicator, and compressed data
            int[] finalData = new int[compressedData.length + 3]; // +1 for size, +1 for indicator

            // Set size as the first element
            finalData[0] = compressedData.length; // Size of the compressed data in ints!
            finalData[1] = uncompressedSize; // Size of the uncompressed data in ints!
            finalData[2] = indicator; // Set the indicator

            // Copy compressed data into the new array starting from index 2
            System.arraycopy(compressedData, 0, finalData, 3, compressedData.length);

            return finalData;
        }

        private int[] compressWithFastPFor(int[] data) {
            // Create output array with a reasonable size
           /*
           The FastPFOR algorithm requires a buffer that is large enough to store the compressed result.
           Since the actual size of the compressed data is not known beforehand, you need to over-allocate an array that can hold the worst-case size of the compressed output.

           Estimation: The size estimation data.length + (data.length / 256) + 16 is a rule of thumb for FastPFOR because:
                    data.length gives the base size.
                    data.length / 256 accounts for potential overhead since FastPFOR operates on blocks of 256 integers.
                    +16 provides extra space to account for any additional overhead or metadata that might be needed during compression.
           If the array is too small, the compression will fail or truncate the data.*/
            IntegerCODEC codec = new Composition(fastPForCompressor, this.variableByte);
            int[] compressedData = new int[data.length + (data.length / 256) + 16];
            //keeps track of how many integers have been read from the input array.
            IntWrapper inPos = new IntWrapper(0);
            //keeps track of how many integers have been written to the output array.
            IntWrapper outPos = new IntWrapper(0);

            // Call the FastPFOR compress method
            codec.compress(data, inPos, data.length, compressedData, outPos);

            // Resize the output array to the actual size
            int compressedSize = outPos.get();
            int[] finalCompressedData = new int[compressedSize];
            System.arraycopy(compressedData, 0, finalCompressedData, 0, compressedSize);

            return finalCompressedData;
        }

        public int[] uncompress(int[] compressedData, int compressionType, int uncompressedSize) {

            if (compressionType == 0) {
                // Uncompress using IntegratedIntCompressor
                return integratedCompressor.uncompress(compressedData);
            } else {
                // Uncompress using FastPFOR
                return uncompressWithFastPFor(compressedData,uncompressedSize);
            }
        }

        private int[] uncompressWithFastPFor(int[] compressedData, int uncompressedSize) {
            IntegerCODEC codec = new Composition(fastPForCompressor, this.variableByte);
            // Prepare output array with the size of the uncompressed data
            int[] uncompressedData = new int[uncompressedSize];

            // Initialize IntWrapper positions for FastPFOR
            IntWrapper inPos = new IntWrapper(0);
            IntWrapper outPos = new IntWrapper(0);

            // Call the FastPFOR uncompress method
            codec.uncompress(compressedData, inPos, compressedData.length, uncompressedData, outPos);

            return uncompressedData;
        }

        // A simple method to determine if delta compression is suitable
        private boolean isSmallDeltas(int[] data) {
            int maxDelta = 0;

            for (int i = 1; i < data.length; i++) {
                int delta = Math.abs(data[i] - data[i - 1]);
                if (delta > maxDelta) {
                    maxDelta = delta;
                }
                if (maxDelta > DELTA_THRESHOLD ) {
                    return false;
                }
            }
            return true;
        }
    }
    CombinedCompressor combinedCompressor = new CombinedCompressor();


    public EvolvingMultiGraph(String graphFile, boolean headers, String basename, long aggregationFactor) {
        super();
        this.graphFile = graphFile;
        this.headers = headers;
        this.basename = basename;
        this.aggregationFactor = aggregationFactor;
        this.combinedCompressor = new CombinedCompressor();
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
        for(Long seconds: currentNeighborsTimestamps) {
            long periodsBetween = TimestampComparerAggregator.timestampsDifference(previousNeighborTimestamp, seconds, aggregationFactor);
            periodsBetween = Fast.int2nat(periodsBetween);
            periodsBetweenList.add((int) periodsBetween);  // Store the period
            previousNeighborTimestamp = seconds;
        }

        // Compress all periods together
        int[] periodsArray = periodsBetweenList.stream().mapToInt(i -> i).toArray(); // Convert to int[]
        int[] compressedData = combinedCompressor.compress(periodsArray);  // Compress all periods at once


        // Convert int[] to byte[] for writing
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        DataOutputStream dos = new DataOutputStream(baos);
        for (int value : compressedData) {
            dos.writeInt(value);
        }
        dos.flush();
        //compressedBytes has compressed data size, compressor flag, compressed data
        byte[] compressedBytes = baos.toByteArray();

        // Write compressed bytes to OutputBitStream bit by bit = compressed data size (shows size in ints), uncompressed data size (shows size in ints), compressor flag , and then compressed data
        for (byte b : compressedBytes) {
            for (int i = 7; i >= 0; i--) {
                obs.writeBit((b >> i) & 1);
            }
        }
        //so finally, we have: compressed size, compressor flag, compressed data

        ret += compressedBytes.length * 8; // Track bits written
        return ret;
    }

    public void store() throws IOException, InterruptedException {
        ExecutorService executorService = Executors.newFixedThreadPool(2);
        executorService.execute(()-> {storeBVMultiGraph();});
        executorService.execute(()-> {storeTimestampsAndIndex();});
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

        // Read the size of the compressed data
        int compressedSizeInInts  = ibs.readInt(32); // Size in bits
        int uncompressedSizeInInts  = ibs.readInt(32); // Size in bits

        // Read the compression type (0 for IntegratedIntCompressor, 1 for FastPFOR)
        int compressionType = ibs.readInt(32);

        // Calculate the number of integers needed to store the compressed data
        int[] compressedData = new int[compressedSizeInInts];

        // Read the compressed data as ints
        for (int i = 0; i < compressedData.length; i++) {
            compressedData[i] = ibs.readInt(32); // Each int is 32 bits
        }

        // Decompress the data
        int[] decompressedData = combinedCompressor.uncompress(compressedData, compressionType, uncompressedSizeInInts);

        // Initialize index for accessing decompressed data
        int index = 0;

        if (from < 0 || to >= decompressedData.length || from > to) {
            ibs.close();
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

        ibs.close();

        return false;
    }

    public SuccessorIterator successors(int node) throws Exception {
        return new SuccessorIterator(node);
    }

    public class SuccessorIterator implements Iterator<Successor> {

        LazyIntIterator neighborsIterator;  // Iterator for neighbors
        InputBitStream ibs;                 // Bit stream for timestamps
        long previous;                       // To keep track of the previous timestamp
        int[] decompressedData;              // To store decompressed timestamps
        int index = 0;                       // Index for decompressedData

        public SuccessorIterator(int node) throws Exception {
            neighborsIterator = graph.successors(node);
            if (timestamps != null) {
                ibs = new InputBitStream(timestamps);
            } else {
                ibs = new InputBitStream(new FastMultiByteArrayInputStream(timestampsStream));
            }
            long position = efindex.getLong(node);
            ibs.position(position);

            // Read the size of the compressed data
            int compressedSizeInInts  = ibs.readInt(32); // compressed size in bits
            int uncompressedSizeInInts  = ibs.readInt(32); // uncompressed size in bits

            // Read the compression type (0 for IntegratedIntCompressor, 1 for FastPFOR)
            int compressionType = ibs.readInt(32);

            // Calculate the number of integers needed to store the compressed data
            int[] compressedData = new int[compressedSizeInInts];

            // Read the compressed data as ints
            for (int i = 0; i < compressedData.length; i++) {
                compressedData[i] = ibs.readInt(32); // Each int is 32 bits
            }

            // Decompress the data
            this.decompressedData = combinedCompressor.uncompress(compressedData, compressionType, uncompressedSizeInInts);
            previous = minTimestamp; // Initialize previous timestamp
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
