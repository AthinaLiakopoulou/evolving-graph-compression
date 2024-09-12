package gr.uoa.di.networkanalysis;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import java.util.NoSuchElementException;
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

        // Returns the number of bits appended to the file
        long ret = 0;
        long previousNeighborTimestamp = minTimestamp;

        for(Long seconds: currentNeighborsTimestamps) {
            double periodsBetween = TimestampComparerAggregator.timestampsDifference(previousNeighborTimestamp, seconds, aggregationFactor);
            periodsBetween = Fast.int2nat(Math.round(periodsBetween));
            previousNeighborTimestamp = seconds;

            // Compress the data
            int[] compressedData = compressor.compress(new int[]{(int) periodsBetween});

            // Convert int[] to byte[]
            ByteArrayOutputStream baos = new ByteArrayOutputStream();
            DataOutputStream dos = new DataOutputStream(baos);
            for (int value : compressedData) {
                dos.writeInt(value);
            }
            dos.flush();
            byte[] compressedBytes = baos.toByteArray();

            // Write the size of the compressed data (in bits)
            int compressedSize = compressedBytes.length * 8;
            obs.writeInt(compressedSize,32); // Write the size

            // Write compressed bytes to OutputBitStream bit by bit
            for (byte b : compressedBytes) {
                for (int i = 7; i >= 0; i--) {
                    obs.writeBit((b >> i) & 1);
                }
            }

            ret += 32 + compressedBytes.length * 8; // Number of bits written
        }
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
        int compressedSize = ibs.readInt(32); // Size in bits
        int[] compressedData = new int[compressedSize / 32]; // Each int is 32 bits

        // Read the compressed data as ints
        for (int i = 0; i < compressedData.length; i++) {
            compressedData[i] = ibs.readInt(32); // Each int is 32 bits
        }

        // Decompress the data
        int[] decompressedData = compressor.uncompress(compressedData);

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

            int compressedSize = ibs.readInt(32); // Read the size of the compressed data
            int[] compressedData = new int[compressedSize / 32]; // Initialize array to hold compressed data

            for (int i = 0; i < compressedData.length; i++) {
                compressedData[i] = ibs.readInt(32); // Read the compressed ints
            }

            // Decompress the data
            decompressedData = compressor.uncompress(compressedData);
            
            previous = minTimestamp;
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
