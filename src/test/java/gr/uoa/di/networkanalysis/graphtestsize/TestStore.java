package gr.uoa.di.networkanalysis.graphtestsize;


import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.*;
import java.util.zip.GZIPInputStream;

import com.google.common.io.Resources;
import me.lemire.integercompression.differential.IntegratedIntCompressor;
import org.junit.Assert;
import org.junit.Test;
import gr.uoa.di.networkanalysis.EvolvingMultiGraph;
import gr.uoa.di.networkanalysis.EvolvingMultiGraph.SuccessorIterator;
import gr.uoa.di.networkanalysis.Successor;

public class TestStore {

	// Flickr
//	private static final String graphFile =  "out.flickr-growth.sorted.gz";
//	private static final String basename =  "flickr";
//	private static final boolean headers = true;
//	private static final int k = 2;
//	private static int aggregation = 24*60*60;

	// Wiki
//	private static final String graphFile =  "out.edit-enwiki.sorted.gz";
//	private static final String basename =  "wiki";
//	private static final boolean headers = true;
//	private static final int k = 2;
//	private static int aggregation = 60*60;

	// Yahoo
//	private static final String graphFile =  "yahoo-G5-sorted.tsv.gz";
//	private static final String basename =  "yahoo";
//	private static final boolean headers = false;
//	private static final int k = 2;
//	private static int aggregation = 15*60;

	// cbtComm
	private static final String graphFile =  "cbtComm-sorted.txt.gz";
	private static String sampleFile = "cbtComm-sample.txt";
	private static final String basename =  "cbtComm";
	private static final boolean headers = false;
	private static final int k = 2;
	private static double aggregation =  0.00001;

	// cbtPow
//	private static final String graphFile =  "cbtPow-sorted.txt.gz";
//	private static final String basename =  "cbtPow";
//	private static final boolean headers = false;
//	private static final int k = 2;
//	private static int aggregation = 1;

	
	@Test
	public void testStore() throws Exception {

		ClassLoader classLoader = getClass().getClassLoader();
		String graphFileResourcePath = classLoader.getResource(graphFile).getPath();

		EvolvingMultiGraph emg = new EvolvingMultiGraph(
				graphFileResourcePath,
				headers,
				basename,
				aggregation
		);
		
		long t1 = System.nanoTime();
		emg.store();
		long t2 = System.nanoTime();
		System.out.println("Compression took: " + (t2-t1) + " nanoseconds");
	}
	
	@Test
	public void assertEqualsFromOriginalFile() throws Exception {

		String graphFileResourcePath = Resources.getResource(graphFile).getPath();

		EvolvingMultiGraph emg = new EvolvingMultiGraph(
				graphFileResourcePath,
				headers,
				basename,
				aggregation
		);
		
		emg.load();

		FileInputStream fileStream = new FileInputStream(graphFileResourcePath);
        GZIPInputStream gzipStream = new GZIPInputStream(fileStream);
        InputStreamReader decoder = new InputStreamReader(gzipStream, "UTF-8");
        BufferedReader buffered = new BufferedReader(decoder);

        int current = 0;
        String line;
        if(headers) {
        	buffered.readLine(); // Get rid of headers
		}
        ArrayList<Successor> list = new ArrayList<>();

        while ((line = buffered.readLine()) != null) {
        	String[] tokens = line.split("\\s+");
        	int node = Integer.parseInt(tokens[0]);
            int neighbor = Integer.parseInt(tokens[1]);
            long timestamp = Long.parseLong(tokens[3]);

            if(node == current) {
            	list.add(new Successor(neighbor, timestamp));
            }
            else {
            	// Check the list so far
				list.sort(Comparator.comparing(Successor::getTimestamp));
            	SuccessorIterator it = emg.successors(current);
            	int i = 0;
        		while(true) {
        			try {
        				Successor s = it.next();
        				Assert.assertEquals(s.getTimestamp(), list.get(i).getTimestamp(), aggregation);
        				i++;
        			}
        			catch(NoSuchElementException e) {
        				break;
        			}
        		}

            	list = new ArrayList<Successor>();
            	list.add(new Successor(neighbor, timestamp));
            	current = node;
            }
        }
		list.sort(Comparator.comparing(Successor::getTimestamp));
        SuccessorIterator it = emg.successors(current);
    	int i = 0;
		while(true) {
			try {
				Successor s = it.next();
				Assert.assertEquals(s.getTimestamp(), list.get(i).getTimestamp(), aggregation);
				i++;
			}
			catch(NoSuchElementException e) {
				break;
			}

		}
	}

	@Test
	public void mytest() {

		IntegratedIntCompressor iic = new IntegratedIntCompressor();
		int[] data = new int[2342351];
		for(int k = 0; k < data.length; ++k)
			data[k] = k;
		//System.out.println("orig: "+Arrays.toString(data));
		System.out.println("Compressing "+data.length+" integers using friendly interface");
		int[] compressed = iic.compress(data);
		//System.out.println("compressed: "+Arrays.toString(compressed));
		int[] recov = iic.uncompress(compressed);
		//System.out.println("uncompress: "+Arrays.toString(recov));
		System.out.println("compressed from "+data.length*4/1024+"KB to "+compressed.length*4/1024+"KB");
		if(!Arrays.equals(recov,data)) throw new RuntimeException("bug");
	}
}
