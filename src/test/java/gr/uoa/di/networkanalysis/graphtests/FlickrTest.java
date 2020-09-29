package gr.uoa.di.networkanalysis.graphtests;


import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

import org.junit.Assert;
import org.junit.Test;

import gr.uoa.di.networkanalysis.EvolvingMultiGraph;
import gr.uoa.di.networkanalysis.EvolvingMultiGraph.SuccessorIterator;
import gr.uoa.di.networkanalysis.TimestampComparer;
import gr.uoa.di.networkanalysis.TimestampComparerAggregator;
import gr.uoa.di.networkanalysis.Successor;

public class FlickrTest {

	private static final int factor = 2*24*60*60;
	
	private static TimestampComparer ic = new TimestampComparerAggregator(factor);
	
	@Test
	public void testAll() throws Exception {
		testStore();
		testLoadAndSuccesors();
	}
	
//	@Test
	public void testStore() throws Exception {
		EvolvingMultiGraph emg = new EvolvingMultiGraph(
				"out.flickr-growth-sorted.gz",
				true,
				4,
				"flickr",
				ic
		);

		emg.store();
	}
	
//	@Test
	public void testLoadAndSuccesors() throws Exception {
		
		EvolvingMultiGraph emg = new EvolvingMultiGraph(
				"out.flickr-growth-sorted.gz",
				true,
				4,
				"flickr",
				ic
		);
		
		emg.load();
		
		FileInputStream fileStream = new FileInputStream("out.flickr-growth-sorted.gz");
        GZIPInputStream gzipStream = new GZIPInputStream(fileStream);
        InputStreamReader decoder = new InputStreamReader(gzipStream, "UTF-8");
        BufferedReader buffered = new BufferedReader(decoder);
		
        int current = 1;
        String line = buffered.readLine(); // Get rid of headers
        ArrayList<Successor> list = new ArrayList<Successor>();
        
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
            	SuccessorIterator it = emg.successors(current);
            	int i = 0;
        		while(true) {
        			try {
        				Successor s = it.next();
        				Assert.assertEquals((double) s.getTimestamp(), (double) list.get(i).getTimestamp(), factor);
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
        
        SuccessorIterator it = emg.successors(current);
    	int i = 0;
		while(true) {
			try {
				Successor s = it.next();
				Assert.assertEquals((double) s.getTimestamp(), (double) list.get(i).getTimestamp(), factor);
				i++;
			}
			catch(NoSuchElementException e) {
				break;
			}
		}
        
		buffered.close();
	}
}
