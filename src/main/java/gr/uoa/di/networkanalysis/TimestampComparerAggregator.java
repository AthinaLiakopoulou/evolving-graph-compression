package gr.uoa.di.networkanalysis;

public class TimestampComparerAggregator {

	
	public static double timestampsDifference(long t1, long t2, double factor) {
		// TODO: Check for a better way to do this. For min timestamps as well!!!
		return (t2 / factor) - (t1 / factor);
	}

	public static long reverse(long previous, long difference, double factor) {
		return previous + Math.round(factor * difference);
	}

	public static long aggregateMinTimestamp(long min, double factor) {
		return Math.round(Math.floor(min / factor) * factor);
	}

}
