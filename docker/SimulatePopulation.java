import java.util.Arrays;
import java.util.Random;
import java.util.PriorityQueue;
import java.io.*;
import java.text.NumberFormat;


/**
 * Creates a set of chr1 instances using a model learnt from gnomAD-SV.
 *
 * Remark: the program does not currently represent a DNA character in just two
 * bits, so it wastes a lot of memory. This is done for simplicity and should be
 * improved in the future.
 */
public class SimulatePopulation {
	/**
	 * Reference chromosome
	 */
	private static StringBuilder reference;
	public static int referenceLength;
	
	/**
	 * $positions[i]$ contains the set of all the distinct positions of the 
	 * reference that belong to a repeat interval of type $i$ (see 
	 * $BuildModel.N_INTERVAL_TYPES$ for details). The same position of the 
	 * reference can belong to intervals of different type. The last row
	 * contains all the positions that do not belong to any repeat interval.
	 *
	 * Remark: the program samples from these positions uniformly at random, so
	 * hotspots of variation that are not explained by repeat content alone 
	 * cannot be captured.
	 */
	private static int[][] positions;
	private static int[] positions_last;
	
	/**
	 * Set of realistic insertion sequences. The program never inserts random
	 * sequences.
	 */
	private static String[] insertions;
	private static int lastInsertion;
	
	/**
	 * The first element of $insertions$ with length in bin $i$, or -1 if bin
	 * $i$ is empty.
	 */
	private static int[] firstInsertionWithLength;
	
	/**
	 * Variants
	 */
	private static Variant[] variants;
	private static int lastVariant, nVariants_frequent, nVariants_rare;
	
	/**
	 * For every chromosome (rows), the sorted list of variants it contains.
	 */
	private static int N_CHROMOSOMES;
	private static int[][] chromosome2variants;
	private static int[] chromosome2variants_last;
	
	/**
	 * Inverted version of $chromosome2variants$.
	 */
	private static int[][] variant2chromosomes;
	private static int[] variant2chromosomes_last;
	
	/**
	 * Used throughout the code
	 */
	private static Random random = new Random();
	private static Pair[] tmpPairs;
	

	public static void main(String[] args) throws IOException {
		final String MODEL_DIR = args[0];
		final String POPULATION_NAME = args[1];
		final String REFERENCE_FILE = args[2];
		final String INSERTION_STRINGS_FILE = args[3];
		final int N_INSERTION_STRINGS = Integer.parseInt(args[4]);
		final double ALLELE_FREQUENCY_THRESHOLD = Double.parseDouble(args[5]);
		final int N_VARIANTS_FREQUENT = Integer.parseInt(args[6]);
		final int N_INDIVIDUALS = Integer.parseInt(args[7]);
		final String OUTPUT_DIR = args[8];
		
		final int MAX_DISTANCE_TO_REPEAT_INTERVAL = 100;  // Arbitrary
		StringBuilder buffer = new StringBuilder(10000);  // Arbitrary
		
		System.err.println("Loading reference...");
		loadReference(REFERENCE_FILE);
		System.err.println("Loading model...");
		BuildModel.deserializeStatistics(MODEL_DIR);
		BuildModel.deserializeRepeatIntervals(MODEL_DIR+"/repeatIntervals.txt");
		BuildModel.deserializeSegdups(MODEL_DIR+"/segdups.txt");
		System.err.println("Loading repetitive positions...");
		loadPositions(referenceLength,BuildModel.MAX_DISTANCE_SR);
		System.err.println("Sampling variants...");
		sampleVariants(ALLELE_FREQUENCY_THRESHOLD,N_VARIANTS_FREQUENT,Arrays.binarySearch(BuildModel.AF_NAMES,POPULATION_NAME));
		System.err.println("DONE  Total variants="+(lastVariant+1)+". Frequent="+nVariants_frequent+" ("+((100.0*nVariants_frequent)/(lastVariant+1))+"%). Rare="+nVariants_rare+" ("+((100.0*nVariants_rare)/(lastVariant+1))+"%).");
		System.err.println("Distributing variants to the population...");
		buildPopulation(N_INDIVIDUALS,ALLELE_FREQUENCY_THRESHOLD,OUTPUT_DIR);
		buildVariant2chromosomes();
		System.err.println("Loading insertion strings...");
		loadInsertions(INSERTION_STRINGS_FILE,N_INSERTION_STRINGS,BuildModel.SV_LENGTH_QUANTUM);
		System.err.println("Assigning insertion strings to variants...");
		assignInsertionStrings(buffer);
		System.err.println("Annotating variants with repeats...");
		annotateVariantsWithRepeats(MAX_DISTANCE_TO_REPEAT_INTERVAL);
		buildVCF(OUTPUT_DIR);
		System.err.println("Saving population to disk...");
		serialize(OUTPUT_DIR+"/chromosome2variants.txt",OUTPUT_DIR+"/variants.txt");
	}
	
	
	public static final void loadReference(String path) throws IOException {
		String str;
		BufferedReader br;
		
		reference = new StringBuilder(250000000);  // Arbitrary
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();  // Skipping header
		str=br.readLine();
		while (str!=null) {
			reference.append(str.trim());
			str=br.readLine();
		}
		br.close();
		referenceLength=reference.length();
		System.err.println("Length of the reference: "+referenceLength);
	}
	
	
	/**
	 * Remark: positions whose character is not in the alphabet are discarded.
	 *
	 * @param maxDistance positions at distance $<=maxDistance$ from an interval
	 * of type X are assigned to type X.
	 */
	private static final void loadPositions(int chromosomeLength, int maxDistance) {
		char c;
		int i, j, k, p;
		int type, last1, last2;
		int[] tmpArray1, tmpArray2;
		
		// Collecting positions from repeat intervals and segdups
		positions = new int[BuildModel.N_INTERVAL_TYPES][chromosomeLength];
		positions_last = new int[BuildModel.N_INTERVAL_TYPES];
		Arrays.fill(positions_last,-1);
		for (i=0; i<=BuildModel.lastRepeatInterval; i++) {
			type=BuildModel.repeatIntervals[i][0];
			for (j=BuildModel.repeatIntervals[i][1]-maxDistance; j<=BuildModel.repeatIntervals[i][2]+maxDistance; j++) positions[type][++positions_last[type]]=j;
		}
		for (i=0; i<=BuildModel.lastSegdup; i++) {
			for (j=BuildModel.segdups[i][0]-maxDistance; j<=BuildModel.segdups[i][1]+maxDistance; j++) positions[3][++positions_last[3]]=j;
		}
		for (i=0; i<positions.length; i++) {
			Arrays.sort(positions[i],0,positions_last[i]+1);
			j=0;
			for (k=1; k<=positions_last[i]; k++) {
				if (positions[i][k]!=positions[i][k-1]) positions[i][++j]=positions[i][k];
			}
			positions_last[i]=j;
		}
		
		// Adding to BOTH the intersection of REP and SAT
		tmpArray1 = new int[chromosomeLength];
		last1=setIntersection(positions[0],0,positions_last[0],positions[1],0,positions_last[1],tmpArray1,0);
		tmpArray2 = new int[chromosomeLength];
		last2=setUnion(tmpArray1,last1,positions[2],positions_last[2],tmpArray2);
		System.arraycopy(tmpArray2,0,positions[2],0,last2+1);
		positions_last[2]=last2;
		
		// Subtracting BOTH from REP and SAT
		last1=setMinus(positions[0],positions_last[0],positions[2],positions_last[2],tmpArray1);
		System.arraycopy(tmpArray1,0,positions[0],0,last1+1);
		positions_last[0]=last1;
		last1=setMinus(positions[1],positions_last[1],positions[2],positions_last[2],tmpArray1);
		System.arraycopy(tmpArray1,0,positions[1],0,last1+1);
		positions_last[1]=last1;
		
		// Subtracting REP,SAT,BOTH from SEGDUP.
		last1=setMinus(positions[3],positions_last[3],positions[0],positions_last[0],tmpArray1);
		last2=setMinus(tmpArray1,last1,positions[1],positions_last[1],tmpArray2);
		last1=setMinus(tmpArray2,last2,positions[2],positions_last[2],tmpArray1);
		System.arraycopy(tmpArray1,0,positions[3],0,last1+1);
		positions_last[3]=last1;
		
		// Building nonrepetitive positions
		last1=setUnion(positions[0],positions_last[0],positions[1],positions_last[1],tmpArray1);
		last2=setUnion(tmpArray1,last1,positions[2],positions_last[2],tmpArray2);
		last1=setUnion(tmpArray2,last2,positions[3],positions_last[3],tmpArray1);
		positions_last[4]=setMinus(chromosomeLength,tmpArray1,last1,positions[4]);
		
		// Filtering positions by character
		for (i=0; i<positions.length; i++) {
			k=-1; last1=positions_last[i];
			for (j=0; j<=last1; j++) {
				p=positions[i][j];
				if (p>=referenceLength) continue;
				c=reference.charAt(p);
				if (c!='A' && c!='C' && c!='G' && c!='T' && c!='a' && c!='c' && c!='g' && c!='t') continue;
				positions[i][++k]=p;
			}
			positions_last[i]=k;
		}
	}
	
	
	/**
	 * Stores in a prefix of $to$ the union of the sorted sets $from1[0..last1]$
	 * and $from2[0..last2]$.
	 *
	 * @return the last element of $to$.
	 */
	private static final int setUnion(int[] from1, int last1, int[] from2, int last2, int[] to) {
		int i1, i2, j;
	
		i1=0; i2=0; j=0;
		while (i1<=last1 && i2<=last2) {
			if (from1[i1]<from2[i2]) to[j++]=from1[i1++];
			else if (from1[i1]>from2[i2]) to[j++]=from2[i2++];
			else { 
				to[j++]=from1[i1];
				i1++; i2++;
			}
		}
		while (i1<=last1) to[j++]=from1[i1++];
		while (i2<=last2) to[j++]=from2[i2++];
		return j-1;
	}
	
	
	/**
	 * Stores in a prefix of $to$ the subtraction of the sorted set $from2[0..
	 * last2]$ from the sorted set $from1[0..last1]$.
	 *
	 * @return the last element of $to$.
	 */
	public static final int setMinus(int[] from1, int last1, int[] from2, int last2, int[] to) {
		int i1, i2, j;
	
		i1=0; i2=0; j=0;
		while (i1<=last1 && i2<=last2) {
			if (from1[i1]<from2[i2]) to[j++]=from1[i1++];
			else if (from1[i1]>from2[i2]) i2++;
			else { i1++; i2++; }
		}
		while (i1<=last1) to[j++]=from1[i1++];
		return j-1;
	}
	
	
	/**
	 * Stores in a prefix of $to$ the subtraction of $from2[0..last2]$ from 
	 * $[0..chromosomeLength-1]$.
	 *
	 * @return the last element of $to$.
	 */
	private static final int setMinus(int chromosomeLength, int[] from2, int last2, int[] to) {
		int i1, i2, j;
	
		i1=0; i2=0; j=0;
		while (i1<chromosomeLength && i2<=last2) {
			if (i1<from2[i2]) to[j++]=i1++;
			else if (i1>from2[i2]) i2++;
			else { i1++; i2++; }
		}
		while (i1<chromosomeLength) to[j++]=i1++;
		return j-1;
	}
	
	
	/**
	 * Stores in $y[fromY..]$ the intersection of the sorted sets $x1[from1..last1]$ and 
	 * $x2[from2..last2]$.
	 *
	 * @return the last element of $y$.
	 */
	private static final int setIntersection(int[] x1, int from1, int last1, int[] x2, int from2, int last2, int[] y, int fromY) {
		int i1, i2, j;
	
		i1=from1; i2=from2; j=fromY;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]<x2[i2]) i1++;
			else if (x1[i1]>x2[i2]) i2++;
			else {
				y[j++]=x1[i1];
				i1++; i2++;
			}
		}
		return j-1;
	}
	
	
	/**
	 * Remark: some bins of $firstInsertionWithLength$ might be empty.
	 *
	 * @param path assumed to contain no duplicates and to be sorted by
	 * increasing length.
	 */
	private static final void loadInsertions(String path, int nInsertions, int quantum) throws IOException {
		int i, j;
		int nBins;
		String str;
		BufferedReader br;
		
		// Loading sequences
		insertions = new String[nInsertions];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); lastInsertion=-1;
		while (str!=null) {
			insertions[++lastInsertion]=str.trim();
			str=br.readLine();
		}
		br.close();
		
		// Building index
		nBins=1+insertions[lastInsertion].length()/quantum;
		firstInsertionWithLength = new int[nBins];
		Arrays.fill(firstInsertionWithLength,-1);
		for (i=0; i<=lastInsertion; i++) {
			j=insertions[i].length()/quantum;
			if (firstInsertionWithLength[j]==-1) firstInsertionWithLength[j]=i;
		}
	}
	
	
	/**
	 * Stores in global variable $variants$ at least $nFrequent$ variants with 
	 * frequency $>= alleleFrequencyThreshold$, where frequency is sampled from 
	 * the $populationID$-th population in $BuildModel.AF_FILE_NAMES$. Variants 
	 * are sorted by their first position in the genome.
	 *
	 * Remark: the procedure does not explicitly model satellite expansion/
	 * contraction, i.e. a set of variants with the same position but with
	 * different lengths. This is left to the future, assuming that a dataset
	 * from which to learn distributions exists.
	 *
	 * Remark: $positions$ is deallocated at the end of the procedure.
	 */
	private static final void sampleVariants(double alleleFrequencyThreshold, int nFrequent, int populationID) {
		Variant variant;
		
		variants = new Variant[nFrequent<<1];  // Arbitrary
		nVariants_frequent=0; nVariants_rare=0; lastVariant=-1;
		while (nVariants_frequent<nFrequent) {
			lastVariant++;
			if (lastVariant==variants.length) {
				Variant[] newArray = new Variant[variants.length<<1];
				System.arraycopy(variants,0,newArray,0,variants.length);
				variants=newArray;
			}
			do { variant=sampleVariant(populationID); }
			while (variant.svType!=2 && variant.position+variant.length>referenceLength);
			variants[lastVariant]=variant;
			if (variants[lastVariant].alleleFrequency>=alleleFrequencyThreshold) nVariants_frequent++;
			else nVariants_rare++;
		}
		Arrays.sort(variants,0,lastVariant+1);
		positions=null; positions_last=null;
		System.gc();
	}
	
	
	/**
	 * Remark: after deciding the length (respectively, frequency) bin according
	 * to the distribution, the variant is assigned a length (resp. frequency)
	 * inside the bin uniformly at random.
	 *
	 * @param populationID index in $BuildModel.AF_FILE_NAMES$.
	 */
	private static final Variant sampleVariant(int populationID) {
		Variant out = new Variant();
		out.svType=Arrays.binarySearch(BuildModel.statistics_svType,random.nextDouble());
		if (out.svType<0) out.svType=-1-out.svType;
		out.intervalType=Arrays.binarySearch(BuildModel.statistics_intervalType[out.svType],random.nextDouble());
		if (out.intervalType<0) out.intervalType=-1-out.intervalType;
		int p = Arrays.binarySearch(BuildModel.statistics_length[out.svType][out.intervalType],random.nextDouble());
		if (p<0) p=-1-p;
		out.length=BuildModel.MIN_SV_LENGTH+BuildModel.SV_LENGTH_QUANTUM*p+random.nextInt(BuildModel.SV_LENGTH_QUANTUM);
		final double min = BuildModel.statistics_alleleFrequency_minmax[populationID][out.svType][0];
		final double quantum = (BuildModel.statistics_alleleFrequency_minmax[populationID][out.svType][1]-min)/(BuildModel.ALLELE_FREQUENCY_NBINS-1);
		p=Arrays.binarySearch(BuildModel.statistics_alleleFrequency[populationID][out.svType],0,BuildModel.ALLELE_FREQUENCY_NBINS,random.nextDouble());
		if (p<0) p=-1-p;
		out.alleleFrequency=min+quantum*p+random.nextDouble()*quantum;
		out.position=positions[out.intervalType][random.nextInt(positions_last[out.intervalType]+1)];
		return out;
	}
	
	
	/**
	 * Sets values $intervalType_start, intervalType_end, intervalSurface$ of 
	 * every variant.
	 *
	 * Remark: the procedure assumes that $variants$ is sorted by $position$.
	 *
	 * @param maxDistance only repeat intervals at distance $<=maxDistance$ from
	 * a variant endpoint are used to annotate the breakpoint.
	 */
	private static final void annotateVariantsWithRepeats(int maxDistance) {
		int i;
		int last;
		int[] tmpArray1 = new int[2];
		int[] tmpArray2 = new int[2];
		int[] tmpArray3 = new int[2];
		
		tmpArray1[0]=0; tmpArray1[1]=0; tmpArray2[0]=0; tmpArray2[1]=0;
		for (i=0; i<=lastVariant; i++) {
			variants[i].intervalType_start=getIntervalType(variants[i].position,maxDistance,tmpArray1,true);
			if (variants[i].svType==2) {
				variants[i].intervalType_end=variants[i].intervalType_start;
				variants[i].intervalSurface=0;
			}
			else {
				last=variants[i].position+variants[i].length-1;
				variants[i].intervalType_end=getIntervalType(last,maxDistance,tmpArray1,false);
				variants[i].intervalSurface=getIntervalSurface(variants[i].position,last,tmpArray2,tmpArray3);
			}
		}
	}
	
	
	/**
	 * @param cursors positions in $repeatIntervals$ (0) and $segdups$ (1);
	 * @param updateCursors TRUE=updates $cursors$ with the current values after
	 * the procedure completes;
	 * @return the union of the types of all repeat intervals that are at most
	 * $maxDistance$ away from $position$.
	 */
	private static final int getIntervalType(int position, int maxDistance, int[] cursors, boolean updateCursors) {
		int i1, i2, j;
		int intervalType;
		
		i1=cursors[0]; i2=cursors[1];
		while (i1<=BuildModel.lastRepeatInterval && BuildModel.repeatIntervals[i1][2]+maxDistance<position) {
			i1++;
			continue;
		}
		if (i1>BuildModel.lastRepeatInterval || BuildModel.repeatIntervals[i1][1]>position+maxDistance) intervalType=4;
		else {
			intervalType=BuildModel.repeatIntervals[i1][0];
			for (j=i1-1; j>=0; j--) {
				if (BuildModel.repeatIntervals[j][2]+maxDistance<position) break;
				intervalType=BuildModel.ADD_TYPE[intervalType][BuildModel.repeatIntervals[j][0]];
			}
			for (j=i1+1; j<=BuildModel.lastRepeatInterval; j++) {
				if (BuildModel.repeatIntervals[j][1]>position+maxDistance) break;
				intervalType=BuildModel.ADD_TYPE[intervalType][BuildModel.repeatIntervals[j][0]];
			}
		}
		if (intervalType==4) {
			while (i2<=BuildModel.lastSegdup && BuildModel.segdups[i2][1]+maxDistance<position) {
				i2++;
				continue;
			}
			if (i2>BuildModel.lastSegdup || BuildModel.segdups[i2][0]>position+maxDistance) intervalType=4;
			else intervalType=3;
		}
		if (updateCursors) { cursors[0]=i1; cursors[1]=i2; }
		return intervalType;
	}
	
	
	/**
	 * Remark: the procedure uses global variable $tmpPairs$.
	 *
	 * @param cursors position in $repeatIntervals$ (0) and $segdups$ (1);
	 * @param tmpArray temporary space of size at least 2;
	 * @return the total number of basepairs inside $[first..last]$ that are 
	 * covered by repeats of any type.
	 */
	private static final int getIntervalSurface(int first, int last, int[] cursors, int[] tmpArray) {
		int i, j;
		int lastPair, sum;
		
		// Collecting intervals
		if (tmpPairs==null) {
			tmpPairs = new Pair[100];  // Arbitrary
			for (i=0; i<tmpPairs.length; i++) tmpPairs[i] = new Pair();
		}
		tmpArray[0]=cursors[0]; tmpArray[1]=-1;
		getIntersectingIntervals(first,last,BuildModel.repeatIntervals,BuildModel.lastRepeatInterval,1,2,tmpArray);
		cursors[0]=tmpArray[0]; 
		tmpArray[0]=cursors[1];
		getIntersectingIntervals(first,last,BuildModel.segdups,BuildModel.lastSegdup,0,1,tmpArray);
		cursors[1]=tmpArray[0];
		lastPair=tmpArray[1];
		if (lastPair!=0) Arrays.sort(tmpPairs,0,lastPair+1);
		
		// Merging intervals
		j=0; sum=0;
		for (i=1; i<=lastPair; i++) {
			if (tmpPairs[i].first<=tmpPairs[j].last) tmpPairs[j].last=Math.max(tmpPairs[j].last,tmpPairs[i].last);
			else {
				sum+=tmpPairs[j].last-tmpPairs[j].first+1;
				j++;
				tmpPairs[j].first=tmpPairs[i].first;
				tmpPairs[j].last=tmpPairs[i].last;
			}
		}
		sum+=tmpPairs[j].last-tmpPairs[j].first+1;
		return sum;
	}
	
	
	/**
	 * Appends to global variable $tmpPairs$ all the intersections between 
	 * $[first..last]$ and an interval in $matrix$.
	 *
	 * @param column_* the column of $matrix$ that stores the first/last 
	 * position of the corresponding interval;
	 * @param cursors 0=position in $matrix$; 1=position in $tmpPairs$; these
	 * are updated when the procedure completes.
	 */
	private static final void getIntersectingIntervals(int first, int last, int[][] matrix, int lastRow, int column_first, int column_last, int[] cursors) {
		int i, j, k;
		
		i=cursors[0]; j=cursors[1];
		while (i<=lastRow && matrix[i][column_last]<first) {
			i++;
			continue;
		}
		if (i<=lastRow && matrix[i][column_first]<=last) {
			j++;
			if (j==tmpPairs.length) {
				Pair[] newArray = new Pair[tmpPairs.length<<1];
				System.arraycopy(tmpPairs,0,newArray,0,tmpPairs.length);
				for (k=tmpPairs.length; k<newArray.length; k++) newArray[k] = new Pair();
				tmpPairs=newArray;
			}
			tmpPairs[j].first=Math.max(first,matrix[i][column_first]);
			tmpPairs[j].last=Math.min(last,matrix[i][column_last]);
			for (k=i-1; k>=0; k--) {
				if (matrix[k][column_last]<first) break;
				j++;
				if (j==tmpPairs.length) {
					Pair[] newArray = new Pair[tmpPairs.length<<1];
					System.arraycopy(tmpPairs,0,newArray,0,tmpPairs.length);
					for (k=tmpPairs.length; k<newArray.length; k++) newArray[k] = new Pair();
					tmpPairs=newArray;
				}
				tmpPairs[j].first=Math.max(first,matrix[k][column_first]);
				tmpPairs[j].last=Math.min(last,matrix[k][column_last]);
			}
			for (k=i+1; k<=lastRow; k++) {
				if (matrix[k][column_first]>last) break;
				j++;
				if (j==tmpPairs.length) {
					Pair[] newArray = new Pair[tmpPairs.length<<1];
					System.arraycopy(tmpPairs,0,newArray,0,tmpPairs.length);
					for (k=tmpPairs.length; k<newArray.length; k++) newArray[k] = new Pair();
					tmpPairs=newArray;
				}
				tmpPairs[j].first=Math.max(first,matrix[k][column_first]);
				tmpPairs[j].last=Math.min(last,matrix[k][column_last]);
			}
		}
		cursors[0]=i; cursors[1]=j;
	}
	
	
	private static class Pair implements Comparable {
		public int first, last;
		
		public int compareTo(Object other) {
			Pair otherPair = (Pair)other;
			if (first<otherPair.first) return -1;
			else if (first>otherPair.first) return 1;
			if (last<otherPair.last) return -1;
			else if (last>otherPair.last) return 1;
			return 0;
		}
	}
	
	
	/**
	 * Builds global variable $chromosome2variants$ by assigning the sampled 
	 * variants to chromosomes at random, and by fixing collisions inside a 
	 * chromosome using a simple heuristic for maximum independent set. Variants
	 * that do not appear in any chromosome after the procedure completes are 
	 * marked in their $isActive$ field, and counts $nVariants_{frequent,rare}$ 
	 * are updated.
	 *
	 * Remark: variants can overlap in the population (but of course not on a
	 * single chromosome).
	 *
	 * Remark: there is no notion of "individual" in this procedure.
	 *
	 * Remark: $variants$ is assumed to be sorted by position in the reference.
	 */
	private static final void buildPopulation(int nIndividuals, double alleleFrequencyThreshold, String outputDir) throws IOException {
		N_CHROMOSOMES=nIndividuals<<1;
		final int COLLISION_DISTANCE = 5;  // Arbitrary, bps.
		final int HISTOGRAM_MAX_LENGTH = 4000;  // Arbitrary
		int i, j, k, h;
		int nActive, nEdges, lastPosition, nChromosomesWithCollision;
		int chromosome2variantsLast, variantJ, variantK, neighbor;
		QueueElement element;
		int[] lastNeighbor, histogram;
		QueueElement[] queueElements;
		int[][] neighbors;
		PriorityQueue<QueueElement> queue;
		BufferedWriter bw;
		
		System.err.println("Assigning each variant to chromosomes in isolation...");
		chromosome2variants = new int[N_CHROMOSOMES][1000];  // Arbitrary
		chromosome2variants_last = new int[N_CHROMOSOMES];
		Arrays.fill(chromosome2variants_last,-1);
		for (j=0; j<=lastVariant; j++) {
			if (j%10000==0) System.err.println("Processed "+j+" variants");
			for (i=0; i<N_CHROMOSOMES; i++) {
				if (random.nextDouble()<=variants[j].alleleFrequency) {
					chromosome2variants_last[i]++;
					if (chromosome2variants_last[i]==chromosome2variants[i].length) {
						int[] newArray = new int[chromosome2variants[i].length<<1];
						System.arraycopy(chromosome2variants[i],0,newArray,0,chromosome2variants[i].length);
						chromosome2variants[i]=newArray;
					}
					chromosome2variants[i][chromosome2variants_last[i]]=j;
				}
			}
		}
		
		System.err.println("Fixing collisions inside each chromosome...");
		neighbors = new int[lastVariant+1][2];  // Arbitrary
		lastNeighbor = new int[lastVariant+1];
		queue = new PriorityQueue<QueueElement>(lastVariant+1);
		queueElements = new QueueElement[lastVariant+1];
		for (i=0; i<=lastVariant; i++) queueElements[i] = new QueueElement(i);
		nChromosomesWithCollision=0;
		for (i=0; i<N_CHROMOSOMES; i++) {
			if (i%1000==0) System.err.println("Processed "+i+" chromosomes");
			Arrays.fill(lastNeighbor,-1);
			nEdges=0; chromosome2variantsLast=chromosome2variants_last[i];
			for (j=0; j<=chromosome2variantsLast; j++) {
				variantJ=chromosome2variants[i][j];
				lastPosition=variants[variantJ].position+(variants[variantJ].svType==2?0:variants[variantJ].length-1);
				for (k=j+1; k<=chromosome2variantsLast; k++) {
					variantK=chromosome2variants[i][k];
					if (variants[variantK].position>lastPosition+COLLISION_DISTANCE) break;
					if (variants[variantK].collidesWith(variants[variantJ],COLLISION_DISTANCE)) {
						addEdge(variantJ,variantK,neighbors,lastNeighbor);
						nEdges++;
					}
				}
			}
			if (nEdges==0) continue;
			nChromosomesWithCollision++;
			queue.clear();
			for (j=0; j<=chromosome2variantsLast; j++) {
				variantJ=chromosome2variants[i][j];
				queueElements[variantJ].degree=lastNeighbor[variantJ]+1;
				queueElements[variantJ].removed=false;
				queue.add(queueElements[variantJ]);
			}
			while (!queue.isEmpty()) {
				element=queue.poll();
				if (element.removed) continue;
				variantJ=element.variantID;
				for (k=0; k<=lastNeighbor[variantJ]; k++) {
					neighbor=neighbors[variantJ][k];
					if (queueElements[neighbor].removed) continue;
					queueElements[neighbor].removed=true;
					for (h=0; h<=lastNeighbor[neighbor]; h++) {
						element=queueElements[neighbors[neighbor][h]];
						if (!queue.remove(element)) continue;
						element.degree--;
						queue.add(element);
					}
				}
			}
			k=-1;
			for (j=0; j<=chromosome2variantsLast; j++) {
				variantJ=chromosome2variants[i][j];
				if (!queueElements[variantJ].removed) chromosome2variants[i][++k]=variantJ;
			}
			chromosome2variants_last[i]=k;
		}
		System.err.println("Fixed "+nChromosomesWithCollision+" chromosomes with collisions ("+((100.0*nChromosomesWithCollision)/N_CHROMOSOMES)+"%)");
		
		// Marking unused variants
		for (i=0; i<=lastVariant; i++) variants[i].isActive=false;
		for (i=0; i<N_CHROMOSOMES; i++) {
			for (j=0; j<=chromosome2variants_last[i]; j++) variants[chromosome2variants[i][j]].isActive=true;
		}
		nVariants_frequent=0; nVariants_rare=0;
		for (i=0; i<=lastVariant; i++) {
			if (!variants[i].isActive) continue;
			if (variants[i].alleleFrequency>=alleleFrequencyThreshold) nVariants_frequent++;
			else nVariants_rare++;
		}
		System.err.println("Total variants after fixing collisions: "+(nVariants_frequent+nVariants_rare)+", frequent="+nVariants_frequent+" rare="+nVariants_rare);
		System.err.println("Variants removed: "+(lastVariant+1-nVariants_frequent-nVariants_rare));
		
		// Printing histogram of nVariants per individual
		histogram = new int[HISTOGRAM_MAX_LENGTH+1];
		Arrays.fill(histogram,0);
		for (i=0; i<N_CHROMOSOMES; i++) {
			j=chromosome2variants_last[i]+1;
			if (j>HISTOGRAM_MAX_LENGTH) j=HISTOGRAM_MAX_LENGTH;
			histogram[j]++;
		}
		bw = new BufferedWriter(new FileWriter(outputDir+"/histogram_nVariantsPerChromosome.txt"));
		for (i=0; i<=HISTOGRAM_MAX_LENGTH; i++) bw.write(i+","+histogram[i]+"\n");
		bw.close();
	}
	
	
	private static final void addEdge(int node1, int node2, int[][] neighbors, int[] lastNeighbor) {
		lastNeighbor[node1]++;
		if (lastNeighbor[node1]==neighbors[node1].length) {
			int[] newArray = new int[neighbors[node1].length<<1];
			System.arraycopy(neighbors[node1],0,newArray,0,neighbors[node1].length);
			neighbors[node1]=newArray;
		}
		neighbors[node1][lastNeighbor[node1]]=node2;
		lastNeighbor[node2]++;
		if (lastNeighbor[node2]==neighbors[node2].length) {
			int[] newArray = new int[neighbors[node2].length<<1];
			System.arraycopy(neighbors[node2],0,newArray,0,neighbors[node2].length);
			neighbors[node2]=newArray;
		}
		neighbors[node2][lastNeighbor[node2]]=node1;
	}
	
	
	private static class QueueElement implements Comparable {
		boolean removed;
		public int variantID, degree;
		
		public QueueElement(int v) { 
			this.variantID=v; 
			this.degree=0;
			this.removed=false;
		}
		
		public boolean equals(Object other) {
			QueueElement otherElement = (QueueElement)other;
			return variantID==otherElement.variantID;
		}
		
		/**
		 * Sorts by increasing degree
		 */
		public int compareTo(Object other) {
			QueueElement otherElement = (QueueElement)other;
			if (degree<otherElement.degree) return -1;
			else if (degree>otherElement.degree) return 1;
			return 0;
		}
	}
	
	
	/**
	 * Assigns an insertion string to every active INS variant. An insertion
	 * string is randomly sampled from $insertions$ according to 
	 * $getInsertionString()$, and it might be reverse-complemented. A variant 
	 * is shortened if its length is greater than any insertion string 
	 * available, otherwise a substring of length exactly equal to the variant
	 * is extracted.
	 *
	 * Remark: the procedure does not ensure that the resulting INS is not 
	 * actually a DUP, i.e. it does not check whether the randomly chosen 
	 * insertion string is a copy of an adjacent substring of the reference.
	 *
	 * Remark: an alternative to sampling from a database of insertions, would
	 * be just copying an adjacent substring, i.e. making every INS a DUP. For
	 * now we assume that SV type in the model is correct and we don't mix DUP
	 * with INS.
	 *
	 * Remark: $insertions$ and $firstInsertionWithLength$ are deallocated after
	 * the procedure completes.
	 */
	private static final void assignInsertionStrings(StringBuilder buffer) {
		int i, j, p;
		int insertionLength, variantLength;
		
		for (i=0; i<=lastVariant; i++) {
			if (i%10000==0) {
				System.err.println("Processed "+i+" variants");
				System.gc();
			}
			if (!variants[i].isActive || variants[i].svType!=2) continue;
			j=getInsertionString(variants[i].length);
			buffer.delete(0,buffer.length());
			insertionLength=insertions[j].length();
			variantLength=variants[i].length;
			if (insertionLength==variantLength) buffer.append(insertions[j]);
			else if (insertionLength>variantLength) {
				p=random.nextInt(insertionLength-variantLength);
				buffer.append(insertions[j],p,p+variantLength);
			}
			else {
				buffer.append(insertions[j]);
				variants[i].length=insertionLength;
			}
			if (random.nextBoolean()) reverseComplement(buffer);
			variants[i].sequence=buffer.toString();
		}
		insertions=null; firstInsertionWithLength=null;
		System.gc();
	}
	
	
	/**
	 * @return if some strings in $insertions$ have length equal to $length$, 
	 * the procedure returns the index of one of them at random; otherwise, if a
	 * string of length $>length$ exists in the same length bin, the procedure
	 * returns the index of one such string at random; otherwise, the procedure 
	 * returns the index of a random string of length $>length$ in the next 
	 * nonempty bin. If everything else fails, the procedure returns the index 
	 * of a random string in the last nonempty bin (so this string is shorter 
	 * than $length$).
	 */
	private static final int getInsertionString(int length) {
		int i, j;
		int bin, firstString, lastString, nBins;
		
		bin=length/BuildModel.SV_LENGTH_QUANTUM;
		nBins=firstInsertionWithLength.length;
		firstString=firstInsertionWithLength[bin];
		if (firstString==-1) {
			while (bin<firstInsertionWithLength.length && firstInsertionWithLength[bin]==-1) bin++;
			if (bin==firstInsertionWithLength.length) {
				bin=length/BuildModel.SV_LENGTH_QUANTUM;
				while (bin>=0 && firstInsertionWithLength[bin]==-1) bin--;
				return firstInsertionWithLength[bin]+random.nextInt(lastInsertion-firstInsertionWithLength[bin]+1);
			}
			else {
				if (bin==firstInsertionWithLength.length-1) lastString=lastInsertion;
				else {
					lastString=-1;
					for (i=bin+1; i<nBins; i++) {
						if (firstInsertionWithLength[i]!=-1) {
							lastString=firstInsertionWithLength[i]-1;
							break;
						}
					}
				}
				return firstInsertionWithLength[bin]+random.nextInt(lastString-firstInsertionWithLength[bin]+1);
			}
		}
		else {
			i=firstString;
			if (bin==firstInsertionWithLength.length-1) lastString=lastInsertion;
			else {
				lastString=-1;
				for (i=bin+1; i<nBins; i++) {
					if (firstInsertionWithLength[i]!=-1) {
						lastString=firstInsertionWithLength[i]-1;
						break;
					}
				}
			}
			while (i<=lastString && insertions[i].length()<length) i++;
			j=lastString;
			while (j>=firstString && insertions[j].length()>length) j--;
			if (i==j) return i;
			else if (i<j) return i+random.nextInt(j-i+1);
			else if (insertions[lastString].length()>=length) return j+1+random.nextInt(lastString-(j+1)+1);
			else return firstString+random.nextInt(lastString-firstString+1);
		}
	}
	
	
	private static final void reverseComplement(StringBuilder buffer) {
		char c, d;
		int i;
		final int length = buffer.length();
		
		buffer.reverse();
		for (i=0; i<length; i++) {
			c=buffer.charAt(i);
			switch (c) {
				case 'a': d='t'; break;
				case 'c': d='g'; break;
				case 'g': d='c'; break;
				case 't': d='a'; break;
				case 'A': d='T'; break;
				case 'C': d='G'; break;
				case 'G': d='C'; break;
				case 'T': d='A'; break;
				default: d='N';
			}
			buffer.setCharAt(i,d);
		}
	}
	
	
	/**
	 * Remark: every variant is applied at exactly its position in the reference
	 * and with exactly its length. In the future we could introduce some 
	 * uncertainty in either or both, although it's not clear what biological
	 * process this would model.
	 *
	 * Remark: the procedure allows an insertion to occur strictly before or 
	 * strictly after the entire reference.
	 *
	 * Remark: the procedure is sequential just for simplicity. In the future, 
	 * we should build a chunk of chromosomes in parallel in a memory buffer,
	 * and flush it periodically to disk.
	 *
	 * @param buffer temporary space.
	 */
	public static final void buildChromosomes(int firstChromosome, int lastChromosome, String outDir, StringBuilder buffer) throws IOException {
		int i, j, k;
		int bin, firstInReference;
		StringBuilder out;
		BufferedWriter bw;
		Variant variant;
		
		out = new StringBuilder(referenceLength);
		for (i=firstChromosome; i<=lastChromosome; i++) {
			out.delete(0,out.length());
			firstInReference=0;
			System.err.println("Building chromosome "+i+" with "+(chromosome2variants_last[i]+1)+" variants...");
			for (j=0; j<=chromosome2variants_last[i]; j++) {
				variant=variants[chromosome2variants[i][j]];
				out.append(reference,firstInReference,variant.position);
				if (variant.svType==0) {  // DEL
					firstInReference=variant.position+variant.length;
				}
				else if (variant.svType==1) {  // DUP
					k=variant.position+variant.length;
					out.append(reference,variant.position,k);
					out.append(reference,variant.position,k);
					firstInReference=k;
				}
				else if (variant.svType==2) {  // INS
					out.append(variant.sequence);
					firstInReference=variant.position;
				}
				else if (variant.svType==3) {  // INV
					buffer.delete(0,buffer.length());
					buffer.append(reference,variant.position,variant.position+variant.length);
					reverseComplement(buffer);
					out.append(buffer);
					firstInReference=variant.position+variant.length;
				}
			}
			if (firstInReference<referenceLength) out.append(reference,firstInReference,referenceLength);
			bw = new BufferedWriter(new FileWriter(outDir+"/chromosome_"+i+".fa"));
			bw.write(">GRCh37_chr1_"+i); bw.newLine();
			bw.write(out.toString()); bw.newLine();
			bw.close();
			System.gc();
		}
	}
	
	
	/**
	 * Remark: the procedure updates variant frequencies with their empirical
	 * values in the simulated chromosomes.
	 */
	private static final void buildVariant2chromosomes() {
		int i, j;
		int variantID;
		
		variant2chromosomes = new int[lastVariant+1][10];  // Arbitrary
		variant2chromosomes_last = new int[lastVariant+1];
		Arrays.fill(variant2chromosomes_last,-1);
		for (i=0; i<N_CHROMOSOMES; i++) {
			for (j=0; j<=chromosome2variants_last[i]; j++) {
				variantID=chromosome2variants[i][j];
				variant2chromosomes_last[variantID]++;
				if (variant2chromosomes_last[variantID]==variant2chromosomes[variantID].length) {
					int[] newArray = new int[variant2chromosomes[variantID].length<<1];
					System.arraycopy(variant2chromosomes[variantID],0,newArray,0,variant2chromosomes[variantID].length);
					variant2chromosomes[variantID]=newArray;
				}
				variant2chromosomes[variantID][variant2chromosomes_last[variantID]]=i;
			}
		}
		for (i=0; i<=lastVariant; i++) variants[i].alleleFrequency=(variant2chromosomes_last[i]+1.0)/N_CHROMOSOMES;
	}
	
	
	private static final String VCF_SEPARATOR = "\t";
	private static final String CHROMOSOME = "1";
	private static final String ID_PREFIX = "simulated_";
	private static final String REF = "N";
	private static final String QUALITY = "999";
	private static final String FILTER = "PASS";
	private static final String INFO_SEPARATOR = ";";
	private static final String GENOTYPE_STR = "GT";
	
	
	/**
	 * Prints a joint VCF that contains the sorted set of distinct variants, and
	 * a VCF for each individual, where an individual consists of two 
	 * consecutive chromosomes in the list of all chromosomes.
	 *
	 * Remark: the procedure does not print the header (could be copied from the
	 * gnomAD-SV file).
	 */
	private static final void buildVCF(String outDir) throws IOException {
		final String SUPPORT_STR = INFO_SEPARATOR+"CHROMOSOMES=";
		int i, j, k;
		int variant, last1, last2;
		BufferedWriter bw;
		NumberFormat formatter = NumberFormat.getInstance();
		formatter.setMinimumFractionDigits(4); formatter.setMaximumFractionDigits(4);
		
		// Joint VCF
		bw = new BufferedWriter(new FileWriter(outDir+"/groundTruth_joint.vcf"));
		printVcfHeader(false,null,bw);
		for (i=0; i<=lastVariant; i++) {
			if (!variants[i].isActive) continue;
			variants[i].writeVCF(i,formatter,bw);
			bw.write(SUPPORT_STR);
			for (j=0; j<variant2chromosomes_last[i]; j++) bw.write(variant2chromosomes[i][j]+",");
			bw.write(variant2chromosomes[i][variant2chromosomes_last[i]]+"");
			bw.newLine();
		}
		bw.close();
		
		// Individuals VCFs
		for (i=0; i<N_CHROMOSOMES; i+=2) {
			bw = new BufferedWriter(new FileWriter(outDir+"/groundTruth_individual_"+i+".vcf"));
			printVcfHeader(true,"individual_"+i,bw);
			last1=chromosome2variants_last[i];
			last2=chromosome2variants_last[i+1];
			j=0; k=0;
			while (j<=last1 && k<=last2) {
				if (chromosome2variants[i][j]<chromosome2variants[i+1][k]) {
					variant=chromosome2variants[i][j];
					variants[variant].writeVCF(variant,formatter,bw);
					bw.write(VCF_SEPARATOR+GENOTYPE_STR+VCF_SEPARATOR+"1|0\n");
					j++;
				}
				else if (chromosome2variants[i+1][k]<chromosome2variants[i][j]) {
					variant=chromosome2variants[i+1][k];
					variants[variant].writeVCF(variant,formatter,bw);
					bw.write(VCF_SEPARATOR+GENOTYPE_STR+VCF_SEPARATOR+"0|1\n");
					k++;
				}
				else {
					variant=chromosome2variants[i][j];
					variants[variant].writeVCF(variant,formatter,bw);
					bw.write(VCF_SEPARATOR+GENOTYPE_STR+VCF_SEPARATOR+"1|1\n");
					j++; k++;
				}
			}
			while (j<=last1) {
				variant=chromosome2variants[i][j];
				variants[variant].writeVCF(variant,formatter,bw);
				bw.write(VCF_SEPARATOR+GENOTYPE_STR+VCF_SEPARATOR+"1|0\n");
				j++;
			}
			while (k<=last2) {
				variant=chromosome2variants[i+1][k];
				variants[variant].writeVCF(variant,formatter,bw);
				bw.write(VCF_SEPARATOR+GENOTYPE_STR+VCF_SEPARATOR+"0|1\n");
				k++;
			}
			bw.close();
		}
	}
	
	
	/**
	 * @param sampleName discarded if $printGenotype=FALSE$.
	 */
	private static final void printVcfHeader(boolean printGenotype, String sampleName, BufferedWriter bw) throws IOException {
		bw.write("##fileformat=VCFv4.2\n");
		bw.write("##contig=<ID=chr1,length="+referenceLength+">\n");
		bw.write("##ALT=<ID=INS,Description=\"Insertion\">\n");
		bw.write("##ALT=<ID=DEL,Description=\"Deletion\">\n");
		bw.write("##ALT=<ID=DUP,Description=\"Duplication\">\n");
		bw.write("##ALT=<ID=INV,Description=\"Inversion\">\n");
		bw.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
		bw.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variation\">\n");
		bw.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
		bw.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">\n");
		bw.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n");
		bw.write("##INFO=<ID=REPEATS_START,Number=1,Type=Integer,Description=\"Repetitive context around the first position\">\n");
		bw.write("##INFO=<ID=REPEATS_END,Number=1,Type=Integer,Description=\"Repetitive context around the last position\">\n");
		bw.write("##INFO=<ID=REPEATS_FRACTION,Number=1,Type=Float,Description=\"Fraction of the SV covered by repeats of any type\">\n");
		bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		bw.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		if (printGenotype) bw.write("\tFORMAT\t"+sampleName);
		bw.newLine();
	}
	
	
	private static final void serialize(String chromosome2variantsFile, String variantsFile) throws IOException {
		int i, j;
		int last;
		BufferedWriter bw;
		
		// $chromosome2variants$
		bw = new BufferedWriter(new FileWriter(chromosome2variantsFile));
		bw.write(N_CHROMOSOMES+"\n");
		for (i=0; i<N_CHROMOSOMES; i++) {
			last=chromosome2variants_last[i];
			bw.write(last+"");
			for (j=0; j<=last; j++) bw.write(","+chromosome2variants[i][j]);
			bw.newLine();
		}
		bw.close();
		
		// $variants$
		bw = new BufferedWriter(new FileWriter(variantsFile));
		bw.write(lastVariant+"\n");
		for (i=0; i<=lastVariant; i++) variants[i].serialize(bw);
		bw.close();
	}
	
	
	public static final void deserialize(String chromosome2variantsFile, String variantsFile) throws IOException {
		int i, j;
		int last;
		String str;
		BufferedReader br;
		String[] tokens;
		
		// $chromosome2variants$
		br = new BufferedReader(new FileReader(chromosome2variantsFile));
		N_CHROMOSOMES=Integer.parseInt(br.readLine());
		chromosome2variants = new int[N_CHROMOSOMES][0];
		chromosome2variants_last = new int[N_CHROMOSOMES];
		for (i=0; i<N_CHROMOSOMES; i++) {
			str=br.readLine();
			tokens=str.split(",");
			last=Integer.parseInt(tokens[0]);
			chromosome2variants_last[i]=last;
			chromosome2variants[i] = new int[last+1];
			for (j=0; j<=last; j++) chromosome2variants[i][j]=Integer.parseInt(tokens[1+j]);
		}
		br.close();
		
		// $variants$
		br = new BufferedReader(new FileReader(variantsFile));
		lastVariant=Integer.parseInt(br.readLine());
		variants = new Variant[lastVariant+1];
		for (i=0; i<=lastVariant; i++) {
			variants[i] = new Variant();
			variants[i].deserialize(br);
		}
		br.close();
	}
	
	
	private static class Variant implements Comparable {
		public int svType, length;
		public double alleleFrequency;
		public int position;  // First position in the reference
		public int intervalType;  // Type of repeat to which $position$ belongs
		public String sequence;  // Allocated only for INS
		public boolean isActive;
		
		/**
		 * Type of repeat intervals near the first/last position of the variant
		 */
		public int intervalType_start, intervalType_end;
		
		/**
		 * Number of bps of the variant that are covered by repeat intervals
		 */
		public int intervalSurface;
		
		
		/**
		 * Sorts by increasing position
		 */
		public int compareTo(Object other) {
			Variant otherVariant = (Variant)other;
			if (position<otherVariant.position) return -1;
			else if (position>otherVariant.position) return 1;
			return 0;
		}
		
		/**
		 * @param maxDistance two variants are assumed to collide iff they are 
		 * up to this far apart from each other.
		 */
		public final boolean collidesWith(Variant other, int maxDistance) {
			int first1, first2, last1, last2;
			
			if (svType==2) {  // INS
				if (other.svType==2) return Math.abs(position-other.position)<=maxDistance;
				else return (position>=other.position-maxDistance && position<=other.position+other.length-1+maxDistance);
			}
			else {
				if (other.svType==2) {  // INS
					return other.position>=position-maxDistance && other.position<=position+length-1+maxDistance;
				}
				else {
					first1=Math.min(position,other.position);
					last1=Math.max(position+length-1,other.position+other.length-1);
					if (first1<=last1) return true;
					first1=position; last1=position+length-1;
					first2=other.position; last2=other.position+other.length-1;
					return (first1>last2 && first1<=last2+maxDistance) || (last1<first2 && last1>=first2-maxDistance);
				}
			}
		}
		
		public String toString() {
			return "svType="+svType+" position="+position+" length="+length+" isActive="+isActive+" intervalType="+intervalType+" alleleFrequency="+alleleFrequency;
		}
		
		public void writeVCF(int id, NumberFormat formatter, BufferedWriter bw) throws IOException {
			bw.write(CHROMOSOME); bw.write(VCF_SEPARATOR);
			bw.write((position+1)+VCF_SEPARATOR);
			bw.write(ID_PREFIX); bw.write(id+VCF_SEPARATOR);
			bw.write(REF); bw.write(VCF_SEPARATOR);
			if (svType==2) bw.write(sequence);
			else bw.write("<"+BuildModel.SV_TYPES[svType]+">");
			bw.write(VCF_SEPARATOR);
			bw.write(QUALITY); bw.write(VCF_SEPARATOR);
			bw.write(FILTER); bw.write(VCF_SEPARATOR);
			bw.write("END="+(svType==2?position:(position+length))); bw.write(INFO_SEPARATOR);
			bw.write("SVTYPE="); bw.write(BuildModel.SV_TYPES[svType]); bw.write(INFO_SEPARATOR);
			bw.write("SVLEN="+length); bw.write(INFO_SEPARATOR);
			bw.write("AF="+alleleFrequency); bw.write(INFO_SEPARATOR);
			bw.write("REPEATS_START="+intervalType_start); bw.write(INFO_SEPARATOR);
			bw.write("REPEATS_END="+intervalType_end); bw.write(INFO_SEPARATOR);
			bw.write("REPEATS_FRACTION="+formatter.format(((double)intervalSurface)/length));
		}
		
		public final void serialize(BufferedWriter bw) throws IOException {
			bw.write(svType+"\n");
			bw.write(intervalType+"\n");
			bw.write(length+"\n");
			bw.write(alleleFrequency+"\n");
			bw.write(position+"\n");
			bw.write(isActive+"\n");
			bw.write(sequence==null?"null":sequence); bw.newLine();
			bw.write(intervalType_start+"\n");
			bw.write(intervalType_end+"\n");
			bw.write(intervalSurface+"\n");
		}
		
		public final void deserialize(BufferedReader br) throws IOException {
			svType=Integer.parseInt(br.readLine());
			intervalType=Integer.parseInt(br.readLine());
			length=Integer.parseInt(br.readLine());
			alleleFrequency=Double.parseDouble(br.readLine());
			position=Integer.parseInt(br.readLine());
			isActive=Boolean.parseBoolean(br.readLine());
			sequence=br.readLine();
			if (sequence.equalsIgnoreCase("null")) sequence=null;
			intervalType_start=Integer.parseInt(br.readLine());
			intervalType_end=Integer.parseInt(br.readLine());
			intervalSurface=Integer.parseInt(br.readLine());
		}
	}

}