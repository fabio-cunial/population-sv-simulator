import java.util.Arrays;
import java.io.*;


/**
 * Derives from gnomAD-SV a model for simulating structural variants on chr1.
 */
public class BuildModel {
	/**
	 * Distance thresholds
	 */
	public static final int MAX_DISTANCE_REPEATMASKER = 25;
	public static final int MIN_OVERLAP_TRF = 10;
	public static final int MIN_OVERLAP_INTERVALS = 25;
	public static final int MAX_DISTANCE_SR = 25;  // Suggested by Ryan Collins
	private static final int MAX_DISTANCE_PE = 250; // Suggested by Ryan Collins
	private static final int MAX_DISTANCE_OTHER = 1000;  // Suggested by Collins
	public static final int MAX_DISTANCE_SEGDUPS = 25;


	public static void main(String[] args) throws IOException {
		final String REPEATMASKER_FILE = args[0];
		final int REPEATMASKER_FILE_NROWS = Integer.parseInt(args[1]);
		final String TRF_FILE = args[2];
		final int TRF_FILE_NROWS = Integer.parseInt(args[3]);
		final String SEGDUPS_FILE = args[4];
		final int SEGDUPS_FILE_NROWS = Integer.parseInt(args[5]);
        final String REFERENCE_FAI = args[6];
        final String REFERENCE_CYTOBANDS = args[7];
		final String GNOMADSV_FILE = args[8];
		final boolean ONLY_PASS = Integer.parseInt(args[9])==1;
		final String OUTPUT_DIR = args[10];
		final boolean VERBOSE = Integer.parseInt(args[11])==1;
		
		loadRepeatMaskerAnnotations(REPEATMASKER_FILE,REPEATMASKER_FILE_NROWS,MAX_DISTANCE_REPEATMASKER,false);
		loadTrfAnnotations(TRF_FILE,TRF_FILE_NROWS,MIN_OVERLAP_TRF,false);
		buildRepeatIntervals(MIN_OVERLAP_INTERVALS);
		loadSegdups(SEGDUPS_FILE,SEGDUPS_FILE_NROWS,MAX_DISTANCE_SEGDUPS);
		//alleleFrequencyVsRepeatContext(GNOMADSV_FILE,MAX_DISTANCE_SR,MAX_DISTANCE_PE,MAX_DISTANCE_OTHER,MAX_DISTANCE_SEGDUPS,ONLY_PASS,"./");
        loadChromosomeLengths(REFERENCE_FAI,REFERENCE_CYTOBANDS);
		collectStatistics(GNOMADSV_FILE,MAX_DISTANCE_SR,MAX_DISTANCE_PE,MAX_DISTANCE_OTHER,MAX_DISTANCE_SEGDUPS,ONLY_PASS,VERBOSE);
		makeStatisticsCumulative();
		serializeStatistics(OUTPUT_DIR);
		serializeRepeatIntervals(OUTPUT_DIR+"/repeatIntervals.txt");
		serializeSegdups(OUTPUT_DIR+"/segdups.txt");
	}
	
	
	
	
	// ---------------------- REPEAT LOADING PROCEDURES ------------------------
	
	/**
	 * Data structures encoding all repeat annotations over all chromosomes. 
	 * Columns are annotations. Rows:
	 *
	 * 0: start in chromosome;
	 * 1: end in chromosome;
	 * 2: start in FWD repeat consensus (if it is known);
	 * 3: end in FWD repeat consensus (if it is known);
	 * 4: orientation (1=FWD, 0=RC) (if it is known).
	 */
	public static int[][] repeatTrack;
	public static int lastRepeatTrack = -1;
	
	/**
	 * For every record of $repeatTrack$.
	 */
	public static String[] chromosomeIDs;
	public static boolean[] isSatellite;
	public static String[][] repeatTrackIDs;  // 0=repeat; 1=family.
	
	/**
	 * Maximal ranges of overlapping repeat annotations. Rows are intervals with
	 * small overlap. Columns: 1=start in chromosome; 2=end in chromosome;
	 * 0=type of interval (0=non-satellite; 1=satellite; 2=both).
	 */
	public static int[][] repeatIntervals;
	public static int lastRepeatInterval;
	private static String[] repeatIntervals_chromosomes;
	private static int[] chromosome2firstInterval;
	
	public static final int[][] ADD_TYPE = new int[][] {{0,2,2},{2,1,2},{2,2,2}};
	
	
	/**
	 * Remark: the procedure loads all repeats from all chromosomes.
	 *
	 * @param path assumed to be a cleaned version of the original RepeatMasker 
	 * file from the UCSC Genome Browser, without header and with runs of spaces
	 * replaced by a single comma (the rows of this file are already grouped by 
	 * contig and sorted by starting position); contigs are assumed to be sorted
	 * as follows:
	 *
	 * chr1,chr10,...,chr19,  chr2,...,chr22,  chr3,...,chr9,  chrX,chrY, other;
	 *
	 * contigs in "other" are not loaded;
	 * @param maxDistance for merging nearby partial matches of the same
	 * repeat.
	 */
	public static final void loadRepeatMaskerAnnotations(String path, int nRecords, int maxDistance, boolean verbose) throws IOException {
		boolean found, orientation;
		int i, j, k;
		int referenceStart, referenceEnd, repeatStart, repeatEnd, length;
		int nSatellites;
		long repeatSurface, satelliteSurface;
		String str, chr;
		BufferedReader br;
		String[] tokens;
		
		chromosomeIDs = new String[nRecords];
		repeatTrack = new int[6][nRecords];
		isSatellite = new boolean[nRecords];
		repeatTrackIDs = new String[nRecords][2];
		System.err.print("Loading "+nRecords+" RepeatMasker annotations... ");
		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); found=false; lastRepeatTrack=-1; 
		nSatellites=0; repeatSurface=0; satelliteSurface=0;
		while (str!=null) {
			tokens=str.split(",");
			if (tokens[4].length()>5 || tokens[4].equalsIgnoreCase("chrM")) break;
			lastRepeatTrack++;
			chromosomeIDs[lastRepeatTrack]=tokens[4];
			referenceStart=Integer.parseInt(tokens[5])-1;
			referenceEnd=Integer.parseInt(tokens[6])-1;
			length=referenceEnd-referenceStart+1;
			repeatSurface+=length;
			repeatTrack[0][lastRepeatTrack]=referenceStart;
			repeatTrack[1][lastRepeatTrack]=referenceEnd;
			orientation=tokens[8].charAt(0)=='+';
			if (orientation) {
				repeatStart=Integer.parseInt(tokens[11])-1;
				repeatEnd=Integer.parseInt(tokens[12])-1;
			}
			else {
				repeatEnd=Integer.parseInt(tokens[12])-1;
				repeatStart=Integer.parseInt(tokens[13])-1;
			}
			repeatTrack[2][lastRepeatTrack]=repeatStart;
			repeatTrack[3][lastRepeatTrack]=repeatEnd;
			repeatTrack[4][lastRepeatTrack]=orientation?1:0;
			repeatTrackIDs[lastRepeatTrack][0]=tokens[9].toLowerCase();
			repeatTrackIDs[lastRepeatTrack][1]=tokens[10].toLowerCase();
			isSatellite[lastRepeatTrack] = repeatTrackIDs[lastRepeatTrack][1].indexOf("satellite")>=0 ||
				                           repeatTrackIDs[lastRepeatTrack][1].indexOf("sat")>=0 ||
										   repeatTrackIDs[lastRepeatTrack][1].indexOf("simple")>=0 ||
										   repeatTrackIDs[lastRepeatTrack][1].indexOf("low_complexity")>=0;
			if (isSatellite[lastRepeatTrack]) {
				nSatellites++;
				satelliteSurface+=length;
			}
			str=br.readLine();
		}
		br.close();
		System.err.println("DONE. "+repeatSurface+" bps annotated as repeats over the entire genome.");
		System.err.println(nSatellites+" satellites ("+((100.0*nSatellites)/(lastRepeatTrack+1))+"% of all records, "+((100.0*satelliteSurface)/repeatSurface)+"% of total repeat surface)");
		
		// Merging adjacent annotations if they are a few bps apart in both the
		// reference and the repeat.
		System.err.print("Compacting "+(lastRepeatTrack+1)+" RepeatMasker annotations... ");
		j=0;
		for (i=1; i<=lastRepeatTrack; i++) {
			if ( repeatTrackIDs[i][0].equals(repeatTrackIDs[j][0]) && 
				 repeatTrackIDs[i][1].equals(repeatTrackIDs[j][1]) && 
			     repeatTrack[4][i]==repeatTrack[4][j] &&
				 chromosomeIDs[i].equalsIgnoreCase(chromosomeIDs[j]) &&
				 Math.abs(repeatTrack[0][i]-repeatTrack[1][j])<=maxDistance &&
				 ( isSatellite[i] ||
				   (repeatTrack[4][i]==1 && Math.abs(repeatTrack[2][i]-repeatTrack[3][j])<=maxDistance) ||
				   (repeatTrack[4][i]==0 && Math.abs(repeatTrack[3][i]-repeatTrack[2][j])<=maxDistance)
				 )
		       ) {
				if (verbose) {
					System.err.println("Merging:");				   
					System.err.println(printRepeatMaskerAnnotation(j));
					System.err.println(printRepeatMaskerAnnotation(i));
				}
				repeatTrack[1][j]=repeatTrack[1][i];
				repeatTrack[2][j]=Math.min(repeatTrack[2][j],repeatTrack[2][i]);
				repeatTrack[3][j]=Math.max(repeatTrack[3][j],repeatTrack[3][i]);	
			}
			else {
				j++;
				for (k=0; k<repeatTrack.length; k++) repeatTrack[k][j]=repeatTrack[k][i];
				isSatellite[j]=isSatellite[i];
				System.arraycopy(repeatTrackIDs[i],0,repeatTrackIDs[j],0,repeatTrackIDs[i].length);
				chromosomeIDs[j]=chromosomeIDs[i];
			}
		}
		lastRepeatTrack=j;
		System.err.println(" DONE. "+(lastRepeatTrack+1)+" annotations after compaction.");		
	}
	
	
	public static final String printRepeatMaskerAnnotation(int i) {		
		return chromosomeIDs[i]+"["+repeatTrack[0][i]+".."+repeatTrack[1][i]+"] =="+(repeatTrack[4][i]==1?"FWD":"REV")+"== ["+repeatTrack[2][i]+".."+repeatTrack[3][i]+"]"+repeatTrackIDs[i][0]+"<-"+repeatTrackIDs[i][1]+" "+(isSatellite[i]?"SATELLITE":"TRANSPOSON");
	}
	
	
	/**
	 * Loads all Tandem Repeats Finder annotations (some of which might be very
	 * short), from all chromosomes, and merges them with existing satellite 
	 * annotations.
	 *
	 * Remark: the "simple repeats" track in [1], taken from the UCSC Genome 
	 * Browser, comes from TRF. The microsatellite track in the UCSC Genome 
	 * Browser is a subset of the simple repeats track.
	 *
	 * [1] Zhao et al. "Expectations and blind spots for structural variation 
	 * detection from long-read assemblies and short-read genome sequencing 
	 * technologies." The American Journal of Human Genetics 2021.
	 * 
	 * @param path a cleaned version of the TRF file from UCSC Genome Browser, 
	 * with runs of tabs replaced by commas (the rows of this file are already
	 * grouped by contig and sorted by starting position); contigs are assumed
	 * to be sorted as follows (contigs in "other" are not loaded):
	 *
	 * chr1,chr10,...,chr19,  chr2,...,chr22,  chr3,...,chr9,  chrX,chrY, other.
	 */
	public static final void loadTrfAnnotations(String path, int nRecords, int minOverlap, boolean verbose) throws IOException {
		int i, j;
		int referenceStart, referenceEnd, length, last;
		long repeatSurface;
		String str, chr;
		BufferedReader br;
		Tuple tmpTuple;
		int[][] trfTrack;
		int lastTrfTrack;
		String[] tokens, chrIDs;
		Tuple[] tuples;

		// Loading annotations
		System.err.print("Loading "+nRecords+" TRF annotations... ");
		trfTrack = new int[nRecords][2]; chrIDs = new String[nRecords];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine(); lastTrfTrack=-1;
		repeatSurface=0;
		while (str!=null) {
			tokens=str.split(",");
			lastTrfTrack++;
			chrIDs[lastTrfTrack]=tokens[0];
			referenceStart=Integer.parseInt(tokens[1])-1;
			referenceEnd=Integer.parseInt(tokens[2])-1;
			repeatSurface+=referenceEnd-referenceStart+1;
			trfTrack[lastTrfTrack][0]=referenceStart;
			trfTrack[lastTrfTrack][1]=referenceEnd;
			str=br.readLine();
		}
		br.close();
		System.err.println("DONE. "+repeatSurface+" bps annotated as tandem over the entire genome.");

		// Merging consecutive annotations if they overlap in the reference,
	 	// regardless of their period. The result of a merge might be short.
		System.err.print("Compacting "+(lastTrfTrack+1)+" TRF annotations... ");
		j=0;
		for (i=1; i<=lastTrfTrack; i++) {
			if (chrIDs[i].equalsIgnoreCase(chrIDs[j]) && trfTrack[i][0]<=trfTrack[j][1]-minOverlap) {
				if (verbose) {
					System.err.println("Merging:");
					System.err.println(printTrfAnnotation(j,trfTrack,chrIDs));
					System.err.println(printTrfAnnotation(i,trfTrack,chrIDs));
				}
				if (trfTrack[i][1]>trfTrack[j][1]) trfTrack[j][1]=trfTrack[i][1];
			}
			else {
				j++;
				System.arraycopy(trfTrack[i],0,trfTrack[j],0,trfTrack[i].length);
				chrIDs[j]=chrIDs[i];
			}
		}
		lastTrfTrack=j;
		System.err.println(" DONE. "+(lastTrfTrack+1)+" TRF annotations after compaction.");

		// Stores in the global data structures, the union of existing and TRF
		// intervals. New intervals might be built by merging existing and TRF
		// intervals that overlap.
		// Remark: $trfTrack$ is null after this section completes.
		System.err.println("Merging "+(lastRepeatTrack+lastTrfTrack+2)+" total repeat annotations...");
		tuples = new Tuple[lastRepeatTrack+lastTrfTrack+2];
		for (i=0; i<=lastRepeatTrack; i++) tuples[i] = new Tuple(chromosomeIDs[i],repeatTrack[0][i],repeatTrack[1][i],repeatTrack[2][i],repeatTrack[3][i],repeatTrack[4][i],isSatellite[i],repeatTrackIDs[i][0],repeatTrackIDs[i][1]);
		last=i-1;
		for (i=0; i<=lastTrfTrack; i++) tuples[++last] = new Tuple(chrIDs[i],trfTrack[i][0],trfTrack[i][1],-1,-1,1,true,Tuple.ARTIFICIAL_SHORT_PERIOD,"");
		trfTrack=null;
		Arrays.sort(tuples,0,last+1);
		j=0;
		for (i=1; i<=last; i++) {
			if ( (tuples[i].repeat.equals(Tuple.ARTIFICIAL_SHORT_PERIOD) || tuples[j].repeat.equals(Tuple.ARTIFICIAL_SHORT_PERIOD)) &&
				 tuples[i].chromosome.equals(tuples[j].chromosome) &&
	             tuples[i].referenceStart<=tuples[j].referenceEnd-minOverlap
			   ) tuples[j].merge(tuples[i]);
			else {
				j++;
				tmpTuple=tuples[j];
				tuples[j]=tuples[i];
				tuples[i]=tmpTuple;
			}
		}
		last=j;
		System.err.println("DONE. "+(last+1)+" annotations after the merge.");	

		// Updating data structures
		if (repeatTrack[0].length<last+1) repeatTrack = new int[5][last+1];
		if (isSatellite.length<last+1) isSatellite = new boolean[last+1];
		if (repeatTrackIDs.length<last+1) repeatTrackIDs = new String[last+1][2];
		if (chromosomeIDs.length<last+1) chromosomeIDs = new String[last+1];
		for (i=0; i<=last; i++) {
			repeatTrack[0][i]=tuples[i].referenceStart;
			repeatTrack[1][i]=tuples[i].referenceEnd;
			repeatTrack[2][i]=tuples[i].repeatStart;
			repeatTrack[3][i]=tuples[i].repeatEnd;
			repeatTrack[4][i]=tuples[i].orientation;
			isSatellite[i]=tuples[i].isSatellite;
			repeatTrackIDs[i][0]=tuples[i].repeat;
			repeatTrackIDs[i][1]=tuples[i].repeatFamily;
			chromosomeIDs[i]=tuples[i].chromosome;
		}
		lastRepeatTrack=last;
		System.err.println((lastRepeatTrack+1)+" merged annotations");
		if (verbose) {
			System.err.println("Merged annotations:");
			for (i=0; i<=lastRepeatTrack; i++) System.err.println(printRepeatMaskerAnnotation(i));
		}
	}
	
	
	public static final String printTrfAnnotation(int i, int[][] trfTrack, String[] chrIDs) {
		return chrIDs[i]+"["+trfTrack[i][0]+".."+trfTrack[i][1]+"]";
	}


	private static class Tuple implements Comparable {
		public static final String ARTIFICIAL_SHORT_PERIOD = "ARTIFICIAL_SHORT_PERIOD";

		public int referenceStart, referenceEnd, repeatStart, repeatEnd, orientation;
		public boolean isSatellite;
		public String chromosome, repeat, repeatFamily;

		public Tuple(String chr, int gs, int ge, int rs, int re, int o, boolean s, String rp, String rf) {
			set(chr,gs,ge,rs,re,o,s,rp,rf);
		}

		public final void set(String chr, int gs, int ge, int rs, int re, int o, boolean s, String rp, String rf) {
			chromosome=chr; referenceStart=gs; referenceEnd=ge;
			repeatStart=rs; repeatEnd=re; orientation=o;
			isSatellite=s;
			repeat=rp; repeatFamily=rf;
		}

		/**
		 * Adds $otherTuple$ to this tuple, assuming that at least one of the
		 * tuples is an artificial short-period interval.
		 */
		public final void merge(Tuple otherTuple) {
			boolean changed = false;

			if (otherTuple.referenceStart<referenceStart) {
				referenceStart=otherTuple.referenceStart;
				changed=true;
			}
			if (otherTuple.referenceEnd>referenceEnd) {
				referenceEnd=otherTuple.referenceEnd;
				changed=true;
			}
			isSatellite|=otherTuple.isSatellite;
			if (changed) {
				repeatStart=-1; repeatEnd=-1;
				if (!repeat.equals(ARTIFICIAL_SHORT_PERIOD) && otherTuple.repeat.equals(ARTIFICIAL_SHORT_PERIOD)) {
					repeat=ARTIFICIAL_SHORT_PERIOD;
					repeatFamily="";
					orientation=1;
				}
				else if (repeat.equals(ARTIFICIAL_SHORT_PERIOD) && !otherTuple.repeat.equals(ARTIFICIAL_SHORT_PERIOD)) { /* NOP */ }
			}
		}

		/**
		 * Order: $chromosome,referenceStart,referenceEnd$, where $chromosome$
		 * is in the same order as $VCF2gfa.string2contig()$.
		 */
		public int compareTo(Object other) {
			Tuple otherTuple = (Tuple)other;
			final int i = chromosome.compareTo(otherTuple.chromosome);
			if (i!=0) return i;
			if (referenceStart<otherTuple.referenceStart) return -1;
			else if (referenceStart>otherTuple.referenceStart) return 1;
			if (referenceEnd<otherTuple.referenceEnd) return -1;
			else if (referenceEnd>otherTuple.referenceEnd) return 1;
			return 0;
		}

		public String toString() {
			return chromosome+"["+referenceStart+".."+referenceEnd+"] =="+(orientation==1?"FWD":"REV")+"== ["+repeatStart+".."+repeatEnd+"]"+repeat+"<-"+repeatFamily;
		}
	}
	
	
	/**
	 * Builds $repeatIntervals$ by merging all maximal ranges of overlapping
	 * elements of $repeatTrack$, regardless of their repeat type.
	 */
	public static final void buildRepeatIntervals(int minOverlap) {
		int i, j;
		
		System.err.println("Building repeat intervals from "+(lastRepeatTrack+1)+" total repeat annotations...");
		repeatIntervals = new int[lastRepeatTrack][3];
		repeatIntervals_chromosomes = new String[lastRepeatTrack];
		repeatIntervals_chromosomes[0]=chromosomeIDs[0];
		repeatIntervals[0][0]=isSatellite[0]?1:0;
		repeatIntervals[0][1]=repeatTrack[0][0];
		repeatIntervals[0][2]=repeatTrack[1][0];
		lastRepeatInterval=0;
		for (i=1; i<=lastRepeatTrack; i++) {
			if ( chromosomeIDs[i].equals(chromosomeIDs[lastRepeatInterval]) &&
	             repeatTrack[0][i]<=repeatIntervals[lastRepeatInterval][2]-minOverlap
			   ) {
				if (repeatTrack[1][i]>repeatIntervals[lastRepeatInterval][2]) repeatIntervals[lastRepeatInterval][2]=repeatTrack[1][i];
				repeatIntervals[lastRepeatInterval][0]=ADD_TYPE[repeatIntervals[lastRepeatInterval][0]][isSatellite[i]?1:0];
			}
			else {
				lastRepeatInterval++;
				repeatIntervals[lastRepeatInterval][0]=isSatellite[i]?1:0;
				repeatIntervals[lastRepeatInterval][1]=repeatTrack[0][i];
				repeatIntervals[lastRepeatInterval][2]=repeatTrack[1][i];
				repeatIntervals_chromosomes[lastRepeatInterval]=chromosomeIDs[i];
			}
		}
		chromosome2firstInterval = new int[24];
		Arrays.fill(chromosome2firstInterval,-1);
		chromosome2firstInterval[Integer.parseInt(repeatIntervals_chromosomes[0].substring(3))-1]=0;
		for (i=1; i<=lastRepeatInterval; i++) {
			if (repeatIntervals_chromosomes[i].length()>5 || !repeatIntervals_chromosomes[i].substring(0,3).equalsIgnoreCase("chr")) break;
			if (!repeatIntervals_chromosomes[i].equalsIgnoreCase(repeatIntervals_chromosomes[i-1])) chromosome2firstInterval[Integer.parseInt(repeatIntervals_chromosomes[i].substring(3))-1]=i;
		}
		j=lastRepeatInterval;
		for (i=chromosome2firstInterval.length-1; i>=0; i--) {
			if (chromosome2firstInterval[i]==-1) chromosome2firstInterval[i]=j;
			else j=chromosome2firstInterval[i];
		}
		System.err.println("DONE. "+(lastRepeatInterval+1)+" repeat intervals built over the entire genome.");
	}
	
	
	/**
	 * Stores only chr1 intervals.
	 */
	public static final void serializeRepeatIntervals(String outFile) throws IOException {
		int i;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outFile));
		for (i=0; i<=lastRepeatInterval; i++) {
			if (!repeatIntervals_chromosomes[i].equalsIgnoreCase("chr1")) break;
		}
		bw.write(i+"\n");
		for (i=0; i<=lastRepeatInterval; i++) {
			if (!repeatIntervals_chromosomes[i].equalsIgnoreCase("chr1")) break;
			bw.write(repeatIntervals[i][0]+","+repeatIntervals[i][1]+","+repeatIntervals[i][2]);
			bw.newLine();
		}
		bw.close();
	}
	
	
	public static final void deserializeRepeatIntervals(String inFile) throws IOException {
		int i, j;
		String str;
		BufferedReader br;
		String[] tokens;
		
		br = new BufferedReader(new FileReader(inFile));
		str=br.readLine();
		lastRepeatInterval=Integer.parseInt(str)-1;
		repeatIntervals = new int[lastRepeatInterval+1][3];
		for (i=0; i<=lastRepeatInterval; i++) {
			str=br.readLine();
			tokens=str.split(",");
			for (j=0; j<3; j++) repeatIntervals[i][j]=Integer.parseInt(tokens[j]);
		}
		br.close();
	}

	
	
	
	// ------------------------ SEGDUPS PROCEDURES -----------------------------
	
	/**
	 * Rows: maximal overlapping or adjacent segdup intervals. Columns: 0=start 
	 * in chromosome; 1=end in chromosome.
	 */
	public static int[][] segdups;
	public static int lastSegdup;
	public static String[] segdups_chromosomeIDs;
	private static int[] chromosome2firstSegdup;
	
	
	/**
	 * Remark: the procedure loads all segdup intervals, regardless of their
	 * edit distance from their mate.
	 *
	 * @param path file downloaded from the UCSC Table Browser (this is not 
	 * sorted); chromosome IDs have the form "chrZ";
	 * @param maxDistance adjacent segdups at this distance or smaller are 
	 * merged. Overlapping segupds are always merged.
	 */
	public static final void loadSegdups(String path, int nSegdups, int maxDistance) throws IOException {
		int i, j;
		String str;
		BufferedReader br;
		String[] tokens;
		Segdup[] array;
		
		// Loading file
		array = new Segdup[nSegdups];
		br = new BufferedReader(new FileReader(path));
		str=br.readLine();  // Skipping header
		str=br.readLine(); i=-1;
		while (str!=null) {
			tokens=str.split("\t");
			array[++i] = new Segdup(tokens[1],Integer.parseInt(tokens[2])-1,Integer.parseInt(tokens[3])-2);
			str=br.readLine();
		}
		br.close();
		
		// Computing maximal overlapping intervals
		Arrays.sort(array);
		j=0; lastSegdup=array.length-1;
		for (i=1; i<=lastSegdup; i++) {
			if (array[i].first>array[j].last+maxDistance) array[++j]=array[i];
			else if (array[i].last>array[j].last) array[j].last=array[i].last;
		}
		lastSegdup=j;
		
		// Saving
		segdups = new int[lastSegdup+1][2];
		segdups_chromosomeIDs = new String[lastSegdup+1];
		for (i=0; i<=lastSegdup; i++) {
			segdups[i][0]=array[i].first;
			segdups[i][1]=array[i].last;
			segdups_chromosomeIDs[i]=array[i].chromosome;
		}
		
		// Building $chromosome2firstSegdup$.
		chromosome2firstSegdup = new int[24];
		Arrays.fill(chromosome2firstSegdup,-1);
		chromosome2firstSegdup[Integer.parseInt(segdups_chromosomeIDs[0].substring(3))-1]=0;
		for (i=1; i<=lastSegdup; i++) {
			if (segdups_chromosomeIDs[i].length()>5 || !segdups_chromosomeIDs[i].substring(0,3).equalsIgnoreCase("chr")) break;
			if (!segdups_chromosomeIDs[i].equalsIgnoreCase(segdups_chromosomeIDs[i-1])) chromosome2firstSegdup[Integer.parseInt(segdups_chromosomeIDs[i].substring(3))-1]=i;
		}
		j=lastSegdup;
		for (i=chromosome2firstSegdup.length-1; i>=0; i--) {
			if (chromosome2firstSegdup[i]==-1) chromosome2firstSegdup[i]=j;
			else j=chromosome2firstSegdup[i];
		}
	}
	
	
	/**
	 * Stores every segdup interval in chr1, regardless of whether its mate is
	 * also in chr1 or not.
	 */
	public static final void serializeSegdups(String outFile) throws IOException {
		int i;
		BufferedWriter bw;
		
		bw = new BufferedWriter(new FileWriter(outFile));
		for (i=0; i<=lastSegdup; i++) {
			if (!segdups_chromosomeIDs[i].equalsIgnoreCase("chr1")) break;
		}
		bw.write(i+"\n");
		for (i=0; i<=lastSegdup; i++) {
			if (!segdups_chromosomeIDs[i].equalsIgnoreCase("chr1")) break;
			bw.write(segdups[i][0]+","+segdups[i][1]);
			bw.newLine();
		}
		bw.close();
	}
	
	
	public static final void deserializeSegdups(String inFile) throws IOException {
		int i, j;
		String str;
		BufferedReader br;
		String[] tokens;
		
		br = new BufferedReader(new FileReader(inFile));
		str=br.readLine();
		lastSegdup=Integer.parseInt(str)-1;
		segdups = new int[lastSegdup+1][2];
		for (i=0; i<=lastSegdup; i++) {
			str=br.readLine();
			tokens=str.split(",");
			segdups[i][0]=Integer.parseInt(tokens[0]);
			segdups[i][1]=Integer.parseInt(tokens[1]);
		}
		br.close();
	}
	
	
	private static class Segdup implements Comparable {
		public String chromosome;
		public int first, last;
		
		public Segdup(String c, int f, int l) {
			this.chromosome=c;
			this.first=f; this.last=l;
		}
		
		/**
		 * Sorts by $chromosome$ (lex order), $first,last$.
		 */
		public int compareTo(Object other) {
			Segdup otherSegdup = (Segdup)other;
			int out = chromosome.compareTo(otherSegdup.chromosome);
			if (out!=0) return out;
			if (first<otherSegdup.first) return -1;
			else if (first>otherSegdup.first) return 1;
			if (last<otherSegdup.last) return -1;
			else if (last>otherSegdup.last) return 1;
			return 0;
		}
	}	
	
    
    
    
    // ----------------------- SV POSITION PROCEDURES --------------------------
    
    private static final int N_CHROMOSOMES = 24;
    private static final int N_AUTOSOMES = 22;
    
    /**
     * Parameters from the gnomAD-SV paper.
     */
    private static final int SV_POSITION_QUANTUM = 100000;
	public static final int SMOOTHING_RADIUS_POSITIONS = 5;  // In bins
	public static final int SMOOTHING_WINDOW_POSITIONS = 1+((SMOOTHING_RADIUS_POSITIONS)<<1);
    public static final double TELOMERE_ARM_FRACTION = 0.05;
    public static final double CENTROMERE_ARM_FRACTION = 0.05;
    
    /**
     * Length of each chromosome and of its P-arm
     */
    private static int[] chromosomeLengths, chromosomePLengths;
    
    /**
     * (SV type, Chromosome, Bin) -> count(SVs that start inside the bin)
     */
    private static double[][][] svs_per_position;
    
    /**
     * (SV type, Bin) -> Prob(SV starts inside the bin)
     */
    private static double[][] statistics_position;
    
    
    /**
     * Loads $chromosomeLengths,chromosomePLengths$.
     *
     * Remark: the procedure initializes also $svs_per_position$.
     *
     * @param cytobandFile the cytoband file distributed with the reference. We 
     * assume that every interval is zero-based and of the form $[a..b)$.
     */
    private static final void loadChromosomeLengths(String faiFile, String cytobandFile) throws IOException {
        char c;
        int i, j;
        int length;
        String str;
        BufferedReader br;
        String[] tokens;
        
        chromosomeLengths = new int[N_CHROMOSOMES];
        br = new BufferedReader(new FileReader(faiFile));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            i=str2chromosomeID(tokens[0]);
            if (i==-1) {
                str=br.readLine();
                continue;
            }
            chromosomeLengths[i]=Integer.parseInt(tokens[1]);
            str=br.readLine();
        }
        br.close();
        svs_per_position = new double[SV_TYPES.length][N_CHROMOSOMES][0];
        for (i=0; i<SV_TYPES.length; i++) {
            for (j=0; j<N_CHROMOSOMES; j++) svs_per_position[i][j] = new double[(chromosomeLengths[j]+SV_POSITION_QUANTUM)/SV_POSITION_QUANTUM];
        }
        
        chromosomePLengths = new int[N_CHROMOSOMES];
        br = new BufferedReader(new FileReader(cytobandFile));
        str=br.readLine();
        while (str!=null) {
            tokens=str.split("\t");
            i=str2chromosomeID(tokens[0]);
            if (i==-1) {
                str=br.readLine();
                continue;
            }
            length=Integer.parseInt(tokens[2]);
            if (tokens[3].charAt(0)=='p' && length>chromosomePLengths[i]) chromosomePLengths[i]=length;
            str=br.readLine();
        }
        br.close();
    }
    
    
    /**
     * @return -1 if $str$ does not represent a non-mitochondrial chromosome;
     * the 0-based integer ID of a chromosome otherwise.
     */
    private static final int str2chromosomeID(String str) {
        char c;
        
        if (str.length()>3) {
            if (!str.substring(0,3).equalsIgnoreCase("chr")) return -1;
            if (str.length()>5) return -1;
            c=str.charAt(3);
        }
        else c=str.charAt(0);
        if (c=='M' || c=='m') return -1;
        if (c=='X' || c=='x') return N_CHROMOSOMES-2;
        else if (c=='Y' || c=='y') return N_CHROMOSOMES-1;
        else return Integer.parseInt(str.length()>3?str.substring(3):str)-1;
    }
    
    
    /**
     * Stretches the P-arm of every autosome onto the P-arm of chr1 and 
     * cumulates SV counts from $svs_per_position$ proportionally in every bin; 
     * does the same for the Q-arm; then, normalizes the total counts over the
     * entire chr1. Using only autosomes comes from the gnomAD-SV paper.
     *
     * Remark: contrary to the gnomAD-SV paper, we do not exclude bins that
     * overlap centromeres.
     */
    private static final void buildStatisticsPosition() {
        int i, j, k, h;
        int first, last, firstP, lastP, firstQ, lastQ, length;
        int projectedFirst, projectedLast, projectedFirstBin, projectedLastBin;
        double ratioP, ratioQ, massP, massQ, projectedMassPerBp, sum;
        
        statistics_position = new double[SV_TYPES.length][svs_per_position[0][0].length];
        for (i=0; i<SV_TYPES.length; i++) {
            System.arraycopy(svs_per_position[i][0],0,statistics_position[i],0,svs_per_position[i][0].length);
            for (j=1; j<N_AUTOSOMES; j++) {
                ratioP=((double)chromosomePLengths[0])/chromosomePLengths[j];
                ratioQ=((double)(chromosomeLengths[0]-chromosomePLengths[0]))/(chromosomeLengths[j]-chromosomePLengths[j]);
                length=svs_per_position[i][j].length;
                for (k=0; k<length; k++) {
                    first=k*SV_POSITION_QUANTUM; 
                    last=Math.min((k+1)*SV_POSITION_QUANTUM,chromosomeLengths[j])-1;
                    if (last<chromosomePLengths[j]) { 
                        firstP=first; lastP=last; firstQ=-1; lastQ=-1; 
                        massP=svs_per_position[i][j][k]; massQ=0.0; 
                    }
                    else if (first<chromosomePLengths[j] && last>=chromosomePLengths[j]) {
                        firstP=first; lastP=chromosomePLengths[j]-1;
                        firstQ=chromosomePLengths[j]; lastQ=last;
                        massP=((double)((lastP-firstP+1)*svs_per_position[i][j][k]))/SV_POSITION_QUANTUM;
                        massQ=((double)((lastQ-firstQ+1)*svs_per_position[i][j][k]))/SV_POSITION_QUANTUM;
                    }
                    else {
                        firstP=-1; lastP=-1; firstQ=first; lastQ=last; 
                        massP=0.0; massQ=svs_per_position[i][j][k];
                    }
                    if (firstP!=-1) {
                        projectedFirst=(int)(firstP*ratioP); projectedLast=(int)(lastP*ratioP);
                        if (projectedFirst<chromosomePLengths[0]) {
                            projectedMassPerBp=massP/(projectedLast-projectedFirst+1);
                            projectedLast=Math.min(projectedLast,chromosomePLengths[0]-1);
                            projectedFirstBin=projectedFirst/SV_POSITION_QUANTUM; projectedLastBin=projectedLast/SV_POSITION_QUANTUM;
                            if (projectedLastBin==projectedFirstBin) statistics_position[i][projectedFirstBin]+=projectedMassPerBp*(projectedLast-projectedFirst+1);
                            else {
                                statistics_position[i][projectedFirstBin]+=projectedMassPerBp*((projectedFirstBin+1)*SV_POSITION_QUANTUM-projectedFirst);
                                for (h=projectedFirstBin+1; h<projectedLastBin; h++) statistics_position[i][h]+=projectedMassPerBp*SV_POSITION_QUANTUM;
                                statistics_position[i][projectedLastBin]+=projectedMassPerBp*(projectedLast+1-projectedLastBin*SV_POSITION_QUANTUM);
                            }
                        }
                    }
                    if (firstQ!=-1) {
                        projectedFirst=chromosomePLengths[0]+(int)((firstQ-chromosomePLengths[j])*ratioQ); 
                        projectedLast=chromosomePLengths[0]+(int)((lastQ-chromosomePLengths[j])*ratioQ);
                        projectedMassPerBp=massQ/(projectedLast-projectedFirst+1);
                        projectedLast=Math.min(projectedLast,chromosomeLengths[0]-1);
                        projectedFirstBin=projectedFirst/SV_POSITION_QUANTUM; projectedLastBin=projectedLast/SV_POSITION_QUANTUM;
                        if (projectedLastBin==projectedFirstBin) statistics_position[i][projectedFirstBin]+=projectedMassPerBp*(projectedLast-projectedFirst+1);
                        else {
                            statistics_position[i][projectedFirstBin]+=projectedMassPerBp*((projectedFirstBin+1)*SV_POSITION_QUANTUM-projectedFirst);
                            for (h=projectedFirstBin+1; h<projectedLastBin; h++) statistics_position[i][h]+=projectedMassPerBp*SV_POSITION_QUANTUM;
                            statistics_position[i][projectedLastBin]+=projectedMassPerBp*(projectedLast+1-projectedLastBin*SV_POSITION_QUANTUM);
                        }
                    }
                }
            }
        }
        for (i=0; i<SV_TYPES.length; i++) {
            sum=0.0; length=statistics_position[i].length;
            for (j=0; j<length; j++) sum+=statistics_position[i][j];
            for (j=0; j<length; j++) statistics_position[i][j]/=sum;
        }
    }

    
	
	
	// ----------------------- STATISTICS PROCEDURES ---------------------------
	
	/**
	 * (SV type, reference interval type) -> Prob(SV length)
	 */
	public static double[][][] statistics_length;
	
	/**
	 * SV type -> Prob(reference interval type)
	 */
	public static double[][] statistics_intervalType;
	
	/**
	 * Prob(SV type)
	 */
	public static double[] statistics_svType;
	
	/**
	 * (Population, SV type) -> P(frequency in the population)
	 */
	public static double[][][] statistics_alleleFrequency;
	private static int[][] statistics_lastAlleleFrequency;
	public static double[][][] statistics_alleleFrequency_minmax;
    
	/**
	 * The only SVs supported for statistics.
	 */
	public static final String[] SV_TYPES = {"DEL","DUP","INS","INV"};
	public static final int MIN_SV_LENGTH = 50;
	private static final int MAX_SV_LENGTH = 1000000;
	public static final int SV_LENGTH_QUANTUM = 50;
	public static final int ALLELE_FREQUENCY_NBINS = 10000;  // Arbitrary
	public static final String[] AF_NAMES = new String[] {"ALL", "AFR", "AMR", "EAS", "EUR", "OTH"};
	public static final String[] AF_FILE_NAMES = new String[] {"statistics_alleleFrequency_ALL", "statistics_alleleFrequency_AFR", "statistics_alleleFrequency_AMR", "statistics_alleleFrequency_EAS", "statistics_alleleFrequency_EUR", "statistics_alleleFrequency_OTH"};
	public static final String[] ALLELE_FREQUENCY_STR = new String[] {";AF=", ";AFR_AF=", ";AMR_AF=", ";EAS_AF=", ";EUR_AF=", ";OTH_AF="};
	
	/**
	 * Types of repeat contexts: 
	 * 0: non-satellite repeat;
	 * 1: satellite repeat;
	 * 2: both of the above;
	 * 3: none of the above, but segmental duplication;
	 * 4: none of the above.
	 */
	public static final int N_INTERVAL_TYPES = 5;
	
	/**
	 * Other constants defining the format of gnomAD-SV
	 */
	public static final String SVLEN_STR = "SVLEN=";
	public static final int SVLEN_STR_LENGTH = SVLEN_STR.length();
	public static final String SVTYPE_STR = "SVTYPE=";
	public static final int SVTYPE_STR_LENGTH = SVTYPE_STR.length();
	public static final String INFO_SEPARATOR = ";";
	public static final String PASS_STR = "PASS";
	public static final String MULTIALLELIC_STR = "MULTIALLELIC";
	public static final String EVIDENCE_STR = "EVIDENCE";
	public static final int EVIDENCE_STR_LENGTH = EVIDENCE_STR.length();
	public static final String SR_STR = "SR";
	public static final String PE_STR = "PE";
	public static final int[] ALLELE_FREQUENCY_STR_LENGTH = new int[] {4,8,8,8,8,8};
	public static final int SMOOTHING_RADIUS = 4;  // In bins. Arbitrary.
	public static final int SMOOTHING_WINDOW = 1+((SMOOTHING_RADIUS)<<1);
	public static final double MIN_FREQUENCY = 0.01;  // Used just for statistics
	
	
	/**
	 * Computes global variables $statistics_*$, i.e. the empirical probability 
	 * that an SV of a specific type in $SV_TYPES$ and of a specific length,
	 * starts inside a repeat interval that is fully satellite, fully non-
	 * satellite, both, or none.
	 *
	 * Remark: statistics are computed over all chromosomes.
	 *
	 * Remark: length counts are smoothed by a simple moving average, as
     * suggested by Ryan Collins. This is done to reduce peaks in DEL and DUP 
	 * distributions that occur at multiples of 1kb, starting from 5kb. Such 
	 * peaks are due to the read-depth CNV caller used in gnomAD-SV, which
	 * generates calls in fixed increments of approx. 1kb, and to the fact that 
	 * gnomAD-SV does not include calls exclusively made from read depth
	 * evidence below 5kb.
     *
     * Remark: position counts are smoothed by a simple moving average, as
     * described in the gnomAD-SV paper.
	 *
	 * Remark: we don't further decompose non-satellite repeats into transposon
	 * classes. This is done for simplicity.
	 *
	 * @param gnomadVCFfile contigs in the file are assumed to be sorted as
	 * follows (contigs in "other" are not loaded):
	 *
	 * chr1,chr10,...,chr19,  chr2,...,chr22,  chr3,...,chr9,  chrX,chrY, other.
	 *
	 * @param maxDistance_* a position is assigned to a repeat interval if it is 
	 * at most this far from it; there is a different threshold for each type of
	 * SV evidence, as suggested by Ryan Collins;
	 * @param onlyPass uses only records whose filter field equals PASS or
	 * MULTIALLELIC.
	 */
	public static final void collectStatistics(String gnomadVCFfile, int maxDistance_SR, int maxDistance_PE, int maxDistance_other, int maxDistance_segdups, boolean onlyPass, boolean verbose) throws IOException {
		boolean hasSR, hasPE, isFrequent;
		int i, j, k, h, p, q, x;
		int length, svType, intervalType, position, maxDistance, n_chr1, nFrequent_chr1;
		double sum, min, max, quantum, boundary, alleleFrequency;
		String str, info, chromosome, evidence;
		BufferedReader br;
		double[] smoothedArray;
		String[] tokens;
		
		statistics_length = new double[SV_TYPES.length][N_INTERVAL_TYPES][1+(MAX_SV_LENGTH-MIN_SV_LENGTH+1)/SV_LENGTH_QUANTUM];
		smoothedArray = new double[statistics_length[0][0].length];
		statistics_intervalType = new double[SV_TYPES.length][N_INTERVAL_TYPES];
		statistics_svType = new double[SV_TYPES.length];
		statistics_alleleFrequency = new double[ALLELE_FREQUENCY_STR.length][SV_TYPES.length][ALLELE_FREQUENCY_NBINS];
		statistics_lastAlleleFrequency = new int[ALLELE_FREQUENCY_STR.length][SV_TYPES.length];
		for (i=0; i<statistics_lastAlleleFrequency.length; i++) {
			for (j=0; j<statistics_lastAlleleFrequency[i].length; j++) statistics_lastAlleleFrequency[i][j]=-1;
		}
		statistics_alleleFrequency_minmax = new double[ALLELE_FREQUENCY_STR.length][SV_TYPES.length][2];
		
		// Collecting counts
		br = new BufferedReader(new FileReader(gnomadVCFfile));
		str=br.readLine(); k=0; x=0; n_chr1=0; nFrequent_chr1=0;
		while (str!=null) {
			if (str.charAt(0)=='#') {
				str=br.readLine();
				continue;
			}
			tokens=str.split("\t");
			if (tokens[0].length()>5 || tokens[0].equalsIgnoreCase("chrM")) break;
			if (onlyPass && (!tokens[6].equalsIgnoreCase(PASS_STR) && !tokens[6].equalsIgnoreCase(MULTIALLELIC_STR))) {
				str=br.readLine();
				continue;
			}
			info=tokens[7];
			p=info.indexOf(SVLEN_STR);
			if (p<0) {
				str=br.readLine();
				continue;
			}
			q=info.indexOf(INFO_SEPARATOR,p+SVLEN_STR_LENGTH);
			length=Integer.parseInt(info.substring(p+SVLEN_STR_LENGTH,q));
			if (length<0) length=-length;
			if (length<MIN_SV_LENGTH) {
				str=br.readLine();
				continue;
			}
			p=info.indexOf(SVTYPE_STR);
			if (p<0) {
				str=br.readLine();
				continue;
			}
			q=info.indexOf(INFO_SEPARATOR,p+SVTYPE_STR_LENGTH);
			svType=Arrays.binarySearch(SV_TYPES,info.substring(p+SVTYPE_STR_LENGTH,q));
			if (svType<0) {
				str=br.readLine();
				continue;
			}
			chromosome="chr"+tokens[0];
			if (tokens[0].equalsIgnoreCase("1")) n_chr1++;
			position=Integer.parseInt(tokens[1]);  // VCF is one-based, but this position is the one just before the variant.
			p=info.indexOf(EVIDENCE_STR);
			q=info.indexOf(INFO_SEPARATOR,p+EVIDENCE_STR_LENGTH);
			evidence=info.substring(p+EVIDENCE_STR_LENGTH,q);
			hasSR=evidence.indexOf(SR_STR)>=0;
			hasPE=evidence.indexOf(PE_STR)>=0;
			if (hasSR) maxDistance=maxDistance_SR;
			else if (hasPE) maxDistance=maxDistance_PE;
			else maxDistance=maxDistance_other;
			while (k<=lastRepeatInterval && (repeatIntervals_chromosomes[k].compareTo(chromosome)<0 || repeatIntervals[k][2]+maxDistance<position)) {
				k++;
				continue;
			}
			if (k>lastRepeatInterval || repeatIntervals_chromosomes[k].compareTo(chromosome)>0 || repeatIntervals[k][1]>position+maxDistance) intervalType=4;
			else {
				intervalType=repeatIntervals[k][0];
				for (h=k-1; h>=0; h--) {
					if (repeatIntervals[h][2]+maxDistance<position) break;
					intervalType=ADD_TYPE[intervalType][repeatIntervals[h][0]];
				}
				for (h=k+1; h<=lastRepeatInterval; h++) {
					if (repeatIntervals[h][1]>position+maxDistance) break;
					intervalType=ADD_TYPE[intervalType][repeatIntervals[h][0]];
				}
			}
			if (intervalType==4) {
				while (x<=lastSegdup && (segdups_chromosomeIDs[x].compareTo(chromosome)<0 || segdups[x][1]+maxDistance_segdups<position)) {
					x++;
					continue;
				}
				if (x>lastSegdup || segdups_chromosomeIDs[x].compareTo(chromosome)>0 || segdups[x][0]>position+maxDistance_segdups) intervalType=4;
				else intervalType=3;
			}
			length=(length-MIN_SV_LENGTH)/SV_LENGTH_QUANTUM;			
			statistics_length[svType][intervalType][length>=statistics_length[svType][intervalType].length?statistics_length[svType][intervalType].length-1:length]++;
			if (!tokens[0].equalsIgnoreCase("X") && !tokens[0].equalsIgnoreCase("Y")) {
				// We don't use allele frequencies of X,Y since they can be
				// greater than one.
				isFrequent=false;
				for (i=0; i<ALLELE_FREQUENCY_STR.length; i++) {
					p=info.indexOf(ALLELE_FREQUENCY_STR[i]);
					q=info.indexOf(INFO_SEPARATOR,p+ALLELE_FREQUENCY_STR_LENGTH[i]);
					alleleFrequency=Double.parseDouble(info.substring(p+ALLELE_FREQUENCY_STR_LENGTH[i],q));
					if (tokens[0].equalsIgnoreCase("1") && alleleFrequency>=MIN_FREQUENCY) isFrequent=true;
					statistics_lastAlleleFrequency[i][svType]++;
					if (statistics_lastAlleleFrequency[i][svType]==statistics_alleleFrequency[i][svType].length) {
						double[] newArray = new double[statistics_alleleFrequency[i][svType].length<<1];
						System.arraycopy(statistics_alleleFrequency[i][svType],0,newArray,0,statistics_alleleFrequency[i][svType].length);
						statistics_alleleFrequency[i][svType]=newArray;
					}
					statistics_alleleFrequency[i][svType][statistics_lastAlleleFrequency[i][svType]]=alleleFrequency;
				}
				if (isFrequent) nFrequent_chr1++;
			}
            svs_per_position[svType][str2chromosomeID(tokens[0])][position/SV_POSITION_QUANTUM]++;
			str=br.readLine();
		}
		br.close();
		System.err.println("Total SVs in chr1: "+n_chr1+". With allele frequency >="+MIN_FREQUENCY+": "+nFrequent_chr1);

		// Smoothing length counts
		for (i=0; i<SV_TYPES.length; i++) {
			for (j=0; j<N_INTERVAL_TYPES; j++) {
				Arrays.fill(smoothedArray,0.0);
				sum=0.0;
				for (k=0; k<SMOOTHING_WINDOW; k++) sum+=statistics_length[i][j][k];
				for (k=SMOOTHING_RADIUS; k<=statistics_length[i][j].length-SMOOTHING_RADIUS-2; k++) {
					smoothedArray[k]=sum/SMOOTHING_WINDOW;
					sum+=statistics_length[i][j][k+SMOOTHING_RADIUS+1]-statistics_length[i][j][k-SMOOTHING_RADIUS];
				}
				smoothedArray[statistics_length[i][j].length-SMOOTHING_RADIUS-1]=sum/SMOOTHING_WINDOW;
				System.arraycopy(smoothedArray,SMOOTHING_RADIUS,statistics_length[i][j],SMOOTHING_RADIUS,statistics_length[i][j].length-((SMOOTHING_RADIUS)<<1));
			}
		}
        
        // Smoothing position counts
        smoothedArray = new double[svs_per_position[0][0].length];
        for (i=0; i<SV_TYPES.length; i++) {
            for (j=0; j<N_CHROMOSOMES; j++) {
				Arrays.fill(smoothedArray,0,svs_per_position[i][j].length,0.0);
				sum=0.0;
				for (k=0; k<SMOOTHING_WINDOW_POSITIONS; k++) sum+=svs_per_position[i][j][k];
				for (k=SMOOTHING_RADIUS_POSITIONS; k<=svs_per_position[i][j].length-SMOOTHING_RADIUS_POSITIONS-2; k++) {
					smoothedArray[k]=sum/SMOOTHING_WINDOW_POSITIONS;
					sum+=svs_per_position[i][j][k+SMOOTHING_RADIUS_POSITIONS+1]-svs_per_position[i][j][k-SMOOTHING_RADIUS_POSITIONS];
				}
				smoothedArray[svs_per_position[i][j].length-SMOOTHING_RADIUS_POSITIONS-1]=sum/SMOOTHING_WINDOW_POSITIONS;
				System.arraycopy(smoothedArray,SMOOTHING_RADIUS_POSITIONS,svs_per_position[i][j],SMOOTHING_RADIUS_POSITIONS,svs_per_position[i][j].length-((SMOOTHING_RADIUS_POSITIONS)<<1));
            }
        }
		
		// Transforming counts into empirical probabilities
		for (i=0; i<SV_TYPES.length; i++) {
			for (j=0; j<N_INTERVAL_TYPES; j++) {
				sum=0.0; length=statistics_length[i][j].length;
				for (k=0; k<length; k++) sum+=statistics_length[i][j][k];
				for (k=0; k<length; k++) statistics_length[i][j][k]/=sum;
				statistics_intervalType[i][j]=sum;
			}
		}
		for (i=0; i<SV_TYPES.length; i++) {
			sum=0.0;
			for (j=0; j<N_INTERVAL_TYPES; j++) sum+=statistics_intervalType[i][j];
			for (j=0; j<N_INTERVAL_TYPES; j++) statistics_intervalType[i][j]/=sum;
			statistics_svType[i]=sum;
		}
		sum=0.0;
		for (i=0; i<SV_TYPES.length; i++) sum+=statistics_svType[i];
		for (i=0; i<SV_TYPES.length; i++) statistics_svType[i]/=sum;
		for (h=0; h<ALLELE_FREQUENCY_STR.length; h++) {
            length=SV_TYPES.length;
			for (i=0; i<length; i++) {
				if (statistics_lastAlleleFrequency[h][i]>0) Arrays.sort(statistics_alleleFrequency[h][i],0,statistics_lastAlleleFrequency[h][i]+1);
				min=statistics_alleleFrequency[h][i][0];
				max=statistics_alleleFrequency[h][i][statistics_lastAlleleFrequency[h][i]];
				quantum=(max-min)/(ALLELE_FREQUENCY_NBINS-1);
				statistics_alleleFrequency_minmax[h][i][0]=min;
				statistics_alleleFrequency_minmax[h][i][1]=max;
				k=0;
				for (j=0; j<ALLELE_FREQUENCY_NBINS; j++) {
					boundary=min+j*quantum; sum=0.0;
					while (k<=statistics_lastAlleleFrequency[h][i] && statistics_alleleFrequency[h][i][k]<boundary) { sum++; k++; }
					statistics_alleleFrequency[h][i][j]=sum/(statistics_lastAlleleFrequency[h][i]+1);
				}
			}
		}
        buildStatisticsPosition();
	}
	
	
	private static final void makeStatisticsCumulative() throws IOException {
		int i, j, k;
		int lengthI, lengthJ, lengthK;
		
		lengthI=statistics_svType.length;
		for (i=1; i<lengthI; i++) statistics_svType[i]+=statistics_svType[i-1];
		lengthI=statistics_intervalType.length;
		for (i=0; i<lengthI; i++) {
			lengthJ=statistics_intervalType[i].length;
			for (j=1; j<lengthJ; j++) statistics_intervalType[i][j]+=statistics_intervalType[i][j-1];
		}
		lengthI=statistics_length.length;
		for (i=0; i<lengthI; i++) {
			lengthJ=statistics_length[i].length;
			for (j=0; j<lengthJ; j++) {
				lengthK=statistics_length[i][j].length;
				for (k=1; k<lengthK; k++) statistics_length[i][j][k]+=statistics_length[i][j][k-1];
			}
		}
		lengthI=statistics_alleleFrequency.length;
		for (i=0; i<lengthI; i++) {
			lengthJ=statistics_alleleFrequency[i].length;
			for (j=0; j<lengthJ; j++) {
				lengthK=statistics_alleleFrequency[i][j].length;
				for (k=1; k<lengthK; k++) statistics_alleleFrequency[i][j][k]+=statistics_alleleFrequency[i][j][k-1];
			}
		}
		lengthI=statistics_position.length;
		for (i=0; i<lengthI; i++) {
			lengthJ=statistics_position[i].length;
			for (j=1; j<lengthJ; j++) statistics_position[i][j]+=statistics_position[i][j-1];
		}
	}
	
	
	private static final void serializeStatistics(String outputDir) throws IOException {
		int i, j, k;
		int lengthI, lengthJ, lengthK;
		BufferedWriter bw;
		
		// statistics_svType
		bw = new BufferedWriter(new FileWriter(outputDir+"/statistics_svType.txt"));
		lengthI=statistics_svType.length;
		bw.write(lengthI+"\n");
		for (i=0; i<lengthI; i++) bw.write(statistics_svType[i]+",");
		bw.newLine();
		bw.close();
		
		// statistics_intervalType
		bw = new BufferedWriter(new FileWriter(outputDir+"/statistics_intervalType.txt"));
		lengthI=statistics_intervalType.length;
		lengthJ=statistics_intervalType[0].length;
		bw.write(lengthI+","+lengthJ+"\n");
		for (i=0; i<lengthI; i++) {
			for (j=0; j<lengthJ; j++) bw.write(statistics_intervalType[i][j]+",");
			bw.newLine();
		}
		bw.close();
		
		// statistics_length
		bw = new BufferedWriter(new FileWriter(outputDir+"/statistics_length.txt"));
		lengthI=statistics_length.length;
		lengthJ=statistics_length[0].length;
		lengthK=statistics_length[0][0].length;
		bw.write(lengthI+","+lengthJ+","+lengthK+"\n");
		for (i=0; i<lengthI; i++) {
			for (j=0; j<lengthJ; j++) {
				for (k=0; k<lengthK; k++) bw.write(statistics_length[i][j][k]+",");
				bw.newLine();
			}
		}
		bw.close();
		
		// statistics_alleleFrequency
		for (k=0; k<statistics_alleleFrequency.length; k++) {
			bw = new BufferedWriter(new FileWriter(outputDir+"/"+AF_FILE_NAMES[k]+".txt"));
			lengthI=statistics_alleleFrequency[k].length;
			lengthJ=ALLELE_FREQUENCY_NBINS;
			bw.write(lengthI+","+lengthJ+"\n");
			for (i=0; i<lengthI; i++) {
				bw.write(statistics_alleleFrequency_minmax[k][i][0]+","+statistics_alleleFrequency_minmax[k][i][1]+",");
				for (j=0; j<lengthJ; j++) bw.write(statistics_alleleFrequency[k][i][j]+",");
				bw.newLine();
			}
			bw.close();
		}
        
		// statistics_position
		bw = new BufferedWriter(new FileWriter(outputDir+"/statistics_position.txt"));
		lengthI=statistics_position.length;
		lengthJ=statistics_position[0].length;
		bw.write(lengthI+","+lengthJ+"\n");
		for (i=0; i<lengthI; i++) {
			for (j=0; j<lengthJ; j++) bw.write(statistics_position[i][j]+",");
			bw.newLine();
		}
		bw.close();
	}
	
	
	public static final void deserializeStatistics(String inputDir) throws IOException {
		int i, j, k;
		int lengthI, lengthJ, lengthK;
		String str;
		BufferedReader br;
		String[] tokens;
		
		// statistics_svType
		br = new BufferedReader(new FileReader(inputDir+"/statistics_svType.txt"));
		str=br.readLine();
		lengthI=Integer.parseInt(str);
		statistics_svType = new double[lengthI];
		str=br.readLine();
		tokens=str.split(",");
		for (i=0; i<lengthI; i++) statistics_svType[i]=Double.parseDouble(tokens[i]);
		br.close();
		
		// statistics_intervalType
		br = new BufferedReader(new FileReader(inputDir+"/statistics_intervalType.txt"));
		str=br.readLine();
		tokens=str.split(",");
		lengthI=Integer.parseInt(tokens[0]);
		lengthJ=Integer.parseInt(tokens[1]);
		statistics_intervalType = new double[lengthI][lengthJ];
		for (i=0; i<lengthI; i++) {
			str=br.readLine();
			tokens=str.split(",");
			for (j=0; j<lengthJ; j++) statistics_intervalType[i][j]=Double.parseDouble(tokens[j]);
		}
		br.close();
		
		// statistics_length
		br = new BufferedReader(new FileReader(inputDir+"/statistics_length.txt"));
		str=br.readLine();
		tokens=str.split(",");
		lengthI=Integer.parseInt(tokens[0]);
		lengthJ=Integer.parseInt(tokens[1]);
		lengthK=Integer.parseInt(tokens[2]);
		statistics_length = new double[lengthI][lengthJ][lengthK];
		for (i=0; i<lengthI; i++) {
			for (j=0; j<lengthJ; j++) {
				str=br.readLine();
				tokens=str.split(",");
				for (k=0; k<lengthK; k++) statistics_length[i][j][k]=Double.parseDouble(tokens[k]);
			}
		}
		br.close();
		
		// statistics_alleleFrequency
		statistics_alleleFrequency = new double[AF_FILE_NAMES.length][0][0];
		statistics_alleleFrequency_minmax = new double[AF_FILE_NAMES.length][0][0];
		for (k=0; k<AF_FILE_NAMES.length; k++) {
			br = new BufferedReader(new FileReader(inputDir+"/"+AF_FILE_NAMES[k]+".txt"));
			str=br.readLine();
			tokens=str.split(",");
			lengthI=Integer.parseInt(tokens[0]);
			lengthJ=Integer.parseInt(tokens[1]);
			statistics_alleleFrequency[k] = new double[lengthI][lengthJ];
			statistics_alleleFrequency_minmax[k] = new double[lengthI][2];
			for (i=0; i<lengthI; i++) {
				str=br.readLine();
				tokens=str.split(",");
				statistics_alleleFrequency_minmax[k][i][0]=Double.parseDouble(tokens[0]);
				statistics_alleleFrequency_minmax[k][i][1]=Double.parseDouble(tokens[1]);
				for (j=0; j<lengthJ; j++) statistics_alleleFrequency[k][i][j]=Double.parseDouble(tokens[2+j]);
			}
			br.close();
		}
        
		// statistics_position
		br = new BufferedReader(new FileReader(inputDir+"/statistics_position.txt"));
		str=br.readLine();
		tokens=str.split(",");
		lengthI=Integer.parseInt(tokens[0]);
		lengthJ=Integer.parseInt(tokens[1]);
		statistics_position = new double[lengthI][lengthJ];
		for (i=0; i<lengthI; i++) {
			str=br.readLine();
			tokens=str.split(",");
			for (j=0; j<lengthJ; j++) statistics_position[i][j]=Double.parseDouble(tokens[j]);
		}
		br.close();
	}
	
	
	/**
	 * For every population and every allele frequency bin, the procedure prints
	 * the number of SV breakpoints that occur inside or near every interval 
	 * type.
	 * 
	 * Remark: the procedure uses the allele frequencies of X and Y as well,
	 * even though they can be greater than one.
	 *
	 * @param maxDistance_* a position is assigned to a repeat interval if it is 
	 * at most this far from it; there is a different threshold for each type of
	 * SV evidence, as suggested by Ryan Collins;
	 * @param onlyPass uses only records whose filter field equals PASS or
	 * MULTIALLELIC.
	 */
	private static final void alleleFrequencyVsRepeatContext(String gnomadVCFfile, int maxDistance_SR, int maxDistance_PE, int maxDistance_other, int maxDistance_segdups, boolean onlyPass, String outputDir) throws IOException {
		final double FREQ_MIN = -4.5;  // log10. Arbitrary.
		final double FREQ_MAX = 0.0;
		final int N_BINS = 100;  // Arbitrary
		final double FREQ_QUANTUM = (FREQ_MAX-FREQ_MIN)/N_BINS;
		boolean hasSR, hasPE;
		int i, j, k, h, p, q, x;
		int length, svType, intervalType, position, maxDistance;
		double alleleFrequency;
		String str, info, chromosome, evidence;
		BufferedReader br;
		BufferedWriter bw;
		String[] tokens;
		double[][][] frequencies;
		
		frequencies = new double[AF_NAMES.length][N_BINS+1][5];  // Arbitrary
		for (i=0; i<frequencies.length; i++) {
			for (j=0; j<frequencies[i].length; j++) Arrays.fill(frequencies[i][j],0);
		}
		br = new BufferedReader(new FileReader(gnomadVCFfile));
		str=br.readLine(); k=0; x=0;
		while (str!=null) {
			if (str.charAt(0)=='#') {
				str=br.readLine();
				continue;
			}
			tokens=str.split("\t");
			if (tokens[0].length()>5 || tokens[0].equalsIgnoreCase("chrM")) break;
			if (onlyPass && (!tokens[6].equalsIgnoreCase(PASS_STR) && !tokens[6].equalsIgnoreCase(MULTIALLELIC_STR))) {
				str=br.readLine();
				continue;
			}
			info=tokens[7];
			p=info.indexOf(SVLEN_STR);
			if (p<0) {
				str=br.readLine();
				continue;
			}
			q=info.indexOf(INFO_SEPARATOR,p+SVLEN_STR_LENGTH);
			length=Integer.parseInt(info.substring(p+SVLEN_STR_LENGTH,q));
			if (length<0) length=-length;
			p=info.indexOf(SVTYPE_STR);
			if (p<0) {
				str=br.readLine();
				continue;
			}
			q=info.indexOf(INFO_SEPARATOR,p+SVTYPE_STR_LENGTH);
			svType=Arrays.binarySearch(SV_TYPES,info.substring(p+SVTYPE_STR_LENGTH,q));
			if (svType<0) {
				str=br.readLine();
				continue;
			}
			chromosome="chr"+tokens[0];
			position=Integer.parseInt(tokens[1]);  // VCF is one-based, but this position is the one just before the variant.
			p=info.indexOf(EVIDENCE_STR);
			q=info.indexOf(INFO_SEPARATOR,p+EVIDENCE_STR_LENGTH);
			evidence=info.substring(p+EVIDENCE_STR_LENGTH,q);
			hasSR=evidence.indexOf(SR_STR)>=0;
			hasPE=evidence.indexOf(PE_STR)>=0;
			if (hasSR) maxDistance=maxDistance_SR;
			else if (hasPE) maxDistance=maxDistance_PE;
			else maxDistance=maxDistance_other;
			
			// Computing the interval type of the first breakpoint
			while (k<=lastRepeatInterval && (repeatIntervals_chromosomes[k].compareTo(chromosome)<0 || repeatIntervals[k][2]+maxDistance<position)) {
				k++;
				continue;
			}
			if (k>lastRepeatInterval || repeatIntervals_chromosomes[k].compareTo(chromosome)>0 || repeatIntervals[k][1]>position+maxDistance) intervalType=4;
			else {
				intervalType=repeatIntervals[k][0];
				for (h=k-1; h>=0; h--) {
					if (repeatIntervals[h][2]+maxDistance<position) break;
					intervalType=ADD_TYPE[intervalType][repeatIntervals[h][0]];
				}
				for (h=k+1; h<=lastRepeatInterval; h++) {
					if (repeatIntervals[h][1]>position+maxDistance) break;
					intervalType=ADD_TYPE[intervalType][repeatIntervals[h][0]];
				}
			}
			if (intervalType==4) {
				while (x<=lastSegdup && (segdups_chromosomeIDs[x].compareTo(chromosome)<0 || segdups[x][1]+maxDistance_segdups<position)) {
					x++;
					continue;
				}
				if (x>lastSegdup || segdups_chromosomeIDs[x].compareTo(chromosome)>0 || segdups[x][0]>position+maxDistance_segdups) intervalType=4;
				else intervalType=3;
			}
			for (i=0; i<ALLELE_FREQUENCY_STR.length; i++) {
				p=info.indexOf(ALLELE_FREQUENCY_STR[i]);
				q=info.indexOf(INFO_SEPARATOR,p+ALLELE_FREQUENCY_STR_LENGTH[i]);
				alleleFrequency=Math.log10(Double.parseDouble(info.substring(p+ALLELE_FREQUENCY_STR_LENGTH[i],q)));
				if (alleleFrequency>=FREQ_MIN && alleleFrequency<=FREQ_MAX) frequencies[i][(int)((alleleFrequency-FREQ_MIN)/FREQ_QUANTUM)][intervalType]++;
			}
			
			// Computing the interval type of the second breakpoint, if any.
			if (svType==2) {
				str=br.readLine(); 
				continue;
			}
			intervalType=getIntervalType(chromosome,position+length-1,maxDistance,maxDistance_segdups);
			for (i=0; i<ALLELE_FREQUENCY_STR.length; i++) {
				p=info.indexOf(ALLELE_FREQUENCY_STR[i]);
				q=info.indexOf(INFO_SEPARATOR,p+ALLELE_FREQUENCY_STR_LENGTH[i]);
				alleleFrequency=Math.log10(Double.parseDouble(info.substring(p+ALLELE_FREQUENCY_STR_LENGTH[i],q)));
				if (alleleFrequency>=FREQ_MIN && alleleFrequency<=FREQ_MAX) frequencies[i][(int)((alleleFrequency-FREQ_MIN)/FREQ_QUANTUM)][intervalType]++;
			}
			str=br.readLine();
		}
		br.close();
		
		// Writing to disk
		for (i=0; i<frequencies.length; i++) {
			bw = new BufferedWriter(new FileWriter(outputDir+"/breakpoints_"+AF_NAMES[i]+".txt"));
			for (j=0; j<frequencies[i].length; j++) {
				bw.write((FREQ_MIN+j*FREQ_QUANTUM)+",");
				for (k=0; k<frequencies[i][j].length; k++) bw.write(frequencies[i][j][k]+",");
				bw.newLine();
			}
			bw.close();
		}
	}
	
	
	/**
	 * Remark: the procedure is doing a horrible linear scan, just because a 
	 * quick implementation was needed for answering to a reviewer. Sorry, 
	 * reader.
	 */
	private static final int getIntervalType(String chr, int position, int maxDistance, int maxDistance_segdups) {
		int i, j;
		int intervalType, chrID;
		char c;
		
		c=chr.charAt(3);
		if (c=='X') chrID=22;
		else if (c=='Y') chrID=23;
		else chrID = Integer.parseInt(chr.substring(3))-1;
		intervalType=4;
		for (i=chromosome2firstInterval[chrID]; i<=lastRepeatInterval; i++) {
			j=chr.compareTo(repeatIntervals_chromosomes[i]);
			if (j<0) continue;
 			else if (j>0) break;
			if (repeatIntervals[i][2]+maxDistance<position) continue;
			else if (repeatIntervals[i][1]>position+maxDistance) break;
			else if (intervalType==4) intervalType=repeatIntervals[i][0];
			else intervalType=ADD_TYPE[intervalType][repeatIntervals[i][0]];
		}
		if (intervalType==4) {
			for (i=chromosome2firstSegdup[chrID]; i<=lastSegdup; i++) {
				j=chr.compareTo(segdups_chromosomeIDs[i]);
				if (j<0) continue;
				else if (j>0) break;
				if (segdups[i][1]+maxDistance_segdups<position) continue;
				else if (segdups[i][0]>position+maxDistance_segdups) break;
				else {
					intervalType=3;
					break;
				}
			}
		}
		return intervalType;
	}

}