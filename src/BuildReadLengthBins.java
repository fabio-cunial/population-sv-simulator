import java.util.Arrays;
import java.util.Random;
import java.io.*;


/**
 * Given a set of FASTQ files, partitions all their reads into bins, permutes 
 * every bin randomly, and stores every bin to a separate file. Every bin file
 * contains one line per read, which is the concatenation of the four rows of a
 * FASTQ read.
 */
public class BuildReadLengthBins {

    public static final String BIN_FILE_SEPARATOR = "__________";


	public static void main(String[] args) throws IOException {
		final String INPUT_FILES_LIST = args[0];
        final int BIN_LENGTH = Integer.parseInt(args[1]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[2]);  // For a bin
        final long GENOME_LENGTH_HAPLOID = Long.parseLong(args[3]);
        final String OUTPUT_PREFIX = args[4];
        
        final int N_BINS = MAX_READ_LENGTH/BIN_LENGTH+1;
        final int HALF_BIN_LENGTH = BIN_LENGTH>>1;
        final int MIN_READS_IN_LOCAL_MAXIMUM = 500;  // Arbitrary
        
        int i, j;
        int bin, last, length, maxLength, nReads, minBin, minValue;
        int nMaxima, previousMaximum, lastMaximum;
        double sum, mass;
        String str, file, header, sequence, separator, quality;
        Random random = new Random();
        BufferedReader br1, br2;
        BufferedWriter bw;
        int[] buffers_last;
        long[] buffers_stringLength;
        String[] tmpArray, tokens;
        BufferedWriter[] buffers;
        
        // Building bins
        buffers = new BufferedWriter[N_BINS];
        for (i=0; i<N_BINS; i++) buffers[i] = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".bin"));
        buffers_last = new int[N_BINS];
        Arrays.fill(buffers_last,-1);
        buffers_stringLength = new long[N_BINS];
        Arrays.fill(buffers_stringLength,0);
        br1 = new BufferedReader(new FileReader(INPUT_FILES_LIST));
        file=br1.readLine();
        while (file!=null) {
            System.err.println("Loading file <"+file+">...");
            br2 = new BufferedReader(new FileReader(file));
            header=br2.readLine(); nReads=0;
            while (header!=null) {
                sequence=br2.readLine(); separator=br2.readLine(); quality=br2.readLine();
                bin=sequence.length()/BIN_LENGTH;
                if (bin>=N_BINS) bin=N_BINS-1;
                buffers_last[bin]++;
                buffers[bin].write(header); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers[bin].write(sequence); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers_stringLength[bin]+=sequence.length();
                buffers[bin].write(separator); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers[bin].write(quality); buffers[bin].newLine();
                nReads++;
                if (nReads%10000==0) System.err.println(nReads+" reads");
                header=br2.readLine();
            }
            br2.close();
            System.err.println("Done");
            file=br1.readLine();
        }
        br1.close();
        for (i=0; i<N_BINS; i++) buffers[i].close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+".histogram"));
        maxLength=0;
        for (i=0; i<N_BINS; i++) {
            length=buffers_last[i]+1;
            if (length>maxLength) maxLength=length;
            bw.write((i*BIN_LENGTH+HALF_BIN_LENGTH)+","+length+"\n");
        }
        bw.close();
        
        // Permuting reads in each bin
        System.err.println("Permuting bins...");
        tmpArray = new String[maxLength];
        for (i=0; i<N_BINS; i++) {
            br1 = new BufferedReader(new FileReader(OUTPUT_PREFIX+i+".bin"));
            str=br1.readLine(); last=-1;
            while (str!=null) {
                tmpArray[++last]=str;
                str=br1.readLine();
            }
            br1.close();
            if (last<=1) continue;
            shuffle(tmpArray,last,random);
            bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".bin"));
            for (j=0; j<=last; j++) {
                bw.write(tmpArray[j]); bw.newLine();
            }
            bw.close();
        }
        System.err.println("Done");
        
        // Reporting high local maxima
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+".maxima"));
        nMaxima=0; lastMaximum=-1; previousMaximum=-1;
        for (i=0; i<N_BINS; i++) {
            if (buffers_last[i]+1>=MIN_READS_IN_LOCAL_MAXIMUM && (i==0 || buffers_last[i]>buffers_last[i-1]) && (i==N_BINS-1 || buffers_last[i]>buffers_last[i+1])) {
                bw.write((i*BIN_LENGTH)+","+((i+1)*BIN_LENGTH-1)+"\n");
                nMaxima++;
                previousMaximum=lastMaximum; lastMaximum=i;
            }
        }
        bw.close();
        System.err.println("The read length distribution contains "+nMaxima+" local maxima with >="+MIN_READS_IN_LOCAL_MAXIMUM+" reads each.");
        if (nMaxima>=2) {
            System.err.println("The two rightmost high maxima are at bins: "+previousMaximum+" (lengths ["+(previousMaximum*BIN_LENGTH)+".."+((previousMaximum+1)*BIN_LENGTH-1)+"]) and "+lastMaximum+" (lengths ["+(lastMaximum*BIN_LENGTH)+".."+((lastMaximum+1)*BIN_LENGTH-1)+"])");
            minBin=-1; minValue=Integer.MAX_VALUE;
            for (i=previousMaximum+1; i<lastMaximum; i++) {
                if (buffers_last[i]<buffers_last[i-1] && buffers_last[i]<buffers_last[i+1] && buffers_last[i]<minValue) {
                    minBin=i; minValue=buffers_last[i];
                }
            }
            System.err.println("Assume to cut the distribution at the "+minBin+"-th bin of smallest local minimum (lengths ["+(minBin*BIN_LENGTH)+".."+((minBin+1)*BIN_LENGTH-1)+"]).");
            sum=0.0;
            for (i=0; i<N_BINS; i++) sum+=buffers_last[i]+1;
            mass=0.0;
            for (i=0; i<=minBin; i++) mass+=buffers_last[i]+1;
            mass/=sum;
            System.err.println("Probability mass up to the "+minBin+"-th bin, included: "+mass);
            System.err.println("Probability mass after the "+minBin+"-th bin: "+(1.0-mass));
            sum=0.0;
            for (i=0; i<=minBin; i++) sum+=buffers_stringLength[i];
            System.err.println("Coverage of one haplotype by bins up to the "+minBin+"-th, included: "+(sum/GENOME_LENGTH_HAPLOID));
            sum=0.0;
            for (i=minBin+1; i<N_BINS; i++) sum+=buffers_stringLength[i];
            System.err.println("Coverage of one haplotype by bins after the "+minBin+"-th: "+(sum/GENOME_LENGTH_HAPLOID));
        }
        else if (nMaxima==1) System.err.println("The only high local maximum is at bin: "+lastMaximum+" (lengths ["+(lastMaximum*BIN_LENGTH)+".."+((lastMaximum+1)*BIN_LENGTH-1)+"])");
        else { /* NOP */ }
	}
    
    
    private static final void shuffle(String[] array, int last, Random random) {
        int i, j;
        String tmp;
        
        for (i=last; i>0; i--) {
            j=random.nextInt(i+1);
            tmp=array[j]; array[j]=array[i]; array[i]=tmp;
        }
    }

}