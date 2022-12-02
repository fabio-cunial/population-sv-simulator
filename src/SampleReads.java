import java.util.Arrays;
import java.util.Random;
import java.io.*;


/**
 * --------------->
 */
public class SampleReads {
    
    private static Random random;
    
    private static boolean[] loaded;
    private static String[][] buffers;
    private static int[] buffers_last;
    private static long[] bufferLength;
    private static long[] lastUsed;
    
    
    
	public static void main(String[] args) throws IOException {
		final String BUCKETS_PREFIX = args[0];
        final int N_BINS = Integer.parseInt(args[1]);
        final int BIN_LENGTH = Integer.parseInt(args[1]);
        final String OUTPUT_FILE = args[3];
        final long TARGET_COVERAGE = Integer.parseInt(args[4]);
        final long MAX_BPS_IN_RAM = Integer.parseInt(args[4]);
        
        final int CAPACITY = 1000;  // Arbitrary
        final String SEPARATOR = "|";
        
        int i;
        int bin;
        long coverage, nBpsInRam;
        String file, header, sequence, separator, quality;
        BufferedReader br1, br2;
        BufferedWriter bw;
        boolean[] loaded, marked;
        
        
        
        
        random = new Random();
        buffers = new String[N_BINS][CAPACITY];
        for (i=0; i<N_BINS; i++) buffers[i] = new StringBuilder();
        buffers_last = new int[N_BINS];
        Arrays.fill(buffers_last,-1);
        bufferLength = new long[N_BINS];
        Arrays.fill(bufferLength,0);
        loaded = new boolean[N_BINS];
        Arrays.fill(loaded,false);
        marked = new boolean[N_BINS];
        lastUsed = new long[N_BINS];
        Arrays.fill(lastUsed,0);
        coverage=0; nBpsInRam=0;
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
        t=-1;
        while (coverage<TARGET_COVERAGE) {
            t++;
            
            // Sampling a nonempty bin
            do { 
                bin=sample();
                if (!loaded[bin]) {
                    buffers_last[bin]=-1;
                    br = new BufferedReader(new FileReader(BUCKETS_PREFIX+bin+".fastq"));
                    header=br.readLine();
                    while (header!=null) {
                        buffers_last[bin]++;
                        if (buffers_last[bin]==buffers[bin].length) {
                            String[] newArray = new String[buffers[bin].length+CAPACITY];
                            System.arraycopy(buffers[bin],0,newArray,0,buffers[bin].length);
                            buffers[bin]=newArray;
                        }
                        sequence=br.readLine();
                        sep=br.readLine();
                        quality=br.readLine();
                        buffers[bin][buffers_last[bin]]=header+SEPARATOR+sequence+SEPARATOR+sep+SEPARATOR+quality;
                        bufferLength[bin]+=sequence.length();
                        nBpsInRam+=bufferLength[bin];
                        header=br.readLine();
                    }
                    br.close();
                    shuffle(bin);
                    Arrays.fill(marked,false);
                    marked[bin]=true;
                    while (nBpsInRam>MAX_BPS_IN_RAM) {
                        min=Long.MAX_VALUE; minIndex=-1;
                        for (i=0; i<N_BINS; i++) {
                            if (!marked[i] && lastUsed[i]<min) { min=lastUsed[i]; minIndex=i; }
                        }
                        nBpsInRam-=bufferLength[minIndex];
                        buffers[minIndex]=null; buffers_last[minIndex]=-1; bufferLength[minIndex]=0;
                        marked[minIndex]=true;
                    }
                }
                if (buffers_last[bin]!=-1) {
                    lastUsed[bin]=t;
                    break;
                }
            }
            
            // Sampling a read inside the bin
            i=random.nextInt(buffers_last[bin]+1);
            
            
            
            
            
            
            br2 = new BufferedReader(new FileReader(file));
            header=br2.readLine();
            while (str2!=null) {
                sequence=br2.readLine(); separator=br2.readLine(); quality=br2.readLine();
                bin=sequence.length()/BIN_LENGTH;
                if (bin>=N_BINS) bin=N_BINS-1;
                bufferLengths[bin]++;
                buffers[bin].write(header); buffers[bin].newLine();
                buffers[bin].write(sequence); buffers[bin].newLine("\n");
                buffers[bin].write(separator); buffers[bin].newLine("\n");
                buffers[bin].write(quality); buffers[bin].newLine("\n");
                header=br2.readLine();
            }
            br2.close();
            file=br1.readLine();
        }
        br1.close();
        for (i=0; i<N_BINS; i++) buffers[i].close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"histogram.txt"));
        for (i=0; i<N_BINS; i++) bw.write(bufferLengths[i]);
        bw.close();
	}
    
    
    
    
    
    

}