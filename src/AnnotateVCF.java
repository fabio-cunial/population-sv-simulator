import java.util.Arrays;
import java.io.*;
import java.text.NumberFormat;


/**
 * Adds fields REPEATS_START, REPEATS_END, REPEATS_FRACTION (defined in 
 * $SimulateHaplotypes.java$) to every VCF record.
 */
public class AnnotateVCF {
    /**
     * VCF constants
     */
    private static final char COMMENT = '#';
    private static final String SEPARATOR = ";";
    private static final String SVLEN_STR = "SVLEN";
    private static final String END_STR = "END";
    
	/**
	 * Confidence intervals of positions.
	 *
	 * Remark: some callers report a standard deviation instead of a confidence
	 * interval. Sniffles1 reports additional interval information in its BEDPE
	 * output (an alternative to VCF), which we disregard for simplicity. Some 
     * callers use CILEN to express a "confidence interval around inserted/
     * deleted material between breakends": we interpret CILEN exactly like 
     * CIEND.
	 */
    private static final String CI_SEPARATOR = ",";
    private static final String CIPOS_STR = "CIPOS";
    private static final int CIPOS_STR_LENGTH = CIPOS_STR.length();
    private static final String STD_START1_STR = "STD_quant_start";
    private static final int STD_START1_STR_LENGTH = STD_START1_STR.length();
    private static final String STD_START2_STR = "STD_POS1";
    private static final int STD_START2_STR_LENGTH = STD_START2_STR.length();
    private static final String CIEND_STR = "CIEND";
    private static final int CIEND_STR_LENGTH = CIEND_STR.length();
    private static final String CILEN_STR = "CILEN";
    private static final int CILEN_STR_LENGTH = CILEN_STR.length();
    private static final String STD_END1_STR = "STD_quant_stop";
    private static final int STD_END1_STR_LENGTH = STD_END1_STR.length();
    private static final String STD_END2_STR = "STD_POS2";
    private static final int STD_END2_STR_LENGTH = STD_END2_STR.length();
    
    
    /**
     * @param args to compute REPEATS_{START,END}, the procedure takes the 
     * union of the repeat types of all the positions in the uncertainty range 
     * declared in the record, excluding positions that are farther than 
     * $args[3]$ from the start/end. If the uncertainty range is too short in a
     * direction, the program forces it to be at least $args[2]$.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILE = args[0];
        final String OUTPUT_FILE = args[1];
        final int MIN_DISTANCE = Integer.parseInt(args[2]);
        final int MAX_DISTANCE = Integer.parseInt(args[3]);
		final String REFERENCE_FILE = args[4];
        final String REPEAT_INTERVALS_FILE = args[5];
        final String SEGDUPS_FILE = args[6];
        
        final int MAX_LENGTH_DIFFERENCE = 100;  // Arbitrary
        int i;
		int position, end, length, intervalFirst, intervalLast, repeatsStart, repeatsEnd;
        double repeatsFraction;
		String str, strPrime;
		BufferedReader br;
        BufferedWriter bw;
        NumberFormat formatter;
        boolean[] tmpBoolean = new boolean[1+((MAX_DISTANCE)<<1)];
        int[] tmpInt = new int[1+((MAX_DISTANCE)<<1)];
		int[] out = new int[2];
		String[] tokens;
        
        formatter = NumberFormat.getInstance();
        formatter.setMinimumFractionDigits(4); formatter.setMaximumFractionDigits(4);
        SimulateHaplotypes.loadReference(REFERENCE_FILE);
		BuildModel.deserializeRepeatIntervals(REPEAT_INTERVALS_FILE);
		BuildModel.deserializeSegdups(SEGDUPS_FILE);
        SimulateHaplotypes.loadPositions(SimulateHaplotypes.referenceLength,BuildModel.MAX_DISTANCE_SR);
		br = new BufferedReader(new FileReader(INPUT_FILE));
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
		str=br.readLine();
		while (str!=null) {
			if (str.charAt(0)==COMMENT) {
                bw.write(str); bw.newLine();
				str=br.readLine();
				continue;
			}
			tokens=str.split("\t");
			position=Integer.parseInt(tokens[1]);
            strPrime=getField(tokens[7],END_STR);
            if (strPrime!=null) {
                end=Integer.parseInt(strPrime)-1;
                length=end-position+1;
                strPrime=getField(tokens[7],SVLEN_STR);
                if (strPrime!=null && Math.abs(Integer.parseInt(strPrime)-length)>MAX_LENGTH_DIFFERENCE) {
                    System.err.println("WARNING: SVLEN != END-POS. Using END. "+str);
                }
            }
            else {
                strPrime=getField(tokens[7],SVLEN_STR);
                if (strPrime!=null) {
                    length=Integer.parseInt(strPrime);
                    end=position+length-1;
                }
                else {
                    System.err.println("ERROR: neither END nor SVLEN are present. "+str);
                    end=Integer.MAX_VALUE; length=Integer.MAX_VALUE;
                }
            }
			getConfidenceInterval(tokens[7],false,out);
			intervalFirst=position+out[0]; intervalLast=position+out[1];
            if (intervalFirst<position-MAX_DISTANCE) intervalFirst=position-MAX_DISTANCE;
            if (intervalFirst>position-MIN_DISTANCE) intervalFirst=position-MIN_DISTANCE;
            if (intervalLast>position+MAX_DISTANCE) intervalLast=position+MAX_DISTANCE;
            if (intervalLast<position+MIN_DISTANCE) intervalLast=position+MIN_DISTANCE;
            repeatsStart=getRepeatType(intervalFirst,intervalLast,tmpBoolean);
            if (end!=Integer.MAX_VALUE) {
                repeatsFraction=getRepeatSurface(position,end,tmpInt)/((double)length);
                getConfidenceInterval(tokens[7],true,out);
                intervalFirst=end+out[0]; intervalLast=end+out[1];
                if (intervalFirst<position-MAX_DISTANCE) intervalFirst=position-MAX_DISTANCE;
                if (intervalFirst>position-MIN_DISTANCE) intervalFirst=position-MIN_DISTANCE;
                if (intervalLast>position+MAX_DISTANCE) intervalLast=position+MAX_DISTANCE;
                if (intervalLast<position+MIN_DISTANCE) intervalLast=position+MIN_DISTANCE;
                repeatsEnd=getRepeatType(intervalFirst,intervalLast,tmpBoolean);
            }
            else {
                repeatsFraction=0.0; repeatsEnd=4;  // Arbitrary
            }
            tokens[7]+=SEPARATOR+"REPEATS_START="+repeatsStart+SEPARATOR+"REPEATS_END="+repeatsEnd+SEPARATOR+"REPEATS_FRACTION="+formatter.format(repeatsFraction);
            bw.write(tokens[0]);
            for (i=1; i<tokens.length; i++) bw.write("\t"+tokens[i]);
            bw.newLine();
			str=br.readLine();
		}
		br.close(); bw.close();
    }
    
    
	/**
	 * @return NULL if $field$ does not occur in $str$.
	 */
	private static final String getField(String str, String field) {
		final int FIELD_LENGTH = field.length()+1;
		int p = str.indexOf(field+"=");
		if (p<0) return null;
		if (field.equalsIgnoreCase(END_STR)) {
			while (p>=2 && str.substring(p-2,p-2+CIEND_STR.length()).equalsIgnoreCase(CIEND_STR)) p=str.indexOf(field+"=",p+1);
			if (p<0) return null;
		}
		final int q = str.indexOf(SEPARATOR,p+FIELD_LENGTH);
		return str.substring(p+FIELD_LENGTH,q<0?str.length():q);
	}
    
    
	/**
	 * @param posID FALSE=first position of a VCF record; TRUE=last position;
	 * @param out 0=quantity to be added to the position to get the first 
	 * position of the interval (typically negative or zero); 1=quantity to be
	 * added to the position to get the last position of the interval (typically
	 * positive or zero).
	 */
	private static final void getConfidenceInterval(String str, boolean posID, int[] out) {
		final int SIGMA_MULTIPLE = 3;  // 2 or 3 captures most of a Gaussian
		char c;
		int p, q;
		int length, quantum;
		
		if (!posID) {
			length=CIPOS_STR_LENGTH;
			p=str.indexOf(CIPOS_STR);
			if (p>=0) {
				p+=length+1; q=str.indexOf(CI_SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c) || c=='-') {
					// Sign already ok
					if (q>=0) out[0]=Integer.parseInt(str.substring(p,q));
					else out[0]=Integer.parseInt(str.substring(p));
				}
				else out[0]=0;
				p=q+1; q=str.indexOf(SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c) || c=='-') {
					// Sign already ok
					if (q>=0) out[1]=Integer.parseInt(str.substring(p,q));
					else out[1]=Integer.parseInt(str.substring(p));
				}
				else out[1]=0;
				return;
			}
			length=STD_START1_STR_LENGTH;
			p=str.indexOf(STD_START1_STR);
			if (p>=0) {
				p+=length+1; q=str.indexOf(SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c)) {
					if (q>=0) quantum=(int)Math.round(Double.parseDouble(str.substring(p,q))*SIGMA_MULTIPLE);
					else quantum=(int)Math.round(Double.parseDouble(str.substring(p))*SIGMA_MULTIPLE);
					out[0]=-quantum; out[1]=quantum;
				}
				else { out[0]=0; out[1]=0; }
				return;
			}
			length=STD_START2_STR_LENGTH;
			p=str.indexOf(STD_START2_STR);
			if (p>=0) {
				p+=length+1; q=str.indexOf(SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c)) {
					if (q>=0) quantum=(int)Math.round(Double.parseDouble(str.substring(p,q))*SIGMA_MULTIPLE);
					else quantum=(int)Math.round(Double.parseDouble(str.substring(p))*SIGMA_MULTIPLE);
					out[0]=-quantum; out[1]=quantum;
				}
				else { out[0]=0; out[1]=0; }
				return;
			}
		}
		else {
			length=CIEND_STR_LENGTH;
			p=str.indexOf(CIEND_STR);
			if (p>=0) {
				p+=length+1; q=str.indexOf(CI_SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c) || c=='-') {
					// Sign already ok
					if (q>=0) out[0]=Integer.parseInt(str.substring(p,q));
					else out[0]=Integer.parseInt(str.substring(p));
				}
				else out[0]=0;
				p=q+1; q=str.indexOf(SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c) || c=='-') {
					// Sign already ok
					if (q>=0) out[1]=Integer.parseInt(str.substring(p,q));
					else out[1]=Integer.parseInt(str.substring(p));
				}
				else out[1]=0;
				return;
			}
			length=CILEN_STR_LENGTH;
			p=str.indexOf(CILEN_STR);
			if (p>=0) {
				p+=length+1; q=str.indexOf(CI_SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c)) {
					// Sign already ok
					if (q>=0) out[0]=Integer.parseInt(str.substring(p,q));
					else out[0]=Integer.parseInt(str.substring(p));
				}
				else out[0]=0;
				p=q+1; q=str.indexOf(SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c)) {
					// Sign already ok
					if (q>=0) out[1]=Integer.parseInt(str.substring(p,q));
					else out[1]=Integer.parseInt(str.substring(p));
				}
				else out[1]=0;
				return;
			}
			length=STD_END1_STR_LENGTH;
			p=str.indexOf(STD_END1_STR);
			if (p>=0) {
				p+=length+1; q=str.indexOf(SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c)) {
					if (q>=0) quantum=(int)Math.round(Double.parseDouble(str.substring(p,q))*SIGMA_MULTIPLE);
					else quantum=(int)Math.round(Double.parseDouble(str.substring(p))*SIGMA_MULTIPLE);
					out[0]=-quantum; out[1]=quantum;
				}
				else { out[0]=0; out[1]=0; }
				return;
			}
			length=STD_END2_STR_LENGTH;
			p=str.indexOf(STD_END2_STR);
			if (p>=0) {
				p+=length+1; q=str.indexOf(SEPARATOR,p);
				c=str.charAt(p);
				if (Character.isDigit(c)) {
					if (q>=0) quantum=(int)Math.round(Double.parseDouble(str.substring(p,q))*SIGMA_MULTIPLE);
					else quantum=(int)Math.round(Double.parseDouble(str.substring(p))*SIGMA_MULTIPLE);
					out[0]=-quantum; out[1]=quantum;
				}
				else { out[0]=0; out[1]=0; }
				return;
			}
		}
		out[0]=0; out[1]=0;
	}
    
    
    /**
     * @return the ID of the union of all repeat types of all positions in
     * $[first..last]$.
     */
    private static final int getRepeatType(int first, int last, boolean[] tmpArray) {
        int i, j;
        
        for (i=0; i<=4; i++) {
            j=Arrays.binarySearch(SimulateHaplotypes.positions[i],0,SimulateHaplotypes.positions_last[i]+1,first);
            if (j<0) j=-1-j;
            tmpArray[i]=nonemptyIntersection(SimulateHaplotypes.positions[i],j,SimulateHaplotypes.positions_last[i],first,last);
        }
        if (tmpArray[2] || (tmpArray[0] && tmpArray[1])) return 2;
        else if (tmpArray[0]) return 0;
        else if (tmpArray[1]) return 1;
        else if (tmpArray[3]) return 3;
        else return 4;
    }
    
    
	/**
	 * @return TRUE iff $x1[from1..last1]$ contains an element in
     * $[from2..last2]$.
	 */
	private static final boolean nonemptyIntersection(int[] x1, int from1, int last1, int from2, int last2) {
		int i1, i2;
	
		i1=from1; i2=from2;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]<i2) i1++;
			else if (x1[i1]>i2) i2++;
			else return true;
		}
		return false;
	}

    
    /**
     * @param tmpArray temporary space of size at least $last-first+1$;
     * @return the number of basepairs inside $[first..last]$ that belong to
     * repeats of any type.
     */
    private static final int getRepeatSurface(int first, int last, int[] tmpArray) {
        int i, j, k;
        
        j=0;
        for (i=0; i<=3; i++) {
            k=Arrays.binarySearch(SimulateHaplotypes.positions[i],0,SimulateHaplotypes.positions_last[i]+1,first);
            if (k<0) k=-1-k;
            j=intersection(SimulateHaplotypes.positions[i],k,SimulateHaplotypes.positions_last[i],first,last,tmpArray,j);
        }
        if (j==0) return 1;
        Arrays.sort(tmpArray,0,j+1);
        k=0;
        for (i=1; i<=j; i++) {
            if (tmpArray[i]!=tmpArray[k]) tmpArray[++k]=tmpArray[i];
        }
        return k+1;
    }
    
    
	/**
	 * Stores in $out[outFrom..x]$ every element of $[from2..last2]$ that occurs
     * in $x1[from1..last1]$. $x$ is returned in output.
	 */
	private static final int intersection(int[] x1, int from1, int last1, int from2, int last2, int[] out, int outFrom) {
		int i1, i2, j;
	
		i1=from1; i2=from2; j=outFrom-1;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]<i2) i1++;
			else if (x1[i1]>i2) i2++;
			else out[++j]=i2;
		}
        return j;
	}

}