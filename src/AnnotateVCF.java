import java.util.Arrays;
import java.io.*;
import java.text.NumberFormat;


/**
 * Adds fields REPEATS_START, REPEATS_END, REPEATS_FRACTION (defined in 
 * $SimulateHaplotypes.java$) to every VCF record.
 *
 * Remark: the program takes ~11 GB of RAM.
 */
public class AnnotateVCF {
    /**
     * VCF constants
     */
    private static final char COMMENT = '#';
    private static final String SEPARATOR = ";";
    private static final String SVLEN_STR = "SVLEN";
    private static final String END_STR = "END";
    private static final String SVTYPE_STR = "SVTYPE";
    
	/**
	 * SV types
	 */
	private static final int TYPE_INSERTION = 1;
	private static final int TYPE_DELETION = 2;
	private static final int TYPE_DEL_INV = 3;
	private static final int TYPE_INVERSION = 4;
	private static final int TYPE_INV_DUP = 5;
	private static final int TYPE_DUPLICATION = 6;
	private static final int TYPE_CNV = 7;
	private static final int TYPE_BREAKEND = 8;
	private static final int TYPE_TRANSLOCATION = 9;
    
	/**
	 * SV types: labels used by callers.
	 */
	private static final String DEL_STR = "DEL";
	private static final String DEL_ME_STR = "DEL:ME";
	private static final String DEL_INV_STR = "DEL/INV";
	private static final String INS_STR = "INS";
	private static final String INS_ME_STR = "INS:ME";
	private static final String INS_NOVEL_STR = "INS:NOVEL";
	private static final String DUP_STR = "DUP";
	private static final String DUP_TANDEM_STR = "DUP:TANDEM";
	private static final String DUP_INT_STR = "DUP:INT";
	private static final String INV_STR = "INV";
	private static final String INV_DUP_STR = "INVDUP";
	private static final String CNV_STR = "CNV";
	private static final String BND_STR = "BND";
	private static final String TRA_STR = "TRA";
    
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
     * @param args 
     * 0: one path per line; the program can process a few thousands of files in
     *    a few minutes;
     * 1: to compute REPEATS_{START,END}, the procedure takes the union of the 
     *    repeat types of all the positions in the uncertainty range reported 
     *    in the record, excluding positions that are farther than $args[2]$ 
     *    from the start/end;
     * 2: if the uncertainty range is too short in one direction, the program 
     *    forces it to be at least $args[1]$.
     */
    public static void main(String[] args) throws IOException {
        final String INPUT_FILES_LIST = args[0];
        final int MIN_DISTANCE = Integer.parseInt(args[1]);
        final int MAX_DISTANCE = Integer.parseInt(args[2]);
		final String REFERENCE_FILE = args[3];
        final String REPEAT_INTERVALS_FILE = args[4];
        final String SEGDUPS_FILE = args[5];

		String str, strPrime;
		BufferedReader br;
        NumberFormat formatter;
        
        formatter = NumberFormat.getInstance();
        formatter.setMinimumFractionDigits(4); formatter.setMaximumFractionDigits(4);
        System.err.println("Loading data structures...");
        SimulateHaplotypes.loadReference(REFERENCE_FILE);
		BuildModel.deserializeRepeatIntervals(REPEAT_INTERVALS_FILE);
		BuildModel.deserializeSegdups(SEGDUPS_FILE);
        SimulateHaplotypes.loadPositions(SimulateHaplotypes.referenceLength,BuildModel.MAX_DISTANCE_SR);
        br = new BufferedReader(new FileReader(INPUT_FILES_LIST));
        str=br.readLine();
        while (str!=null) {
            strPrime=str.substring(0,str.indexOf(".vcf"))+"_annotated.vcf";
            annotateVCF(str,strPrime,MIN_DISTANCE,MAX_DISTANCE,formatter);
            str=br.readLine();
        }
        br.close();
    }
    
    
    private static final void annotateVCF(String inputFile, String outputFile, int minDistance, int maxDistance, NumberFormat formatter) throws IOException {
        final int MAX_LENGTH_DIFFERENCE = 100;  // Arbitrary
        int i, p;
		int position, end, length, lengthPrime, type;
        int intervalFirst, intervalLast, repeatsStart, repeatsEnd;
        double repeatsFraction;
		String str, strPrime;
		BufferedReader br;
        BufferedWriter bw;
        boolean[] tmpBoolean = new boolean[(1+maxDistance*2)*5];
        int[] tmpInt = new int[10000];  // Arbitrary
		int[] out = new int[2];
		String[] tokens;
        
        System.err.println("Annotating file "+inputFile);
		br = new BufferedReader(new FileReader(inputFile));
        bw = new BufferedWriter(new FileWriter(outputFile));
		str=br.readLine(); p=0;
		while (str!=null) {
            p++;
            if (p%1000==0) System.err.println("Processed "+p+" records");
			if (str.charAt(0)==COMMENT) {
                bw.write(str); bw.newLine();
				str=br.readLine();
				continue;
			}
			tokens=str.split("\t");
			type=getType_infoField(getField(tokens[7],SVTYPE_STR));
			if (type==-1) type=getType_altField(tokens[4]);
			if (type!=TYPE_DELETION && type!=TYPE_DUPLICATION && type!=TYPE_INSERTION && type!=TYPE_INVERSION) {
                System.err.println("WARNING: SV type not handled: "+str);
                bw.write(str); bw.newLine();
				str=br.readLine();
				continue;
			}
			position=Integer.parseInt(tokens[1]);
            strPrime=getField(tokens[7],END_STR);
            if (strPrime!=null) {
                end=Integer.parseInt(strPrime)-1;
                if (type!=TYPE_INSERTION) {
                    length=end-position+1;
                    strPrime=getField(tokens[7],SVLEN_STR);
                    if (strPrime!=null) {
                        lengthPrime=Integer.parseInt(strPrime);
                        if (lengthPrime<0) lengthPrime=-lengthPrime;
                        if (Math.abs(lengthPrime-length)>MAX_LENGTH_DIFFERENCE) {
                            System.err.println("WARNING: SVLEN != END-POS. Using END. "+str);
                        }
                    }
                }
                else length=Integer.MAX_VALUE;
            }
            else if (type!=TYPE_INSERTION) {
                strPrime=getField(tokens[7],SVLEN_STR);
                if (strPrime!=null) {
                    length=Integer.parseInt(strPrime);
                    if (length<0) length=-length;
                    end=position+length-1;
                }
                else {
                    System.err.println("ERROR: neither END nor SVLEN are present. "+str);
                    end=Integer.MAX_VALUE; length=Integer.MAX_VALUE;
                }
            }
            else { end=position; length=Integer.MAX_VALUE; }
			getConfidenceInterval(tokens[7],false,out);
			intervalFirst=position+out[0]; intervalLast=position+out[1];
            if (intervalFirst<position-maxDistance) intervalFirst=position-maxDistance;
            if (intervalFirst>position-minDistance) intervalFirst=position-minDistance;
            if (intervalLast>position+maxDistance) intervalLast=position+maxDistance;
            if (intervalLast<position+minDistance) intervalLast=position+minDistance;
            repeatsStart=getRepeatType(intervalFirst,intervalLast,tmpBoolean);
            if (end!=Integer.MAX_VALUE) {
                getConfidenceInterval(tokens[7],true,out);
                intervalFirst=end+out[0]; intervalLast=end+out[1];
                if (intervalFirst<position-maxDistance) intervalFirst=position-maxDistance;
                if (intervalFirst>position-minDistance) intervalFirst=position-minDistance;
                if (intervalLast>position+maxDistance) intervalLast=position+maxDistance;
                if (intervalLast<position+minDistance) intervalLast=position+minDistance;
                repeatsEnd=getRepeatType(intervalFirst,intervalLast,tmpBoolean);
                if (length!=Integer.MAX_VALUE) {
                    if ((end-position+1)*4>tmpInt.length) tmpInt = new int[(end-position+1)*4];
                    repeatsFraction=getRepeatSurface(position,end,tmpInt)/((double)length);
                }
                else repeatsFraction=0.0;
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
        System.err.println("Done annotating file "+inputFile);
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
	 * @return -1 iff the type cannot be determined.
	 */
	private static final int getType_infoField(String type) {
		if (type==null || type.length()==0) return -1;
		if ( type.equalsIgnoreCase(DEL_STR) || 
			 type.equalsIgnoreCase(DEL_ME_STR)
		   ) return TYPE_DELETION;
		else if (type.equalsIgnoreCase(DEL_INV_STR)) return TYPE_DEL_INV;
		else if ( type.equalsIgnoreCase(INS_STR) || 
			      type.equalsIgnoreCase(INS_ME_STR) || 
				  type.equalsIgnoreCase(INS_NOVEL_STR)
				) return TYPE_INSERTION;
		else if ( type.equalsIgnoreCase(DUP_STR) ||
			      type.equalsIgnoreCase(DUP_TANDEM_STR) ||
				  type.equalsIgnoreCase(DUP_INT_STR)
			    ) return TYPE_DUPLICATION;
		else if (type.equalsIgnoreCase(INV_STR)) return TYPE_INVERSION;
		else if (type.equalsIgnoreCase(INV_DUP_STR)) return TYPE_INV_DUP;
		else if (type.equalsIgnoreCase(CNV_STR)) return TYPE_CNV;
		else if (type.equalsIgnoreCase(BND_STR)) return TYPE_BREAKEND;
		else if (type.equalsIgnoreCase(TRA_STR)) return TYPE_TRANSLOCATION;
		else return -1;
	}
	
	
	private static final int getType_altField(String str) {
		if (str==null || str.length()==0) return -1;
		if ( str.equalsIgnoreCase("<"+DEL_STR+">") || 
			 str.equalsIgnoreCase("<"+DEL_ME_STR+">")
		   ) return TYPE_DELETION;
		else if (str.equalsIgnoreCase("<"+DEL_INV_STR+">")) return TYPE_DEL_INV;
		else if ( str.equalsIgnoreCase("<"+INS_STR+">") ||
			      str.equalsIgnoreCase("<"+INS_ME_STR+">") ||
				  str.equalsIgnoreCase("<"+INS_NOVEL_STR+">")
			    ) return TYPE_INSERTION;
		else if ( str.equalsIgnoreCase("<"+DUP_STR+">") ||
			      str.equalsIgnoreCase("<"+DUP_TANDEM_STR+">") ||
				  str.equalsIgnoreCase("<"+DUP_INT_STR+">")
			    ) return TYPE_DUPLICATION;
		else if (str.equalsIgnoreCase("<"+INV_STR+">")) return TYPE_INVERSION;
		else if (str.equalsIgnoreCase("<"+INV_DUP_STR+">")) return TYPE_INV_DUP;
		else if (str.equalsIgnoreCase("<"+CNV_STR+">")) return TYPE_CNV;
		else if (str.equalsIgnoreCase("<"+BND_STR+">")) return TYPE_BREAKEND;
		else if (str.equalsIgnoreCase("<"+TRA_STR+">")) return TYPE_TRANSLOCATION;
		else return -1;
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
        else if (tmpArray[4]) return 4;
        else {
            System.err.println("getRepeatType> WARNING: ["+first+".."+last+"] does not intersect with any position.");
            return 4;
        }
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
        
        j=-1;
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
	 * Stores in $out[lastOut+1..x]$ every element of $[from2..last2]$ that
     * occurs in $x1[from1..last1]$. $x$ is returned in output.
	 */
	private static final int intersection(int[] x1, int from1, int last1, int from2, int last2, int[] out, int lastOut) {
		int i1, i2, j;
    
		i1=from1; i2=from2; j=lastOut;
		while (i1<=last1 && i2<=last2) {
			if (x1[i1]<i2) i1++;
			else if (x1[i1]>i2) i2++;
			else {
                out[++j]=i2;
                i1++; i2++;
            }
		}
        return j;
	}

}