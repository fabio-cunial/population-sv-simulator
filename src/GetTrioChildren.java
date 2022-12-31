import java.util.Arrays;
import java.io.*;


/**
 * Given a trio table describing families, the program keeps all individuals 
 * from consistent families (both parents and children), and for every such 
 * individual it creates a text file that lists the addresses of all the FASTQ 
 * files of the individual taken from a Terra table (so this is just the join of
 * a trio table and a Terra table). The program prints also: (1) for each trio 
 * child, a file that lists the IDs of its parents (one per line); (2) a list of
 * all IDs of trio children (one per line).
 */
public class GetTrioChildren {
    /**
     * Trio table constants
     */
    private static final String TRIO_TABLE_DATE_SEPARATOR = "-";
    private static final String TRIO_TABLE_GENDER_1 = "Female";
    private static final String TRIO_TABLE_GENDER_2 = "Male";
    private static final String TRIO_TABLE_FAMILY_ID_HEADER = "family_id";
	private static final String TRIO_TABLE_INDIVIDUAL_ID_HEADER = "research_id";
    private static final String TRIO_TABLE_DATE_HEADER = "date_of_birth";
    private static final String TRIO_TABLE_SEX_AT_BIRTH_HEADER = "sex_at_birth";
    private static final int TRIO_TABLE_GENDER_FEMALE = 0;
    private static final int TRIO_TABLE_GENDER_MALE = 1;
    private static final int TRIO_TABLE_GENDER_OTHER = 2;
    
    /**
     * Terra table constants
     */
    private static final String TERRA_TABLE_FASTQ_HEADER = "hifi_fq";
    
    /**
     * Constants of both tables
     */
    private static final String SEPARATOR = "\t";
    private static final String TERRA_TABLE_ID_HEADER = "research_id";
    
    /**
     * All children that belong to consistent families in the trio table, with 
     * their IDs in the Terra table and pointers to their parents.
     */
    private static Individual[] children;
    private static int children_last;
    
    /**
     * The FASTQ files (columns) of each individual in a valid trio (rows, 
     * sorted in increasing order).
     */
    private static String[] individualIDs;
    private static int individualIDs_last;
    private static String[][] individual2fastqs;
    private static int[] individual2fastqs_last;
    
    
	public static void main(String[] args) throws IOException {
		final String TRIO_TABLE = args[0];
        final String TERRA_TABLE = args[1];
        final String OUTPUT_PREFIX = args[2];
        
        System.err.println("Loading trio table...");
        loadChildren(TRIO_TABLE);
        System.err.println("Done");
        System.err.println("Loading Terra table...");
        loadFastqs(TERRA_TABLE);
        System.err.println("Done");
        printOutput(OUTPUT_PREFIX);
	}
    
    
    private static final void loadChildren(String trioTable) throws IOException {
        final int MIN_YEAR = 1800;
        int i, j;
        int sexAtBirth, year, familyID, currentFamilyID, family_last;
        int FAMILY_ID_FIELD, INDIVIDUAL_ID_FIELD, DATE_FIELD, SEX_AT_BIRTH_FIELD, TERRA_TABLE_ID_FIELD;
        String str;
        BufferedReader br;
        String[] tokens;
        Individual[] family;
        
        FAMILY_ID_FIELD=-1; INDIVIDUAL_ID_FIELD=-1; DATE_FIELD=-1; SEX_AT_BIRTH_FIELD=-1; TERRA_TABLE_ID_FIELD=-1;
        br = new BufferedReader(new FileReader(trioTable));
        str=br.readLine();
        tokens=str.split(SEPARATOR);
        for (i=0; i<tokens.length; i++) {
            if (tokens[i].equalsIgnoreCase(TRIO_TABLE_FAMILY_ID_HEADER)) FAMILY_ID_FIELD=i;
            if (tokens[i].equalsIgnoreCase(TRIO_TABLE_INDIVIDUAL_ID_HEADER)) INDIVIDUAL_ID_FIELD=i;
            if (tokens[i].equalsIgnoreCase(TRIO_TABLE_DATE_HEADER)) DATE_FIELD=i;
            if (tokens[i].equalsIgnoreCase(TRIO_TABLE_SEX_AT_BIRTH_HEADER)) SEX_AT_BIRTH_FIELD=i;
            if (tokens[i].equalsIgnoreCase(TERRA_TABLE_ID_HEADER)) TERRA_TABLE_ID_FIELD=i;            
        }
        if (FAMILY_ID_FIELD==-1 || INDIVIDUAL_ID_FIELD==-1 || DATE_FIELD==-1 || SEX_AT_BIRTH_FIELD==-1 || TERRA_TABLE_ID_FIELD==-1) {
            System.err.println("ERROR: invalid header in the trio table.");
            System.err.println("FAMILY_ID_FIELD="+FAMILY_ID_FIELD+" INDIVIDUAL_ID_FIELD="+INDIVIDUAL_ID_FIELD+" DATE_FIELD="+DATE_FIELD+" SEX_AT_BIRTH_FIELD="+SEX_AT_BIRTH_FIELD+" TERRA_TABLE_ID_FIELD="+TERRA_TABLE_ID_FIELD);
            System.err.println("Header: "+str);
            System.exit(1);
        }
        children = new Individual[10];  // Arbitrary, resizable.
        children_last=-1;
        family = new Individual[100];  // Arbitrary, never resized.
        for (i=0; i<family.length; i++) family[i] = new Individual();
        str=br.readLine(); currentFamilyID=-1; family_last=-1;
        while (str!=null) {
            tokens=str.split(SEPARATOR);
            familyID=Integer.parseInt(tokens[FAMILY_ID_FIELD]);
            if (tokens[SEX_AT_BIRTH_FIELD].equalsIgnoreCase(TRIO_TABLE_GENDER_1)) sexAtBirth=TRIO_TABLE_GENDER_FEMALE;
            else if (tokens[SEX_AT_BIRTH_FIELD].equalsIgnoreCase(TRIO_TABLE_GENDER_2)) sexAtBirth=TRIO_TABLE_GENDER_MALE;
            else sexAtBirth=TRIO_TABLE_GENDER_OTHER;
            if (tokens[DATE_FIELD].indexOf(TRIO_TABLE_DATE_SEPARATOR)>=0) year=Integer.parseInt(tokens[DATE_FIELD].substring(0,tokens[DATE_FIELD].indexOf(TRIO_TABLE_DATE_SEPARATOR)));
            else year=MIN_YEAR;
            if (familyID==currentFamilyID) family_last++;
            else {
                processFamily(family,family_last,currentFamilyID);
                currentFamilyID=familyID;
                family_last=0;
            }
            family[family_last].id=Integer.parseInt(tokens[INDIVIDUAL_ID_FIELD]);
            family[family_last].gender=sexAtBirth;
            family[family_last].year=year;
            family[family_last].terraTable_field=tokens[TERRA_TABLE_ID_FIELD];
            str=br.readLine();
        }
        br.close();
        processFamily(family,family_last,currentFamilyID);
        if (children_last>0) Arrays.sort(children,0,children_last+1);
        j=0;
        for (i=1; i<=children_last; i++) {
            if (children[i].terraTable_field==children[j].terraTable_field) continue;
            children[++j]=children[i];
        }
        children_last=j;
    }
    
    
    /**
     * Adds to $children$ all and only the children of a consistent family.
     * Warns and continues if $family[0..last]$ is not consistent.
     */
    private static final void processFamily(Individual[] family, int last, int familyID) {
        final int MIN_CONSISTENT_GAP = 12;  // Arbitrary
        int i, j;
        int output_last, gap, maxGap;
                
        if (last<2) return;
        Arrays.sort(family,0,last+1);
        maxGap=0;
        for (i=2; i<=last; i++) {
            gap=family[i].year-family[1].year;
            if (gap>maxGap) maxGap=gap;
        }
        if (family[0].id==family[1].id || family[0].gender==family[1].gender || family[0].gender==TRIO_TABLE_GENDER_OTHER || family[1].gender==TRIO_TABLE_GENDER_OTHER || maxGap<MIN_CONSISTENT_GAP) {
            System.err.println("WARNING: family "+familyID+" is not consistent.");
            for (i=0; i<=last; i++) System.err.println(family[i]);
            return;
        }
        for (i=2; i<=last; i++) {
            children_last++;
            if (children_last==children.length) {
                Individual[] newArray = new Individual[children.length<<1];
                System.arraycopy(children,0,newArray,0,children.length);
                children=newArray;
            }
            children[children_last]=family[i].clone();
            children[children_last].parent1=family[0].clone();
            children[children_last].parent2=family[1].clone();
        }
    }

    
    private static class Individual implements Comparable {
        public static int order;
        public int id, gender, year;
        public String terraTable_field;
        public Individual parent1, parent2;
        public String[] fastqs;
        
        public String toString() {
            return "id="+id+" gender="+gender+" year="+year;
        }
        
        public Individual clone() {
            Individual out = new Individual();
            out.id=id; out.gender=gender; out.year=year;
            out.terraTable_field=terraTable_field;
            out.parent1=parent1; out.parent2=parent2;
            out.fastqs=fastqs;
            return out;
        }
        
        /**
         * By increasing year
         */
        public int compareTo(Object other) {
            Individual otherIndividual = (Individual)other;
            if (year<otherIndividual.year) return -1;
            else if (year>otherIndividual.year) return 1;
            else return 0;
        }
    }
    
    
    /**
     * Loads $individualIDs$ and $individual2fastqs$.
     */
    private static final void loadFastqs(String terraTable) throws IOException {
        final int CAPACITY = 2;  // Arbitrary
        int i, j;
        int INDIVIDUAL_ID_FIELD, FASTQ_FIELD;
        String str, childID;
        BufferedReader br;
        String[] tokens;
        
        // Sorting and compacting the distinct IDs of every child and parent
        individualIDs = new String[children_last*3];
        individualIDs_last=-1;
        for (i=0; i<=children_last; i++) {
            individualIDs[++individualIDs_last]=children[i].terraTable_field;
            individualIDs[++individualIDs_last]=children[i].parent1.terraTable_field;
            individualIDs[++individualIDs_last]=children[i].parent2.terraTable_field;
        }
        Arrays.sort(individualIDs,0,individualIDs_last+1);
        j=0;
        for (i=1; i<=individualIDs_last; i++) {
            if (!individualIDs[i].equalsIgnoreCase(individualIDs[j])) individualIDs[++j]=individualIDs[i];
        }
        individualIDs_last=j;
        
        // Building FASTQ lists
        br = new BufferedReader(new FileReader(terraTable));
        str=br.readLine();
        tokens=str.split(SEPARATOR);
        INDIVIDUAL_ID_FIELD=-1; FASTQ_FIELD=-1;
        for (i=0; i<tokens.length; i++) {
            if (tokens[i].equalsIgnoreCase(TERRA_TABLE_ID_HEADER)) INDIVIDUAL_ID_FIELD=i;
            if (tokens[i].equalsIgnoreCase(TERRA_TABLE_FASTQ_HEADER)) FASTQ_FIELD=i;
        }
        if (INDIVIDUAL_ID_FIELD==-1 || FASTQ_FIELD==-1) {
            System.err.println("ERROR: invalid header in Terra table.");
            System.exit(1);
        }
        individual2fastqs = new String[individualIDs_last+1][CAPACITY];
        individual2fastqs_last = new int[individualIDs_last+1];
        Arrays.fill(individual2fastqs_last,-1);
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(SEPARATOR);
            childID=tokens[INDIVIDUAL_ID_FIELD];
            i=Arrays.binarySearch(individualIDs,0,individualIDs_last+1,childID);
            if (i<0) {
                str=br.readLine();
                continue;
            }
            individual2fastqs_last[i]++;
            if (individual2fastqs_last[i]==individual2fastqs[i].length) {
                String[] newArray = new String[individual2fastqs[i].length<<1];
                System.arraycopy(individual2fastqs[i],0,newArray,0,individual2fastqs[i].length);
                individual2fastqs[i]=newArray;
            }
            individual2fastqs[i][individual2fastqs_last[i]]=tokens[FASTQ_FIELD];
            str=br.readLine();
        }
        br.close();
        for (i=0; i<=individualIDs_last; i++) {
            if (individual2fastqs_last[i]==-1) System.err.println("WARNING: individual "+individualIDs[i]+" has no FASTQ.");
        }
    }
    
    
    /**
     * Remark: only individuals with some FASTQ file are reported in output, and
     * only children such that both them and their parents have some FASTQ file.
     */
    private static final void printOutput(String prefix) throws IOException {
        int i, j, k, h;
        BufferedWriter bw1, bw2;
        
        // List of trio children, and parents of each child
        bw1 = new BufferedWriter(new FileWriter(prefix+"children.txt"));
        for (i=0; i<=children_last; i++) {
            j=Arrays.binarySearch(individualIDs,0,individualIDs_last+1,children[i].id);
            k=Arrays.binarySearch(individualIDs,0,individualIDs_last+1,children[i].parent1.id);
            h=Arrays.binarySearch(individualIDs,0,individualIDs_last+1,children[i].parent2.id);
            if (individual2fastqs_last[j]!=-1 && individual2fastqs_last[k]!=-1 && individual2fastqs_last[h]!=-1) {
                bw1.write(children[i].id+"\n");
                bw2 = new BufferedWriter(new FileWriter(prefix+children[i].id+".parents"));
                bw2.write(children[i].parent1.id+"\n");
                bw2.write(children[i].parent2.id+"\n");
                bw2.close();
            }
        }
        bw1.close();
        
        // FASTQs of each individual
        for (i=0; i<=individualIDs_last; i++) {
            if (individual2fastqs_last[i]!=-1) {
                bw1 = new BufferedWriter(new FileWriter(prefix+individualIDs[i]+".fastqs"));
                for (j=0; j<=individual2fastqs_last[i]; j++) bw1.write(individual2fastqs[i][j]+"\n");
                bw1.close();
            }
        }
    }
    
}
