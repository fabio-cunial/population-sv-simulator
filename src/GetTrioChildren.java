import java.util.Arrays;
import java.io.*;


/**
 * For each trio child in a trio table, creates a text file that lists the 
 * addresses of all its FASTQ files taken from a Terra table (so this is just
 * the join of a trio table and a Terra table).
 */
public class GetTrioChildren {
    /**
     * Trio table constants
     */
    private static final String SEPARATOR = "\t";
    private static final String DATE_SEPARATOR = "-";
    private static final String GENDER_1 = "Female";
    private static final String GENDER_2 = "Male";
    private static final String FAMILY_ID_HEADER = "family_id";
	private static final String INDIVIDUAL_ID_HEADER = "research_id";
    private static final String DATE_HEADER = "date_of_birth";
    private static final String SEX_AT_BIRTH_HEADER = "sex_at_birth";
    private static final int GENDER_MALE = 0;
    private static final int GENDER_FEMALE = 1;
    private static final int GENDER_REST = 2;
    
    /**
     * Terra table constants
     */
    private static final String TERRA_TABLE_FASTQ_HEADER = "hifi_fq";
    
    /**
     * Constants of both tables
     */
    private static final String TERRA_TABLE_ID_HEADER = "research_id";
    
    /**
     * All children that belong to consistent families in the trio table, with 
     * their IDs in the Terra table and pointers to their parents.
     */
    private static Individual[] children;
    private static int children_last;
    
    /**
     * The FASTQ files (columns) of each trio child in $children$ (rows).
     */
    private static String[][] child2fastqs;
    private static int[] child2fastqs_last;
    
    
	public static void main(String[] args) throws IOException {
		final String TRIO_TABLE = args[0];
        final String TERRA_TABLE = args[1];
        final String OUTPUT_PREFIX = args[2];
        
        int i, j;
        int last, nChildren;
        BufferedWriter bw;
        
        System.err.println("Loading trio table...");
        loadChildrenIDs(TRIO_TABLE);
        System.err.println("Done");
        System.err.println("Loading Terra table...");
        loadChild2fastqs(TERRA_TABLE);
        System.err.println("Done");
        for (i=0; i<=childrenIDs_last; i++) {
            last=child2fastqs_last[i];
            if (last==-1) continue;
            bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+childrenIDs[i]+".list"));
            for (j=0; j<=last; j++) {
                bw.write(child2fastqs[i][j]);
                bw.newLine();
            }
            bw.close();
        }
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
            if (tokens[i].equalsIgnoreCase(FAMILY_ID_HEADER)) FAMILY_ID_FIELD=i;
            if (tokens[i].equalsIgnoreCase(INDIVIDUAL_ID_HEADER)) INDIVIDUAL_ID_FIELD=i;
            if (tokens[i].equalsIgnoreCase(DATE_HEADER)) DATE_FIELD=i;
            if (tokens[i].equalsIgnoreCase(SEX_AT_BIRTH_HEADER)) SEX_AT_BIRTH_FIELD=i;
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
            if (tokens[SEX_AT_BIRTH_FIELD].equalsIgnoreCase(GENDER_1)) sexAtBirth=0;
            else if (tokens[SEX_AT_BIRTH_FIELD].equalsIgnoreCase(GENDER_2)) sexAtBirth=1;
            else sexAtBirth=2;
            if (tokens[DATE_FIELD].indexOf(DATE_SEPARATOR)>=0) year=Integer.parseInt(tokens[DATE_FIELD].substring(0,tokens[DATE_FIELD].indexOf(DATE_SEPARATOR)));
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
            if (children[i].terraTable_field==childrenIDs[j].terraTable_field) continue;
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
        if (family[0].id==family[1].id || family[0].gender==family[1].gender || family[0].gender==GENDER_REST || family[1].gender==GENDER_REST || maxGap<MIN_CONSISTENT_GAP) {
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
    
    
    private static final void loadFastqs(String terraTable) throws IOException {
        final int CAPACITY = 2;  // Arbitrary
        int i, j;
        int ids_last;
        int CHILD_ID_FIELD, FASTQ_FIELD;
        String str, childID;
        BufferedReader br;
        String[] tokens, ids;
        
        // Sorting and compacting the distinct IDs of every child and parent
        ids = new String[children_last*3];
        ids_last=-1;
        for (i=0; i<=children_last; i++) {
            ids[++ids_last]=children[i].terraTable_field;
            ids[++ids_last]=children[i].parent1.terraTable_field;
            ids[++ids_last]=children[i].parent2.terraTable_field;
        }
        Arrays.sort(ids,0,ids_last+1);
        j=0;
        for (i=1; i<=ids_last; i++) {
            if (!ids[i].equalsIgnoreCase(ids[j])) ids[++j]=ids[i];
        }
        ids_last=j;
        
        // Building list files
        br = new BufferedReader(new FileReader(terraTable));
        str=br.readLine();
        tokens=str.split(SEPARATOR);
        CHILD_ID_FIELD=-1; FASTQ_FIELD=-1;
        for (i=0; i<tokens.length; i++) {
            if (tokens[i].equalsIgnoreCase(TERRA_TABLE_ID_HEADER)) CHILD_ID_FIELD=i;
            if (tokens[i].equalsIgnoreCase(TERRA_TABLE_FASTQ_HEADER)) FASTQ_FIELD=i;
        }
        if (CHILD_ID_FIELD==-1 || FASTQ_FIELD==-1) {
            System.err.println("ERROR: invalid header in Terra table.");
            System.exit(1);
        }
        child2fastqs = new String[childrenIDs_last+1][CAPACITY];
        child2fastqs_last = new int[childrenIDs_last+1];
        Arrays.fill(child2fastqs_last,-1);
        str=br.readLine();
        while (str!=null) {
            tokens=str.split(SEPARATOR);
            childID=tokens[CHILD_ID_FIELD];
            i=Arrays.binarySearch(childrenIDs,0,childrenIDs_last+1,childID);
            if (i<0) {
                str=br.readLine();
                continue;
            }
            child2fastqs_last[i]++;
            if (child2fastqs_last[i]==child2fastqs[i].length) {
                String[] newArray = new String[child2fastqs[i].length<<1];
                System.arraycopy(child2fastqs[i],0,newArray,0,child2fastqs[i].length);
                child2fastqs[i]=newArray;
            }
            child2fastqs[i][child2fastqs_last[i]]=tokens[FASTQ_FIELD];
            str=br.readLine();
        }
        br.close();
        for (i=0; i<=childrenIDs_last; i++) {
            if (child2fastqs_last[i]==-1) System.err.println("WARNING: child "+childrenIDs[i]+" has no FASTQ.");
        }
    }
    
}
