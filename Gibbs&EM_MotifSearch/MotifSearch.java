import java.util.*;
import java.io.*;
import java.text.*;  // DecimalFormat class

/******************************************************************
 * An instance of the <code>MotifSearch</code> class represents a
 * search for a motif common to a set of genomic sequences.
 ******************************************************************/
public class MotifSearch {

    /**************************************************************
     ********************** INSTANCE VARIABLES ********************
     **************************************************************/

    /**
     * A collection of genomic sequences
     */
    protected Vector<String> sequences;

    /**
     * The number of genomic sequences
     */
    protected int numSequences;

    /**
     * The length of the motif being searched for
     */
    protected int motifLength;

    /**
     * A 2D array corresponding to a position-specific scoring matrix (PSSM)
     * that is a model of the motif being searched for
     */
    protected double[][] matrix;

    /**
     * A collection of integers where each integer represents the index
     * of the start location of a motif instance in a genomic sequence
     */
    protected Vector<Integer> instanceLocations;

    /**
     * A dictionary of IUPAC symbols
     */
    private Hashtable<String, Character> IUPAC;



    /**************************************************************
     ********************** CONSTRUCTOR ***************************
     **************************************************************/

    /** 
     * Creates an initially empty <code>MotifSearch</code>.
     * <p>
     * A set of genomic sequences is determined from the specified <code>String</code> representing the name of a file.
     * Genomic sequences are read-in from the FASTA file. The <code>integer</code> parameter
     * represents the desired length of the motif being searched for.
     * Initially, the constructed
     * <code>MotifSearch</code> is empty.
     *
     * @param   fileName   the name of a FASTA file containing one or more genomic sequences
     * @param   motifLength  the length of the desired motif
     */
    public MotifSearch(String fileName, int motifLength) {
	sequences = readInMultipleSequencesFromFastaFile(fileName);
	numSequences = sequences.size();
	this.motifLength = motifLength;
	matrix = new double[4][motifLength];
	instanceLocations = new Vector<Integer>(numSequences);
	for (int i=0; i<numSequences; i++)
	    instanceLocations.add(-1);
	IUPAC = initialize_IUPAC();
    }



    /**************************************************************
     ********************** PUBLIC INSTANCE METHODS ***************
     **************************************************************/

    /**
     * Returns a <code>String</code> representation of the matrix that models the motif.
     *
     * @return   a <code>String</code> represenation of the matrix modeling the motif
     */
    public String matrixToString() {
	StringBuilder sb = new StringBuilder();
	DecimalFormat df = new DecimalFormat("0.00");  // Output decimal numbers to two decimal places
	for (int i=0; i<4; i++) {
	    for (int j=0; j<motifLength; j++) {
		sb.append(df.format(matrix[i][j]) + "  ");
	    }
	    sb.append("\n");
	}
	return sb.toString();
    }

    /**
     * Returns a <code>String</code> representation of motif instances found in the genomic sequences.
     *
     * @return   a <code>String</code> represenation of motif instances found in the genomic sequences
     */
    public String motifInstancesToString() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<numSequences; i++) {
	    int startIndexOfMotifInstance = instanceLocations.get(i);
	    sb.append("\t" + sequences.get(i).substring(startIndexOfMotifInstance, startIndexOfMotifInstance + motifLength) + "\n");
	}
	return sb.toString();
    }

    /**
     * Replaces values of 0.0 in the matrix (motif model) with small pseudocount values (0.0001).
     */
    public void addPseudocountsToMatrix() {
	double PSEUDOCOUNT = 0.0001;
	for (int i=0; i<4; i++) {
	    for (int j=0; j<motifLength; j++) {
		if (matrix[i][j] == 0.0)
		    matrix[i][j] = PSEUDOCOUNT;
	    }
	}
    }

    /**
     * Returns a <code>double</code> representing the frequency that the specified nucleotide character
     * occurs in the genomic sequences.
     *
     * @param   c   a nucleotide character (e.g., A or C or G or T)
     * @return   a <code>double</code> representing the frequency that <code>c</code> occurs in the sequences
     */
    public double getNucleotideContent(char c) {
	int count = 0;
	int total = 0;
	for (int i=0; i<numSequences; i++) {
	    for (int j=0; j<sequences.get(i).length(); j++) {
		if (sequences.get(i).charAt(j) == c)
		    count++;
	    }
	    total += sequences.get(i).length();
	}
	return ((double)count) / ((double)total);
    }

    /**
     * Returns a <em>copy</em> of the matrix that models the motif.
     *
     * @return   a 2D array represenation of the matrix modeling the motif
     */
    public double[][] getMatrix() {
	double[][] m = new double[4][motifLength];
	for (int i=0; i<4; i++) {
	    for (int j=0; j<motifLength; j++)
		m[i][j] = matrix[i][j];
	}
	return m;
    }

    /**
     * Returns a <em>copy</em> of the indices of start locations of motif instances in the sequences.
     *
     * @return   a collection of indices of start locations of motif instances in the sequences
     */
    public Vector<Integer> getInstanceLocations() {
	Vector<Integer> v = new Vector<Integer>(numSequences);
	for (int i=0; i<numSequences; i++) {
	    v.add(instanceLocations.get(i));
	}
	return v;
    }

    /**
     * Returns the consensus sequence for the motif.
     * <a href="http://www.bioinformatics.org/sms/iupac.html">IUPAC symbols</a>
     * are used as the alphabet for the consensus sequence.
     *
     * @return   a <code>String</code> representing the consensus sequence for the motif
     */
    public String getConsensusSequence() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<motifLength; i++) {
	    Vector<Double> column = new Vector<Double>(4);
	    for (int j=0; j<4; j++)
		column.add(matrix[j][i]);
	    sb.append(get_IUPAC(column));
	}
	return sb.toString();
    }



    /***************************************************************
     ********************** PRIVATE INSTANCE METHODS ***************
     ***************************************************************/

    /**
     * Returns the appropriate IUPAC symbol for the <code>Vector</code> of 4 nucleotide frequencies.
     *
     * @param   frequencies   a collection of 4 values, each representing the frequency of one of the 4 nucleotides
     * @return   an IUPAC character
     */
    private char get_IUPAC(Vector<Double> frequencies) {
	StringBuilder sb = new StringBuilder();
	double max = getMaximum(frequencies);
	if (frequencies.get(0) == max) sb.append("A");
	if (frequencies.get(1) == max) sb.append("C");
	if (frequencies.get(2) == max) sb.append("G");
	if (frequencies.get(3) == max) sb.append("T");
	return IUPAC.get(sb.toString());
    }

    /**
     * Returns the maximum value in the specified <code>Vector</code>.
     *
     * @param   v   a <code>Vector<code> or decimal numbers
     * @return   the maximum value in the specified <code>Vector</code>
     */
    private double getMaximum(Vector<Double> v) {
	double max = -1.0;
	for (int i=0; i<v.size(); i++) {
	    if (v.get(i) > max) max = v.get(i);
	}
	return max;
    }

    /**
     * Returns a dictionary of IUPAC symbols.
     *
     * @return   a <code>Hashtable</code> with keys corresponding to <code>Strings</code> of nucleotides and values corresponding to IUPAC characters
     */
    private Hashtable<String, Character> initialize_IUPAC() {
	Hashtable<String, Character> t = new Hashtable<String, Character>();
	t.put("A", 'A');    // adenine
	t.put("C", 'C');    // cytosine
	t.put("G", 'G');    // guanine
	t.put("T", 'T');    // thymine
	t.put("U", 'U');    // uracil
	t.put("AG", 'R');   // purine
	t.put("CT", 'Y');   // pyrimidine
	t.put("GT", 'K');   // keto
	t.put("AC", 'M');   // amino
	t.put("CG", 'S');
	t.put("AT", 'W');
	t.put("CGT", 'B');
	t.put("AGT", 'D');
	t.put("ACT", 'H');
	t.put("ACG", 'V');
	t.put("ACGT", 'N');
	return t;
    }



    /**************************************************************
     ********************** PRIVATE CLASS METHODS *****************
     **************************************************************/

    /**
     * Returns a collection of genomic sequences found in a FASTA file.
     *
     * @param   fileName   the name of a FASTA file containing one or more genomic sequences
     * @return   a collection of genomic sequences
     */
    private static Vector<String> readInMultipleSequencesFromFastaFile(String fileName) {
	Vector<String> seqs = new Vector<String>();
	try {
	    Scanner reader = new Scanner(new File(fileName));
	    StringBuilder sequence = new StringBuilder();
	    while (reader.hasNext()) {  // continue until we reach end of file
		String line = reader.nextLine();
		if ((line.length() == 0) || (line.charAt(0) == '>')) {
		    if (sequence.length() > 0)  // Add previous (completed) sequence to Vector
			seqs.add(processString(sequence.toString()));
		    sequence = new StringBuilder();  // Begin a new sequence
		} else {
		    sequence.append(line);  // Add line to growing sequence
		}
	    }
	    if (sequence.length() > 0)  // Add final sequence in file to Vector
		seqs.add(processString(sequence.toString()));
	    reader.close();
	} catch (FileNotFoundException e) {
	    System.err.println("Error - the file " + fileName + " could not be found and opened.");
	    System.exit(0);
	}
	return seqs;
    }

    /**
     * Returns a <em>cleaned</em> version of the specified input <code>String</code>.
     * A sequence is <em>clean</em> if it contains all capital letters and if it contains
     * only the 4 nucleotides: A or C or G or T. 
     *
     * @param   s   a <code>String</code> representing a genomic sequence
     * @return   a <em>cleaned</em> version of the input <code>String</code>
     *
     */
    private static String processString(String s) {
	String processed = s.toUpperCase();  // Convert to upper case
	processed = s.replace('U', 'T');  // Replace U with T

	// Check if String contains any non-nucleotide characters
	for (int i=0; i<s.length(); i++) {
	    if ((s.charAt(i) != 'A') && (s.charAt(i) != 'C') &&
		(s.charAt(i) != 'G') && (s.charAt(i) != 'T')) {
		System.out.println("Error - sequence contains non-nucleotide characters");
		System.out.println("\n" + s + "\n");
	    }
	}

	return processed;
    }



    /*********************************************************
     ****************** MAIN METHOD **************************
     *********************************************************/

    /**
     * The <code>main</code> method creates an initially empty <code>MotifSearch</code>
     * for a motif of the specified length in the genomic sequences found in the specified
     * FASTA file.
     * 
     * @param   args   an array of <code>Strings</code> representing any command line arguments
     */
    public static void main(String[] args) {

	if (args.length < 2) {
	    System.err.println("\nTo execute this program, two command line arguments are required corresponding to a FASTA file name and the length of the desired motif. For example,\n");
	    System.err.println("\tjava MotifSearching modA.txt 16\n\n");
	    System.exit(0);
	}

	MotifSearch m = new MotifSearch(args[0], Integer.parseInt(args[1]));

    }

}