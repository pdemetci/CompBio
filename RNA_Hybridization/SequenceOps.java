import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Random;

public class SequenceOps {

    /*********************************************************
     ****************** PRIVATE CLASS VARIABLES **************
     *********************************************************/
    
    /*
     * AAA, AAC, AAG, AAU, ACA, ACC, ACG, ACU,
     * AGA, AGC, AGG, AGU, AUA, AUC, AUG, AUU,
     * CAA, CAC, CAG, CAU, CCA, CCC, CCG, CCU,
     * CGA, CGC, CGG, CGU, CUA, CUC, CUG, CUU,
     * GAA, GAC, GAG, GAU, GCA, GCC, GCG, GCU,
     * GGA, GGC, GGG, GGU, GUA, GUC, GUG, GUU,
     * UAA, UAC, UAG, UAU, UCA, UCC, UCG, UCU,
     * UGA, UGC, UGG, UGU, UUA, UUC, UUG, UUU
     *
     *  K    N    K    N    T    T    T    T
     *  R    S    R    S    I    I    M    I
     *  Q    H    Q    H    P    P    P    P
     *  R    R    R    R    L    L    L    L
     *  E    D    E    D    A    A    A    A
     *  G    G    G    G    V    V    V    V
     *  *    Y    *    Y    S    S    S    S
     *  *    C    W    C    L    F    L    F
     */
    private static char[] codonTable = {'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T',
                                        'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I',
                                        'Q', 'H', 'Q', 'H', 'P', 'P', 'P', 'P',
                                        'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L',
                                        'E', 'D', 'E', 'D', 'A', 'A', 'A', 'A',
                                        'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V',
                                        '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S',
                                        '*', 'C', 'W', 'C', 'L', 'F', 'L', 'F'};



    /*********************************************************
     ****************** PUBLIC CLASS METHODS *****************
     *********************************************************/

    /* Returns a String corresponding to the genomic sequence found 
     * in the specified FASTA file.
     */
    public static String sequenceFromFastaFile(String fileName) {
	StringBuilder sequence = new StringBuilder();
	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String header = reader.nextLine();  // Header line of FASTA file
	    if ((header.length() == 0) || (header.charAt(0) != '>')) {
		System.err.println("Error - first line of file " + fileName + " is not in FASTA format.");
		return sequence.toString();
	    }
	    while (reader.hasNext()) {  // continue until we reach end of file
		sequence.append(reader.nextLine());
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
	    System.err.println("Error - the file " + fileName + " could not be found and opened.");
	    return sequence.toString();
	}
	return sequence.toString();
    }

    /* Takes the specified genomic sequence s and returns a String representing 
     * the sequence in FASTA format. Each line in the returned String should be 
     * 60 characters in length, except possibly the last line. The returned 
     * String need not include the FASTA header line beginning with the character '>'.
     */
    public static String sequenceToFastaFormat(String s) {
	int fastaLineLength = 60;
	StringBuilder sb = new StringBuilder();
	//sb.append("> FASTA sequence\n");
	int i = 0;
	while (i < s.length()) {
	    if (i+fastaLineLength > s.length())  // Last line of FASTA format
		sb.append(s.substring(i) + "\n");
	    else
		sb.append(s.substring(i, i+fastaLineLength) + "\n");
	    i += fastaLineLength;
	}
	return sb.toString();
    }

    /* Returns a String corresponding to the reversed version of the specified 
     * sequence s. For example, if s represents the sequence ACGGACTGC then 
     * the method should return the string CGTCAGGCA.
     */
    public static String reverse(String s) {
	// StringBuilders are more efficient than Strings or StringBuffers
	StringBuilder sb = new StringBuilder(s);
	return sb.reverse().toString();
    }

    /* Returns a String corresponding to the complement of the specified 
     * nucleotide sequence s. For example, if s represents the sequence 
     * ACGGACTGC then the method should return the string TGCCTGACG.
     */
    public static String complement(String s) {
	// StringBuilders are more efficient than Strings or StringBuffers
	StringBuilder sb = new StringBuilder(s);
	boolean isRNA = (sb.indexOf("u") >= 0) || (sb.indexOf("U") >= 0);
	for (int i=0; i<sb.length(); i++) {
	    if (sb.charAt(i) == 'a') {
		if (isRNA) sb.setCharAt(i, 'u');
		else sb.setCharAt(i, 't');
	    }
	    else if (sb.charAt(i) == 'c') sb.setCharAt(i, 'g');
	    else if (sb.charAt(i) == 'g') sb.setCharAt(i, 'c');
	    else if (sb.charAt(i) == 't') sb.setCharAt(i, 'a');
	    else if (sb.charAt(i) == 'u') sb.setCharAt(i, 'a');
	    else if (sb.charAt(i) == 'A') {
		if (isRNA) sb.setCharAt(i, 'U');
		else sb.setCharAt(i, 'T');
	    }
	    else if (sb.charAt(i) == 'C') sb.setCharAt(i, 'G');
	    else if (sb.charAt(i) == 'G') sb.setCharAt(i, 'C');
	    else if (sb.charAt(i) == 'T') sb.setCharAt(i, 'A');
	    else if (sb.charAt(i) == 'U') sb.setCharAt(i, 'A');
	}
	return sb.toString();
    }

    /* Returns a String corresponding to the reverse complement of the 
     * specified nucleotide sequence s. For example, if s represents 
     * the sequence ACGGACTGC then the method should return the string GCAGTCCGT.
     */
    public static String reverseComplement(String s) {
	return reverse(complement(s));
    }

    /* Returns the GC content of the specified nucleotide sequence s. 
     * For example, if s represents the sequence ACGGACTGC then the 
     * method should return 0.6666666666666666.
     */
    public static double GC_content(String s) {
	double GCs = 0.0;
	for (int i=0; i<s.length(); i++) {
	    if ((s.charAt(i) == 'g') || (s.charAt(i) == 'c') || (s.charAt(i) == 'G') || (s.charAt(i) == 'C'))
		GCs += 1.0;
	}
	return GCs / (double)(s.length());
    }

    /* Returns a random permutation of the specified sequence s. The 
     * returned String should have all the same characters as the input 
     * String s, only the order of the characters should be determined 
     * randomly. Each invocation of the method on a particular String s 
     * should result in a (almost certainly) different permutation of s.
     */
    public static String randomPermutation(String s) {
	// StringBuilders are more efficient than Strings or StringBuffers
	StringBuilder sb = new StringBuilder(s);
	Random r = new Random();
	int length = sb.length();
	for (int i=length-1; i>0; i--) {
	    int randomIndex = r.nextInt(i+1);  // Generates random number [0..i+1)
	    swap(sb, i, randomIndex); // Swap characters at indices "i" and "randomIndex"
	}
	return sb.toString();
    }

    /* Returns a random genomic sequence of the specified length with the 
     * expected specified GC content. Each nucleotide in the random 
     * sequence can be generated independently. The probability that each 
     * nucleotide is a G or C (as opposed to A or T) is specified by the 
     * GC_content parameter. Each nucleotide has an equal probability of 
     * being A or T. Each nucleotide has an equal probability of being 
     * G or C. Note: each invocation of the randomSequence method may not 
     * generate a random sequence with a GC content equal to that 
     * specified by the GC_content parameter. However, the expected GC 
     * content will be equal to the GC_content  parameter, i.e., if you 
     * invoked the method infinitely many times then the GC content of 
     * all randomly generated sequences should equal the GC_content parameter.
     */
    public static String randomSequence(int length, double GC_content) {
	// Calculate probability of each nucleotide
	double prob_A = (1.0 - GC_content) / 2.0;
	double prob_C = GC_content / 2.0;
	double prob_G = GC_content / 2.0;
	double prob_T = (1.0 - GC_content) / 2.0;

	// Calulate range for each nucleotide between 0.0 and 1.0
	prob_A = prob_A;  // [0.0, prob_A)
	prob_C = prob_A + prob_C;  // [prob_A, prob_A + prob_C)
	prob_G = prob_C + prob_G;  // [prob_C, prob_C + prob_G)
	prob_T = prob_G + prob_T;  // [prob_G, prob_G + prob_T)

	StringBuilder sb = new StringBuilder();
	for (int i=0; i<length; i++) {
	    double r = Math.random();
	    if (r < prob_A) sb.append("A");
	    else if (r < prob_C) sb.append("C");
	    else if (r < prob_G) sb.append("G");
	    else sb.append("T");
	}
	return sb.toString();
    }

    /* Returns a random genomic sequence with the same length as s and 
     * the same expected GC content as s. Note: each invocation of the 
     * randomSampling method may not generate a random sequence with a 
     * GC content equal to that of s. However, the expected GC content 
     * will be equal to that of s, i.e., if you invoked the method 
     * infinitely many times then the GC content of all randomly 
     * generated sequences should equal that of s.
     */
    public static String randomSampling(String s) {
	return randomSequence(s.length(), GC_content(s));
    }

    /* Assuming the input sequence s is a genomic sequence of exactly 3 
     * nucleotides, the method translates the codon, i.e., it returns a 
     * character representing the amino acid corresponding to the codon. 
     * For example, if the input sequence s is CUA then the method should 
     * return the character L corresponding to the amino acid Leucine. 
     * Note: the expected running time of this method should not be linear 
     * in the size of the table, i.e., the method should not compare each 
     * of the 64 codons. Instead, the expected running time of the method 
     * should be constant. Consider what data structures would be appropriate 
     * to store the table in order to support a constant expected running time 
     * of this method.
     */
    public static char translateCodon(String s) {
	char aa = '?';
	if (s.length() != 3) {
	    System.err.println("Cannot translate codon because its length is not 3.");
	    return aa;
	}
	int index = nucleotideToNumber(s.charAt(0))*16 + nucleotideToNumber(s.charAt(1))*4 + nucleotideToNumber(s.charAt(2))*1;
	if ((index >=0) && (index < 64)) aa = codonTable[index];
	return aa;
    }

    /* Assuming the length of the input sequence s is divisible by 3, the 
     * method returns a translated version of the sequence, i.e., every 3 
     * nucleotides in s are translated into an amino acid, and the sequence 
     * of amino acids is returned. For example, if s represents the 
     * sequence AUGGCCUUUCGAUAG then the method should return the string MAFR*.
     */
    public static String translateORF(String s) {
	if (s.length() % 3 != 0) {  // Length of sequence must be divisible by 3
	    System.err.println("Cannot translate sequence because its length is not divisable by 3 and therefore the sequence is not an ORF.");
	    return "";
	}
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<s.length(); i+=3) {
	    String codon = s.substring(i, i+3);
	    sb.append(translateCodon(codon));
	}
	return sb.toString();
    }



    /*********************************************************
     ****************** PRIVATE CLASS METHODS ****************
     *********************************************************/

    // Swaps characters at indices i and j in given StringBuilder
    private static void swap(StringBuilder sb, int i, int j) {
	char temp = sb.charAt(i);
	sb.setCharAt(i, sb.charAt(j));
	sb.setCharAt(j, temp);
    }

    // Convert nucleotide to number: A=0, C=1, G=2, T/U=3, OTHER=100
    private static int nucleotideToNumber(char nt) {
	if ((nt == 'a') || (nt == 'A')) return 0;
	if ((nt == 'c') || (nt == 'C')) return 1;
	if ((nt == 'g') || (nt == 'G')) return 2;
	if ((nt == 't') || (nt == 'T')) return 3;
	if ((nt == 'u') || (nt == 'U')) return 3;
	return 100;
    }

    // Tests public class methods
    private static void testMethodsOnSequence(String s) {
	System.out.println();
	System.out.println("Sequence:\t" + s);
	System.out.println("Reverse:\t" + reverse(s));
	System.out.println("Complement:\t" + complement(s));
	System.out.println("RevComp:\t" + reverseComplement(s));
	System.out.println("GC content:\t" + GC_content(s));
	System.out.println("Translation:\t" + translateORF(s));

	// Test random permutation
	System.out.println();
	for (int i=0; i<3; i++) {
	    String rp = randomPermutation(s);
	    System.out.println("Rand Perm " + (i+1) + ":\t" + rp);
	    System.out.println("GC of RandPerm:\t" + GC_content(rp));
	}

	// Test random sampling
	System.out.println();
	for (int i=0; i<3; i++) {
	    String rs = randomSampling(s);
	    System.out.println("Rand Samp " + (i+1) + ":\t" + rs);
	    System.out.println("GC of RandSamp:\t" + GC_content(rs));
	}
	System.out.println();
    }



    /*********************************************************
     ****************** MAIN METHOD **************************
     *********************************************************/

    public static void main(String [] args) {
	testMethodsOnSequence("ATGTCACCGCAGCTAGACTcgattcgaggcaagctgtaa");
	testMethodsOnSequence("...................................G.....................");
    }

}

