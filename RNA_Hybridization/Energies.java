import java.util.*;
import java.io.*;

/*****************************************************
 * An instance of the <code>Energies</code> class represents
 * the thermodynamic energies of various RNA secondary
 * structures.
 *****************************************************/
public class Energies {

    /*********************************************************
     ****************** INSTANCE VARIABLES *******************
     *********************************************************/

    private String dir;
    private String sequence1;
    private String sequence2;
    private Vector<Double> interiorLoopLengths;
    private Vector<Double> bulgeLoopLengths;
    private Vector<Double> hairpinLoopLengths;
    private Hashtable<String, Double> hairpinTerminals;
    private Hashtable<String, Double> hairpinBonuses;
    private Hashtable<String, Double> stacking;
    private Hashtable<String, Double> interiorTerminals;
    private Hashtable<String, Double> interior_1_1;
    private Hashtable<String, Double> interior_1_2;
    private Hashtable<String, Double> interior_2_2;
    


    /***************************************************
     ****************** CONSTRUCTORS *******************
     ***************************************************/

    /**
     * Creates an <code>Energies</code> object that contains
     * information about the thermodynamic energies of various
     * RNA secondary structures that may form as a result of
     * hybridization (i.e., <em>inter</em>molecular basepairing)
     * between the two specified RNA nucleotide sequences.
     *
     * @param   dir         the folder where the files of energy data reside
     * @param   sequence1   the sequence of nucleotides in an RNA molecule
     * @param   sequence2   the sequence of nucleotides in an RNA molecule
     */
    public Energies(String dir, String sequence1, String sequence2) {
	if (dir.charAt(dir.length()-1) != '/') dir = dir + "/";
	this.dir = dir;
	this.sequence1 = sequence1;
	this.sequence2 = sequence2;
	readInEnergyData();
    }

    /**
     * Creates an <code>Energies</code> object that contains
     * information about the thermodynamic energies of various
     * RNA secondary structures that may form as a result of
     * folding (i.e., <em>intra</em>molecular basepairing)
     * the specified RNA nucleotide sequence.
     *
     * @param   dir        the folder where the files of energy data reside
     * @param   sequence   the sequence of nucleotides in an RNA molecule
     */
    public Energies(String dir, String sequence) {
	this(dir, sequence, sequence);
    }



    /*******************************************************
     ************** PUBLIC INSTANCE METHIDS ****************
     *******************************************************/

    /**
     * Returns the energy of a hairpin loop whose closing basepair
     * involves the two nucleotides at indices <code>i</code>
     * and <code>j</code>.
     *
     * @param   i   index of one of the nucleotides in the hairpin loop's closing basepair
     * @param   j   index of one of the nucleotides in the hairpin loop's closing basepair
     * @return      the energy of the hairpin loop
     */
    public double hairpin(int i, int j) {
 
	// Ensure loop has closing basepair
	if (!isPair(sequence1.charAt(i), sequence1.charAt(j))) return Double.POSITIVE_INFINITY;
	
	// Get contribution to energy of loop length
	double lengthEnergy = Double.POSITIVE_INFINITY;
	int length = j - i - 1;
	if (length <= 0) return Double.POSITIVE_INFINITY;
	if (length < hairpinLoopLengths.size())
	    lengthEnergy = hairpinLoopLengths.get(length);
	else
	    lengthEnergy = hairpinLoopLengths.get(hairpinLoopLengths.size()-1) + (1.079*Math.log(((double)length)/(double)(hairpinLoopLengths.size()-1)));
	
	// Get terminal stacking pairs contribution
	double terminalEnergy = 0.0;;
	String s = sequence1.charAt(i) + "" + sequence1.charAt(j) + "" + sequence1.charAt(i+1) + "" + sequence1.charAt(j-1);
	if (length > 3)  // Only apply for loops with length greater than 3
	    terminalEnergy = hairpinTerminals.get(s);
	
	// Get bonus energy contribution
	double bonusEnergy = 0.0;
	if (hairpinBonuses.containsKey(sequence1.substring(i, j+1)))
	    bonusEnergy = hairpinBonuses.get(sequence1.substring(i, j+1));
	
	// Get miscellaneous energy contribution
	double miscellaneousEnergy = 0.0;  // Ignore these contributions (miscloop.dat file)
	if ((i-2 >= 0) && (sequence1.charAt(i-2) == 'G') && (sequence1.charAt(i-1) == 'G') &&
	    (sequence1.charAt(i) == 'G') && (sequence1.charAt(j) == 'U'))  // GGG hairpin
	    miscellaneousEnergy += -2.20;
	String ss_loop = sequence1.substring(i+1, j);  // Single stranded loop sequence
	if (ss_loop.replaceAll("C", "").length() == 0) {  // Single strand of loop contains all C's
	    if (length == 3) miscellaneousEnergy += 1.40;
	    else miscellaneousEnergy += 0.30*length + 1.60;
	}
	return lengthEnergy + terminalEnergy + bonusEnergy + miscellaneousEnergy;
    }

    /**
     * Returns the energy of the stacking region defined by two pairs of 
     * basepairing nucleotides. One pair of basepairing nucleotides
     * corresponds to nucleotides at indices <code>i</code> and
     * <code>j</code>. The other pair of basepairing nucleotides
     * corresponds to nucleotides at indices <code>x</code> and
     * <code>y</code>. Note that <code>i</code> and <code>x</code>
     * must be consecutive indices in a sequence and that
     * <code>j</code> and <code>y</code> must be consecutive 
     * indices in a sequence, or else the structure would not be a 
     * <em>stacking region</em> but rather a <em>bulge loop</em> or 
     * <em>interior loop</em>.
     *
     * @param   i   index of a nucleotide participating in the first basepairing
     * @param   j   index of a nucleotide participating in the first basepairing
     * @param   x   index of a nucleotide participating in the second basepairing
     * @param   y   index of a nucleotide participating in the second basepairing
     * @return      the energy of the stacking region
     */
    public double stacking(int i, int j, int x, int y) {
	String s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y);
	if (stacking.containsKey(s)) {
	    return stacking.get(s);
	}
	return Double.POSITIVE_INFINITY;
    }

    /**
     * Returns the energy of a bulge loop structure defined by two pairs of 
     * basepairing nucleotides. One pair of basepairing nucleotides
     * corresponds to nucleotides at indices <code>i</code> and
     * <code>j</code>. The other pair of basepairing nucleotides
     * corresponds to nucleotides at indices <code>x</code> and
     * <code>y</code>. Note that, for the structure to be a
     * <em>bulge loop</em>, either <code>i</code> and <code>x</code>
     * must be consecutive indices in a sequence or 
     * <code>j</code> and <code>y</code> must be consecutive 
     * indices in a sequence, but not both. If neither pair correspond to consecutive
     * indices then the structure is an <em>interior loop</em>. If both pairs
     * correspond to consecutive indices then the structure is a <em>stacking 
     * region</em>.
     *
     * @param   i   index of a nucleotide participating in the first basepairing
     * @param   j   index of a nucleotide participating in the first basepairing
     * @param   x   index of a nucleotide participating in the second basepairing
     * @param   y   index of a nucleotide participating in the second basepairing
     * @return      the energy of the bulge loop structure
     */
    public double bulge(int i, int j, int x, int y) {
	int minBulgeLength = x - i - 1;
	int maxBulgeLength = Math.max(j - y - 1, y - j - 1);  // One sequence or two
	if (minBulgeLength > maxBulgeLength) {
	    int temp = minBulgeLength;
	    minBulgeLength = maxBulgeLength;
	    maxBulgeLength = temp;
	}
	
	if ((minBulgeLength != 0) || (maxBulgeLength <= 0)) {
	    //System.err.println("Error - not a bulge loop.");
	    return Double.POSITIVE_INFINITY;
	}
	
	double lengthEnergy = Double.POSITIVE_INFINITY;
	if (maxBulgeLength < bulgeLoopLengths.size()) 
	    lengthEnergy = bulgeLoopLengths.get(maxBulgeLength);
	else
	    lengthEnergy = bulgeLoopLengths.get(bulgeLoopLengths.size()-1) + (1.079*Math.log(((double)maxBulgeLength)/(double)(bulgeLoopLengths.size()-1)));

	double stackingEnergy = 0.0;  // include stacking energy only for bulge of size 1
	if (maxBulgeLength == 1) {
	    String s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y);
	    stackingEnergy = stacking.get(s);
	}

	// Miscellaneous energy
	double miscellaneousEnergy = 0.0;
	if (((maxBulgeLength > 1) && (sequence1.charAt(i) == 'A') && (sequence2.charAt(j) == 'U')) ||
	    ((maxBulgeLength > 1) && (sequence1.charAt(i) == 'U') && (sequence2.charAt(j) == 'A')))
	    miscellaneousEnergy = 0.5;
	return lengthEnergy + stackingEnergy + miscellaneousEnergy;
    }

    /**
     * Returns the energy of an interior loop structure defined by two pairs of 
     * basepairing nucleotides. One pair of basepairing nucleotides
     * corresponds to nucleotides at indices <code>i</code> and
     * <code>j</code>. The other pair of basepairing nucleotides
     * corresponds to nucleotides at indices <code>x</code> and
     * <code>y</code>. Note that, for the structure to be an
     * <em>interior loop</em>, <code>i</code> and <code>x</code>
     * may not be consecutive indices in a sequence and  
     * <code>j</code> and <code>y</code> may not be consecutive 
     * indices in a sequence. If one pair corresponds to consecutive
     * indices then the structure is a <em>bulge loop</em>. If both pairs 
     * correspond to consecutive indices then the structure is a 
     * <em>stacking region</em>.
     *
     * @param   i   index of a nucleotide participating in the first basepairing
     * @param   j   index of a nucleotide participating in the first basepairing
     * @param   x   index of a nucleotide participating in the second basepairing
     * @param   y   index of a nucleotide participating in the second basepairing
     * @return      the energy of the interior loop structure
     */
    public double interior(int i, int j, int x, int y) {
	int minInteriorLength = x - i - 1;
	int maxInteriorLength = Math.max(j - y - 1, y - j - 1);  // One sequence or two
	if (minInteriorLength > maxInteriorLength) {
	    int temp = minInteriorLength;
	    minInteriorLength = maxInteriorLength;
	    maxInteriorLength = temp;
	}

	if (minInteriorLength == 0) {
	    //System.err.println("Error - not an interior loop.");
	    return Double.POSITIVE_INFINITY;
	}

	// 1x1 interior loop (special case)
	if ((minInteriorLength == 1) && (maxInteriorLength == 1)) {
	    String s = "";
	    if (sequence1.equals(sequence2))  // We have one sequence
		s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j-1);
	    else  // We have two sequences
		s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j+1);
	    if (interior_1_1.containsKey(s)) return interior_1_1.get(s);
	    else return Double.POSITIVE_INFINITY;
	}

	// 1x2 interior loop (special case)
	if ((minInteriorLength == 1) && (maxInteriorLength == 2)) {

	    // Assume, initially, that we have a 1x2 loop as opposed to a 2x1 loop
	    String s = "";
	    if (sequence1.equals(sequence2)) {  // We have one sequence
		s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j-1) + "" + sequence2.charAt(j-2);
		if (x-i > j-y)  // we have a 2x1 loop
		    s = sequence2.charAt(y) + "" + sequence1.charAt(x) + "" + sequence2.charAt(j) + "" + sequence1.charAt(i) + "" + sequence2.charAt(j-1) + "" + sequence1.charAt(x-1) + "" + sequence1.charAt(x-2);
	    } else {  // We have two sequences
		s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j+1) + "" + sequence2.charAt(j+2);
		if (x-i > y-j)  // we have a 2x1 loop
		    s = sequence2.charAt(y) + "" + sequence1.charAt(x) + "" + sequence2.charAt(j) + "" + sequence1.charAt(i) + "" + sequence2.charAt(j+1) + "" + sequence1.charAt(x-1) + "" + sequence1.charAt(x-2);
	    }
	    if (interior_1_2.containsKey(s)) return interior_1_2.get(s);
	    else return Double.POSITIVE_INFINITY;
	}

	// 2x2 interior loop (special case)
	if ((minInteriorLength == 2) && (maxInteriorLength == 2)) {

	    String s = "";
	    if (sequence1.equals(sequence2)) {  // We have one sequence
		s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j-1) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j-2);
	    } else {  // We have two sequences
		s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j+1) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j+2);
	    }
	    if (interior_2_2.containsKey(s)) return interior_2_2.get(s);
	    else return Double.POSITIVE_INFINITY;
	}

	// Length energy
	int loopLength = maxInteriorLength + minInteriorLength;
	double lengthEnergy = Double.POSITIVE_INFINITY;
	if (loopLength < interiorLoopLengths.size()) 
	    lengthEnergy = interiorLoopLengths.get(loopLength);
	else
	    lengthEnergy = interiorLoopLengths.get(interiorLoopLengths.size()-1) + (1.079*Math.log(((double)loopLength)/(double)(interiorLoopLengths.size()-1)));

	// Terminal mismatch energies
	double terminalEnergy = 0.0;
	String s = "";
	if (sequence1.equals(sequence2))  // We have one sequence
	    s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j-1);
	else  // We have two sequences
	    s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + sequence1.charAt(i+1) + "" + sequence2.charAt(j+1);
	if ((minInteriorLength == 1) && (maxInteriorLength > 2))  // GAIL rule (substitute AA mismatch)
	    s = sequence1.charAt(i) + "" + sequence2.charAt(j) + "" + "A" + "" + "A";
	terminalEnergy += interiorTerminals.get(s);
	if (sequence1.equals(sequence2))  // We have one sequence
	    s = sequence2.charAt(y) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y+1) + "" + sequence1.charAt(x-1);
	else  // We have two sequences
	    s = sequence2.charAt(y) + "" + sequence1.charAt(x) + "" + sequence2.charAt(y-1) + "" + sequence1.charAt(x-1);
	if ((minInteriorLength == 1) && (maxInteriorLength > 2))  // GAIL rule (substitute AA mismatch)
	    s = sequence2.charAt(y) + "" + sequence1.charAt(x) + "" + "A" + "" + "A";
	terminalEnergy += interiorTerminals.get(s);
      
	// Asymmetry energy
	// Asymmetry energy with branches N1 and N2 is the minumum of 3.0 or N*f(M),
	// where N = |N1-N2|, M is the minimum of 4, N1, and N2,
	// and f(1) = 0.4, f(2) = 0.3, f(3) = 0.2, and f(4) = 0.1
	// Taken from Jaeger, Turner, and Zuker (PNAS 1989)
	double asymmetryEnergy = 0.0;
	int N = maxInteriorLength - minInteriorLength;
	int M = Math.min(4, minInteriorLength);
	if (M == 1) asymmetryEnergy = N*0.4;
	if (M == 2) asymmetryEnergy = N*0.3;
	if (M == 3) asymmetryEnergy = N*0.2;
	if (M == 4) asymmetryEnergy = N*0.1;
	asymmetryEnergy = Math.min(asymmetryEnergy, 3.0);
      
	return lengthEnergy + terminalEnergy + asymmetryEnergy;
    }
    
    

    /*******************************************************
     ************** PRIVATE INSTANCE METHIDS ***************
     *******************************************************/

    // Reads in energy data from files
    private void readInEnergyData() {
	readInLoopLengths(dir + "loop.dat");
	readInTerminalHairpin(dir + "tstackh.dat");
	readInBonusHairpin(dir + "tloop.dat");
	readInStacking(dir + "stack.dat");
	readInTerminalInterior(dir + "tstacki.dat");
	readInInterior_1_1(dir + "sint2.dat");
	readInInterior_1_2(dir + "asint1x2.dat");
	readInInterior_2_2(dir + "sint4.dat");
    }

    private void readInLoopLengths(String fileName) {
	interiorLoopLengths = new Vector<Double>();
	interiorLoopLengths.add(Double.POSITIVE_INFINITY);
	bulgeLoopLengths = new Vector<Double>();
	bulgeLoopLengths.add(Double.POSITIVE_INFINITY);
	hairpinLoopLengths = new Vector<Double>();
	hairpinLoopLengths.add(Double.POSITIVE_INFINITY);

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = reader.nextLine();  // Ignore header line
	    line = reader.nextLine();  // Ignore header line
	    line = reader.nextLine();  // Ignore header line
	    line = reader.nextLine();  // Ignore header line
	    while (reader.hasNextLine()) {  // continue until we reach end of file
		String[] splitLine = reader.nextLine().split("\\s+");
		if (splitLine.length != 4) {
		    System.err.println("Error - expecting 4 tokens per line in file " + fileName);
		    System.exit(0);
		}
		if (splitLine[1].charAt(0) == '.') interiorLoopLengths.add(Double.POSITIVE_INFINITY);
		else interiorLoopLengths.add(Double.parseDouble(splitLine[1]));
		if (splitLine[2].charAt(0) == '.') bulgeLoopLengths.add(Double.POSITIVE_INFINITY);
		else bulgeLoopLengths.add(Double.parseDouble(splitLine[2]));
		if (splitLine[3].charAt(0) == '.') hairpinLoopLengths.add(Double.POSITIVE_INFINITY);
		else hairpinLoopLengths.add(Double.parseDouble(splitLine[3]));
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }

    private void readInTerminalHairpin(String fileName) {
	hairpinTerminals = new Hashtable<String, Double>();

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = "";
	    for (int i=0; i<17; i++)  // Ignore 17 header lines
		line = reader.nextLine();

	    for (int i=0; i<4; i++) {  // Each of the 4 rows of 4x4 matrices

		for (int j=0; j<9; j++)  // Ignore 9 header lines
		    line = reader.nextLine();

		for (int k=0; k<4; k++) {  // Each row of 4x4 matrices
		    String[] splitLine = reader.nextLine().trim().split("\\s+");

		    for (int m=0; m<4; m++) {  // Set of 4 values on a given row

			for (int n=0; n<4; n++) {  // One particular value in a matrix

			    String sequence = "" + intToNucleotide(i) + "" + intToNucleotide(m) + "" + intToNucleotide(k) + "" + intToNucleotide(n);
			    double energy = Double.POSITIVE_INFINITY;
			    if (splitLine[4*m + n].charAt(0) != '.')
				energy = Double.parseDouble(splitLine[4*m + n]);
			    hairpinTerminals.put(sequence, energy);
			}
		    }
		}
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }

    private void readInBonusHairpin(String fileName) {
	hairpinBonuses = new Hashtable<String, Double>();

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = "";
	    for (int i=0; i<2; i++)  // Ignore 2 header lines
		line = reader.nextLine();

	    while (reader.hasNextLine()) {  // continue until we reach end of file
		String[] splitLine = reader.nextLine().trim().split("\\s+");
		hairpinBonuses.put(splitLine[0], Double.parseDouble(splitLine[1]));
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }

    private void readInStacking(String fileName) {
	stacking = new Hashtable<String, Double>();

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = "";
	    for (int i=0; i<17; i++)  // Ignore 17 header lines
		line = reader.nextLine();

	    for (int i=0; i<4; i++) {  // Each of the 4 rows of 4x4 matrices

		for (int j=0; j<9; j++)  // Ignore 9 header lines
		    line = reader.nextLine();

		for (int k=0; k<4; k++) {  // Each row of 4x4 matrices
		    String[] splitLine = reader.nextLine().trim().split("\\s+");

		    for (int m=0; m<4; m++) {  // Set of 4 values on a given row

			for (int n=0; n<4; n++) {  // One particular value in a matrix

			    String sequence = "" + intToNucleotide(i) + "" + intToNucleotide(m) + "" + intToNucleotide(k) + "" + intToNucleotide(n);
			    double energy = Double.POSITIVE_INFINITY;
			    if (splitLine[4*m + n].charAt(0) != '.')
				energy = Double.parseDouble(splitLine[4*m + n]);
			    stacking.put(sequence, energy);
			}
		    }
		}
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }

    private void readInTerminalInterior(String fileName) {
	interiorTerminals = new Hashtable<String, Double>();

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = "";
	    for (int i=0; i<17; i++)  // Ignore 17 header lines
		line = reader.nextLine();

	    for (int i=0; i<4; i++) {  // Each of the 4 rows of 4x4 matrices

		for (int j=0; j<9; j++)  // Ignore 9 header lines
		    line = reader.nextLine();

		for (int k=0; k<4; k++) {  // Each row of 4x4 matrices
		    String[] splitLine = reader.nextLine().trim().split("\\s+");

		    for (int m=0; m<4; m++) {  // Set of 4 values on a given row

			for (int n=0; n<4; n++) {  // One particular value in a matrix

			    String sequence = "" + intToNucleotide(i) + "" + intToNucleotide(m) + "" + intToNucleotide(k) + "" + intToNucleotide(n);
			    double energy = Double.POSITIVE_INFINITY;
			    if (splitLine[4*m + n].charAt(0) != '.')
				energy = Double.parseDouble(splitLine[4*m + n]);
			    interiorTerminals.put(sequence, energy);
			}
		    }
		}
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }

    // A 1x1 interior loop with the form
    //                        GXA
    //                        | |
    //                        CYU
    //
    // would have the sequence GCAUXY
    private void readInInterior_1_1(String fileName) {
	interior_1_1 = new Hashtable<String, Double>();

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = "";
	    for (int i=0; i<18; i++)  // Ignore 18 header lines
		line = reader.nextLine();

	    for (int i=0; i<6; i++) {  // Each of the 6 rows of 4x4 matrices

		for (int j=0; j<11; j++)  // Ignore 11 header lines
		    line = reader.nextLine();

		for (int k=0; k<4; k++) {  // Each row of 4x4 matrices
		    String[] splitLine = reader.nextLine().trim().split("\\s+");

		    for (int m=0; m<6; m++) {  // Set of 6 values on a given row

			for (int n=0; n<4; n++) {  // One particular value in a matrix

			    String sequence = "" + intToPair(i) + "" + intToPair(m) + "" + intToNucleotide(k) + "" + intToNucleotide(n);
			    double energy = Double.POSITIVE_INFINITY;
			    if (splitLine[4*m + n].charAt(0) != '.')
				energy = Double.parseDouble(splitLine[4*m + n]);
			    interior_1_1.put(sequence, energy);
			}
		    }
		}
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }

    // A 1x2 interior loop with the form
    //                        GX A
    //                        |  |
    //                        CYAU
    //
    // would have the sequence GCAUXYA
    private void readInInterior_1_2(String fileName) {
	interior_1_2 = new Hashtable<String, Double>();

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = "";
	    for (int i=0; i<18; i++)  // Ignore 18 header lines
		line = reader.nextLine();

	    for (int i=0; i<24; i++) {  // Each of the 24 rows of 4x4 matrices

		for (int j=0; j<11; j++)  // Ignore 11 header lines
		    line = reader.nextLine();

		for (int k=0; k<4; k++) {  // Each row of 4x4 matrices
		    String[] splitLine = reader.nextLine().trim().split("\\s+");

		    for (int m=0; m<6; m++) {  // Set of 6 values on a given row

			for (int n=0; n<4; n++) {  // One particular value in a matrix

			    String sequence = "" + intToPair(i/4) + "" + intToPair(m) + "" + intToNucleotide(k) + "" + intToNucleotide(n) + "" + intToNucleotide(i%4);
			    double energy = Double.POSITIVE_INFINITY;
			    if (splitLine[4*m + n].charAt(0) != '.')
				energy = Double.parseDouble(splitLine[4*m + n]);
			    interior_1_2.put(sequence, energy);
			}
		    }
		}
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }

    // A 2x2 interior loop with the form
    //                        GPRA
    //                        |  |
    //                        CQSU
    //
    // would have the sequence GCAUPQRS
    private void readInInterior_2_2(String fileName) {
	interior_2_2 = new Hashtable<String, Double>();

	try {
	    Scanner reader = new Scanner(new File(fileName));
	    String line = "";
	    for (int i=0; i<31; i++)  // Ignore 31 header lines
		line = reader.nextLine();

	    for (int i=0; i<36; i++) {  // Each of the 36 rows of 16x16 matrices

		for (int j=0; j<10; j++)  // Ignore 10 header lines
		    line = reader.nextLine();

		for (int k=0; k<16; k++) {  // Each row of 16x16 matrices
		    String[] splitLine = reader.nextLine().trim().split("\\s+");

		    for (int n=0; n<16; n++) {  // One particular value in a matrix

			String sequence = "" + intToPair(i/6) + "" + intToPair(i%6) + "" + intToNucleotide(k/4) + "" + intToNucleotide(k%4) + "" + intToNucleotide(n/4) + "" + intToNucleotide(n%4);
			double energy = Double.POSITIVE_INFINITY;
			if (splitLine[n].charAt(0) != '.')
			    energy = Double.parseDouble(splitLine[n]);
			interior_2_2.put(sequence, energy);
		    }
		}
	    }
	    reader.close();
	} catch (FileNotFoundException e) {
            System.err.println("Error - the file " + fileName + " could not be found and opened.");
            System.exit(0);
	}
    }
    
    // Do we have a canonical or wobble basepair?
    private boolean isPair(char a, char b) {
	if ((a == 'G') && (b == 'C')) return true;
	if ((a == 'G') && (b == 'U')) return true;
	if ((a == 'C') && (b == 'G')) return true;
	if ((a == 'A') && (b == 'U')) return true;
	if ((a == 'U') && (b == 'A')) return true;
	if ((a == 'U') && (b == 'G')) return true;
	return false;
    }

    // Convert an integer to a nucleotide: 0=A, 1=C, 2=G, 3=U
    private char intToNucleotide(int index) {
	if (index == 0) return 'A';
	if (index == 1) return 'C';
	if (index == 2) return 'G';
	return 'U';
    }

    // Convert an integer to a basepair: 0=AU, 1=CG, 2=GC, 3=UA, 4=GU, 5=UG
    private String intToPair(int index) {
	if (index == 0) return "AU";
	if (index == 1) return "CG";
	if (index == 2) return "GC";
	if (index == 3) return "UA";
	if (index == 4) return "GU";
	return "UG";
    }



    // Main method
    public static void main(String[] args) {

        System.err.println("\nTo use this class, create an Energies object based on one sequence (to compute energies of secondary structures in a folded RNA sequence) or based on two sequences (to compute energies of secondary structures in two hybridizing RNA sequences). Then invoke one of the four instance methods:");
        System.err.println("\thairpin(int, int)");
        System.err.println("\tstacking(int, int, int, int)");
        System.err.println("\tbulge(int, int, int, int)");
        System.err.println("\tinterior(int, int, int, int)\n");


	/* TESTING CODE

	Energies e = new Energies("energyData", "UACGUAGCUUAGUUACGUACGU");

	// Test loop lengths
	System.out.println("Energy of hairpin loop of length 3 is:\t" + e.hairpinLoopLengths.get(3));
	System.out.println("Energy of hairpin loop of length 10 is:\t" + e.hairpinLoopLengths.get(10));
	System.out.println("Energy of hairpin loop of length 30 is:\t" + e.hairpinLoopLengths.get(30));
	System.out.println("Energy of interior loop of length 3 is:\t" + e.interiorLoopLengths.get(3));
	System.out.println("Energy of interior loop of length 10 is:\t" + e.interiorLoopLengths.get(10));
	System.out.println("Energy of interior loop of length 30 is:\t" + e.interiorLoopLengths.get(30));
	System.out.println("Energy of bulge loop of length 3 is:\t" + e.bulgeLoopLengths.get(3));
	System.out.println("Energy of bulge loop of length 10 is:\t" + e.bulgeLoopLengths.get(10));
	System.out.println("Energy of bulge loop of length 30 is:\t" + e.bulgeLoopLengths.get(30));
	System.out.println();

	// Test hairpin terminal nucleotides
	System.out.println("Energy of terminal hairpin nucleotides G:C and A:A is:\t" + e.hairpinTerminals.get("GCAA"));
	System.out.println("Energy of terminal hairpin nucleotides A:A and A:U is:\t" + e.hairpinTerminals.get("AAAU"));
	System.out.println("Energy of terminal hairpin nucleotides A:U and U:A is:\t" + e.hairpinTerminals.get("AUUA"));
	System.out.println("Energy of terminal hairpin nucleotides G:C and G:C is:\t" + e.hairpinTerminals.get("GCGC"));
	System.out.println();
 
	// Test hairpin bonuses
	System.out.println("Energy of hairpin bonus GGGGAC is:\t" + e.hairpinBonuses.get("GGGGAC"));
	System.out.println("Energy of hairpin bonus AGAAAU is:\t" + e.hairpinBonuses.get("AGAAAU"));
	System.out.println("Energy of hairpin bonus UGGAAA is:\t" + e.hairpinBonuses.get("UGGAAA"));
	System.out.println();
 
	// Test stacking basepairs
	System.out.println("Energy of stacking pairs G:C and A:A is:\t" + e.stacking.get("GCAA"));
	System.out.println("Energy of stacking pairs A:A and A:U is:\t" + e.stacking.get("AAAU"));
	System.out.println("Energy of stacking pairs A:U and U:A is:\t" + e.stacking.get("AUUA"));
	System.out.println("Energy of stacking pairs G:C and G:C is:\t" + e.stacking.get("GCGC"));
	System.out.println();
 
	// Test bulge loops
	System.out.println("Energy of bulge loop UA and GC is:\t" + e.bulge(0,3,1,2));
	System.out.println("Energy of bulge loop CG and GC is:\t" + e.bulge(2,16,3,15));
	System.out.println("Energy of bulge loop AGC and UG is:\t" + e.bulge(5,17,7,16));
	System.out.println("Energy of bulge loop AG and UGC is:\t" + e.bulge(5,17,6,15));
	System.out.println("Energy of bulge loop UACGUA and AU is:\t" + e.bulge(0,10,5,9));
	System.out.println("Energy of bulge loop UA and AUUCGAU is:\t" + e.bulge(0,10,1,4));
	System.out.println("Energy of bulge loop UAC and GUA is:\t" + e.bulge(0,18,2,16));
	*/
    }

}
