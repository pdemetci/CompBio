public class RNA_Energy {

    // Instance variables
    private String sequence;
    private Energies e;
    private boolean energyDetails;  // Output energy details for debugging
    private double[][] V;
    private double[][] W;



    // Constructor
    public RNA_Energy(String fileName) {
	sequence = SequenceOps.sequenceFromFastaFile(fileName).toUpperCase().replace('T', 'U');
	e = new Energies("energyData", sequence);
	V = new double[sequence.length()][sequence.length()];
	W = new double[sequence.length()][sequence.length()];
    }



    // Public instance methods

    // For the RNA sequence, fill in the dynamic programming tables, V and W, to compute 
    // the structure that has the lowest energy
    public void computeStructureWithMinimumEnergy() {

	for (int i=sequence.length()-1; i>=0; i--) {
	    for (int j=i; j<sequence.length(); j++) {

		// Fill in V table entry
		if (j-i <= 3)  // If subsequence has length 4 or less then V entry is infinity
		    V[i][j] = Double.POSITIVE_INFINITY;  // Fill in first four diagonals of V matrix with infinity

		else {  // Compute optimal energy for each subsequence
		    double caseV_1 = e.hairpin(i, j);  // hairpin loop
		    double caseV_2 = e.stacking(i, j, i+1, j-1) + V[i+1][j-1];  // stacking region

		    // interior loop
		    double caseV_3 = Double.POSITIVE_INFINITY;
		    // Add your code here to handle caseV_3, interior loop








		    // bulge loop (right side)
		    double caseV_4 = Double.POSITIVE_INFINITY;
		    for (int y=i+2; y<j-1; y++) {
			double temp = e.bulge(i, j, i+1, y) + V[i+1][y];
			if (temp < caseV_4) caseV_4 = temp;
		    }

		    // bulge loop (left side)
		    double caseV_5 = Double.POSITIVE_INFINITY;
		    for (int x=i+2; x<j-1; x++) {
			double temp = e.bulge(i, j, x, j-1) + V[x][j-1];
			if (temp < caseV_5) caseV_5 = temp;
		    }

		    // bifurcation
		    double caseV_6 = Double.POSITIVE_INFINITY;
		    for (int x=i+2; x<j-2; x++) {
			double temp = W[i+1][x] + W[x+1][j-1];
			if (temp < caseV_6) caseV_6 = temp;
		    }

		    V[i][j] = min(caseV_1, caseV_2, caseV_3, caseV_4, caseV_5, caseV_6);  // Fill in V table entry
		}


		// Fill in W table entry
		if (j-i <= 4)  // If subsequence has length 5 or less then energy is infinity
		    W[i][j] = Double.POSITIVE_INFINITY;  // Fill in first five diagonals of W matrix with infinity

		else {  // Compute optimal energy for each subsequence
		    double caseW_1 = W[i+1][j];  // i doesn't pair
		    double caseW_2 = W[i][j-1];  // j doesn't pair
		    double caseW_3 = V[i][j];    // i and j pair with each other

		    // i and j pair, but not with each other
		    double caseW_4 = Double.POSITIVE_INFINITY;
		    for (int x=i+1; x<j-1; x++) {
			double temp = W[i][x] + W[x+1][j];
			if (temp < caseW_4) caseW_4 = temp;
		    }

		    W[i][j] = min(caseW_1, caseW_2, caseW_3, caseW_4);  // Fill in W table entry
		}
	    }
	}
    }

    // Return the energy of the most favorable (lowest energy) structure
    public double getEnergy() {
	return W[0][sequence.length()-1];
    }

    // Backtrack through dynamic programming tables to create the *structure*
    // with the minimum energy.
    // Return the structure as a String
    public String getStructure() {
        return getStructure(W, 0, sequence.length()-1);
    }




    // Private instance methods

    // Returns minimum of four numbers
    private double min(double a, double b, double c, double d) {
	return Math.min(Math.min(a,b), Math.min(c,d));
    }

    // Returns minimum of six numbers
    private double min(double a, double b, double c, double d, double e, double f) {
	return Math.min(min(a,b,c,d), Math.min(e,f));
    }

    // Return a String of the specified length consisting only of "period" characters
    private String periods(int length) {
        StringBuilder sb = new StringBuilder();
        for (int i=0; i<length; i++)
            sb.append('.');
        return sb.toString();
    }

    // Backtrack through dynamic programming tables to create the *structure*
    // with the minimum energy.
    // Recursively calculate the optimal structure for the subsequence between i and j, inclusive.
    // Return the structure as a String
    private String getStructure(double[][] M, int i, int j) {
        if (M[i][j] > (Double.MAX_VALUE-1.0)) return periods(j-i+1);
	java.text.DecimalFormat df = new java.text.DecimalFormat("0.00");

	// V matrix
        // Determine where current table entry came from (i.e., backtrack)
	if (M == V) {
	    double caseV_1 = e.hairpin(i, j);  // hairpin loop
	    double caseV_2 = e.stacking(i, j, i+1, j-1) + V[i+1][j-1];  // stacking region

	    // interior loop
	    double caseV_3 = Double.POSITIVE_INFINITY;
	    int interiorX = -1;
	    int interiorY = -1;
	    for (int x=i+1; x<j-2; x++) {
		for (int y=x+1; y<j; y++) {
		    double temp = e.interior(i, j, x, y) + V[x][y];
		    if (temp < caseV_3) {
			caseV_3 = temp;
			interiorX = x;
			interiorY = y;
		    }
		}
	    }

	    // bulge loop (right side)
	    double caseV_4 = Double.POSITIVE_INFINITY;
	    int bulgeY = -1;
	    for (int y=i+2; y<j-1; y++) {
		double temp = e.bulge(i, j, i+1, y) + V[i+1][y];
		if (temp < caseV_4) {
		    caseV_4 = temp;
		    bulgeY = y;
		}
	    }

	    // bulge loop (left side)
	    double caseV_5 = Double.POSITIVE_INFINITY;
	    int bulgeX = -1;
	    for (int x=i+2; x<j-1; x++) {
		double temp = e.bulge(i, j, x, j-1) + V[x][j-1];
		if (temp < caseV_5) {
		    caseV_5 = temp;
		    bulgeX = x;
		}
	    }

	    // bifurcation
	    double caseV_6 = Double.POSITIVE_INFINITY;
	    int bifurcationX = -1;
	    for (int x=i+2; x<j-2; x++) {
		double temp = W[i+1][x] + W[x+1][j-1];
		if (temp < caseV_6) {
		    caseV_6 = temp;
		    bifurcationX = x;
		}
	    }

	    if (caseV_1 == V[i][j]) {
		if (energyDetails) System.out.println("Hairpin    " + "\t" + df.format(caseV_1) + "\t\t" + "Closing pair " + sequence.charAt(i) + ":" + sequence.charAt(j));
		return "(" + periods(j-i-1) + ")";
	    }
	    if (caseV_2 == V[i][j]) {
		if (energyDetails) System.out.println("Stacking   " + "\t" + df.format(caseV_2-V[i+1][j-1]) + "\t\t" + "Closing pair " + sequence.charAt(i) + ":" + sequence.charAt(j));
		return "(" + getStructure(V, i+1, j-1) + ")";
	    }
	    if (caseV_3 == V[i][j]) {
		if (energyDetails) System.out.println("Interior   " + "\t" + df.format(caseV_3-V[interiorX][interiorY]) + "\t\t" + "Closing pair " + sequence.charAt(i) + ":" + sequence.charAt(j));
		return "(" + periods(interiorX-i-1) + getStructure(V, interiorX, interiorY) + periods(j-interiorY-1) + ")";
	    }
	    if (caseV_4 == V[i][j]) {
		if (energyDetails) System.out.println("Bulge      " + "\t" + df.format(caseV_4-V[i+1][bulgeY]) + "\t\t" + "Closing pair " + sequence.charAt(i) + ":" + sequence.charAt(j));
		return "(" + getStructure(V, i+1, bulgeY) + periods(j-interiorY-1) + ")";
	    }
	    if (caseV_5 == V[i][j]) {
		if (energyDetails) System.out.println("Bulge      " + "\t" + df.format(caseV_5-V[bulgeX][j-1]) + "\t\t" + "Closing pair " + sequence.charAt(i) + ":" + sequence.charAt(j));
		return "(" + periods(bulgeX-i-1) + getStructure(V, bulgeX, j-1) + ")";
	    }
	    if (caseV_6 == V[i][j]) {
		if (energyDetails) System.out.println("Bifurcation" + "\t" + df.format(W[i+1][bifurcationX]) + "," + df.format(W[bifurcationX+1][j-1]) + "\t\t" + "Closing pair " + sequence.charAt(i) + ":" + sequence.charAt(j));
		return "(" + getStructure(W, i+1, bifurcationX) + getStructure(W, bifurcationX+1, j-1) + ")";
	    }
	}

	// W matrix
        // Determine where current table entry came from (i.e., backtrack)
	if (M == W) {
	    double caseW_1 = W[i+1][j];  // i doesn't pair
	    double caseW_2 = W[i][j-1];  // j doesn't pair
	    double caseW_3 = V[i][j];    // i and j pair with each other

	    // i and j pair, but not with each other
	    double caseW_4 = Double.POSITIVE_INFINITY;
	    int pairX = -1;
	    for (int x=i+1; x<j-1; x++) {
		double temp = W[i][x] + W[x+1][j];
		if (temp < caseW_4) {
		    caseW_4 = temp;
		    pairX = x;
		}
	    }

	    if (caseW_1 == W[i][j]) return "." + getStructure(W, i+1, j);
	    if (caseW_2 == W[i][j]) return getStructure(W, i, j-1) + ".";
	    if (caseW_3 == W[i][j]) return getStructure(V, i, j);
	    if (caseW_4 == W[i][j]) return getStructure(W, i, pairX) + getStructure(W, pairX+1, j);
	}

	return "";
    }

    // Returns a String represenation of a 2D double array
    private String arrayToString(double[][] a) {
	if (a.length == 0) return "";
	if (a[0].length == 0) return "";
	if (a.length != sequence.length()) {
	    System.err.println("Error - when converting 2D array to String, the array length must be the same as the length of the RNA sequence.");
	}

	StringBuilder sb = new StringBuilder();
	java.text.DecimalFormat df = new java.text.DecimalFormat("0.0");
	for (int i=sequence.length()-1; i>=0; i--) {
	    sb.append(sequence.charAt(i));
	    for (int j=0; j<sequence.length(); j++) {
		if (j < i) sb.append("\t" + " X ");
		else if (a[i][j] >= (Double.MAX_VALUE-1.0)) sb.append("\t" + " * ");  // Positive Infinity
		else sb.append("\t" + df.format(a[i][j]));
	    }
	    sb.append("\n");
	}
	for (int i=0; i<sequence.length(); i++)
	    sb.append("\t" + sequence.charAt(i));
	sb.append("\n");
	return sb.toString();
    }



    // Main method
    public static void main(String[] args) {

        if (args.length < 1) {
            System.err.println("\nTo execute this program, one command line argument is required corresponding to the name of a FASTA file containing an RNA sequence. For example,\n");
            System.err.println("\tjava RNA_Energy sequences/MiR-96.txt\n");
	    return;
        }

	RNA_Energy r = new RNA_Energy(args[0]);
	r.computeStructureWithMinimumEnergy();
	//System.out.println(r.arrayToString(r.V));
	//System.out.println(r.arrayToString(r.W));
	r.energyDetails = false;
	System.out.println("\nStructure with lowest energy:\n");
	System.out.println("\t" + r.getStructure());
	System.out.println("\t" + r.sequence + "\n");
	System.out.println("Energy of most favorable structure:\t" + r.getEnergy() + "\n");
    }

}
