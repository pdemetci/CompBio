import java.util.Vector;
// Needed for Vector class

/******************************************************************
 * An instance of the <code>Gene<code> class represents (1) the
 * name of a gene (2) the reported function of the gene and (3) the
 * gene's expression values from a set of experiments.
 ******************************************************************/
public class Gene {

    /**************************************************************
     ********************** INSTANCE VARIABLES ********************
     **************************************************************/

    private String name;
    private String function;
    private Vector<Double> expressionValues;



    /**************************************************************
     ********************** CONSTRUCTOR ***************************
     **************************************************************/

    /** 
     * Creates a <code>Gene</code> from a tab-delimited String.
     *
     * @param   lineFromFile   a tab-delimited <code>String</code>, possibly read-in from a file
     */
    public Gene(String lineFromFile) {
	fromString(lineFromFile);
    }



    /**************************************************************
     ********************** PUBLIC INSTANCE METHODS ***************
     **************************************************************/

    /**
     * Returns a <code>String</code> representation of this <code>Gene</code>, including the gene's name and function.
     *
     * @return   a <code>String</code> represenation of this <code>Gene</code>
     */
    public String toString() {
	return name + "\t" + function + "\n";
    }

    /**
     * Returns a <code>String</code> representation of the gene's name.
     *
     * @return   a <code>String</code> represenation of this <code>Gene</code>'s name
     */
    public String getGeneName() {
	return name;
    }

    /**
     * Returns a <code>String</code> representation of the gene's function.
     *
     * @return   a <code>String</code> represenation of this <code>Gene</code>'s function
     */
    public String getGeneFunction() {
	return function;
    }

    /**
     * Returns the expression value of this <code>Gene</code> in the specified experiment.
     *
     * @param   experiment   the index of a particular experiment
     * @return   a <code>double</code> representing this <code>Gene</code>'s expression value
     *
     */
    public double getExpressionInOneExperiment(int experiment) {
	if ((experiment < 0) || (experiment >= expressionValues.size())) {
	    System.err.println("Error - experiment index " + experiment + " is invalid.");
	    return 0.0;
	}
	return expressionValues.get(experiment);
    }

    /**
     * Returns a copy of this <code>Gene</code>'s expression values in all experiments.
     *
     * @return   a <code>Vector</code> of <code>doubles</code> corresponding to this <code>Gene</code>'s expression values
     */
    public Vector<Double> getExpressionVector() {
	Vector<Double> temp = new Vector<Double>();
	for (int j=0; j<expressionValues.size(); j++) {
	    double value = expressionValues.get(j);
	    temp.add(value);
	}
	return temp;
    }

    /**
     * Returns the Euclidean distance between this <code>Gene</code>'s expression values and the specified collection of expression values.
     *
     * @param   vector   a collection of expression values
     * @return   a <code>double</code> corresponding to the distance between two sets of expression values
     */
    public double distanceToExpressionVector(Vector<Double> vector) {
	if (expressionValues.size() != vector.size()) {
	    System.err.println("Error - cannot compute the distance between two vectors as vectors differ in length.");
	    System.out.println(toString());
	    System.out.println("\n" + vector);
	    System.exit(0);
	    return Double.MAX_VALUE;
	}

	double sum = 0.0;
	for (int j=0; j<expressionValues.size(); j++) {
	    sum += (expressionValues.get(j) - vector.get(j)) * (expressionValues.get(j) - vector.get(j));
	}
	return Math.sqrt(sum);
    }



    /**************************************************************
     ********************** PRIVATE INSTANCE METHODS **************
     **************************************************************/

    /**
     * Initializes the state of this gene based the specified tab-delimited <code>String</code>.
     * <p>
     * The first token in the <code>String</code> should correspond
     * to the gene's name. The second token in the <code>String</code> should correspond to the 
     * gene's function. The remaining token in the <code>String</code> represent the
     * expression values of the gene in the corresponding experiments. Thus, if there 
     * are "m" experiments assaying the gene's level of expression, the <code>String</code>
     * should contain "m+2" tokens, with the final "m" tokens corresponding to 
     * decimal numbers representing the gene's exprssion values in the "m" experiments.
     *
     * @param   lineFromFile   a tab-delimited <code>String</code>, possibly read-in from a file
     */
    private void fromString(String lineFromFile) {
	String[] lineSplit = lineFromFile.split("\t", -1);
	if (lineSplit.length < 3) {
	    System.err.println("Error - line does not contain gene name, gene function, and gene expression data.");
	    return;
	}
	name = lineSplit[0];
	function = lineSplit[1];
	expressionValues = new Vector<Double>();
	for (int j=2; j<lineSplit.length; j++) {
	    // Attempt to convert each value to a decimal (i.e., double)
	    double value = 0.0;
	    try {
		value = Double.parseDouble(lineSplit[j]);
	    } catch (NumberFormatException e) {
		//System.err.println("Error - could not convert " + lineSplit[j] + " to a *double*.");
	    }
	    expressionValues.add(value);
	}
    }



}