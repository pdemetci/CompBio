import java.util.Vector;
// Vector class

/******************************************************************
 * An instance of the <code>Cluster</code> class represents a
 * group of genes.
 ******************************************************************/
public class Cluster {

    /**************************************************************
     ********************** INSTANCE VARIABLES ********************
     **************************************************************/

    private Vector<Gene> genesInCluster;
    private Vector<Double> mean;  // Mean expression vector for cluster
    private boolean clusterMeanIsAccurate;  // Indicates if the "mean" is up-to-date. The mean
                                            // is not up-to-date if genes have been added/removed 
                                            // from the cluster since the last time the "mean" 
                                            // was computed.



    /**************************************************************
     ********************** CONSTRUCTOR ***************************
     **************************************************************/

    /** 
     * Creates a <code>Cluster</code> initially containing zero genes.
     */
    public Cluster() {
	genesInCluster = new Vector<Gene>();
	mean = new Vector<Double>();
	clusterMeanIsAccurate = false;
    }



    /**************************************************************
     ********************** PUBLIC INSTANCE METHODS ***************
     **************************************************************/

    /**
     * Clears all genes from this <code>Cluster</code> so that it contains zero genes.
     */
    public void initialize() {
	genesInCluster.clear();
	mean.clear();
	clusterMeanIsAccurate = false;
    }

    /**
     * Returns the number of genes in this <code>Cluster</code>.
     *
     * @return an integer representing the size of this <code>Cluster</code>
     */
    public int getSizeOfCluster() {
	return genesInCluster.size();
    }

    /**
     * Adds the specified <code>Gene</code> to this <code>Cluster</code>.
     * <p>
     * If the specified <code>Gene</code> is already in this
     * <code>Cluster</code> then the <code>Gene</code> is not
     * added again (i.e,. there should be no duplicate
     * <code>Genes</code> in the <code>Cluster</code>).
     *
     * @param   g   a <code>Gene</code> to be added to this <code>Cluster</code>
     */
    public void addGene(Gene g) {
	if (!genesInCluster.contains(g))
	    genesInCluster.add(g);
	clusterMeanIsAccurate = false;
    }

    /**
     * Removes the specified <code>Gene</code> from this <code>Cluster</code>.
     * <p>
     * If the specified <code>Gene</code> is not in this <code>Cluster</code>
     * then this <code>Cluster</code> is not changed.
     *
     * @param   g   a <code>Gene</code> to be removed from this <code>Cluster</code>
     */
    public void removeGene(Gene g) {
	genesInCluster.remove(g);
	clusterMeanIsAccurate = false;
    }

    /**
     * Returns a <code>Gene</code> in this <code>Cluster</code>.
     * <p>
     * The ith <code>Gene</code> in this <code>Cluster</code> is returned.
     *
     * @param   i   the index of the <code>Gene</code> in this <code>Cluster</code> to be returned
     * @return   the <code>Gene</code> at index <code>i</code> in this <code>Cluster</code>
     */
    public Gene getGene(int i) {
	return genesInCluster.get(i);
    }

    /**
     * Returns a <code>boolean</code> value indicating if the specified
     * <code>Gene</code> is in this <code>Cluster</code>.
     *
     * @param   g   a <code>Gene</code> that may or may not be in this <code>Cluster</code>
     * @return      <code>true</code> if <code>g</code> is in this <code>Cluster</code>, false otherwise
     */
    public boolean isGeneInCluster(Gene g) {
	for (int i=0; i<getSizeOfCluster(); i++) {
	    if (getGene(i) == g) return true;
	}
	return false;
    }

    /**
     * Returns the mean (i.e., average) of the expression vectors of all <code>Genes</code> in this <code>Cluster</code>.
     *
     * @return   a collection of <code>doubles</code> corresponding to the average expression of each <code>Gene</code> in this <code>Cluster</code> in each experiment
     */
    public Vector<Double> getClusterMean() {
	if (getSizeOfCluster() == 0) return new Vector<Double>();

	// The length of the any gene's expression vector corresponds to the number of experiments
	int numExperiments = genesInCluster.get(0).getExpressionVector().size();

	// Current mean is accurate, so create copy of current mean and return it
	if (clusterMeanIsAccurate) {
	    Vector<Double> newMean = new Vector<Double>();
	    for (int j=0; j<numExperiments; j++)
		newMean.add(mean.get(j));
	    return newMean;
	}

	// Initialize vector of expression values. The vector should contain
	// one expression value for each experiment.
	mean.clear();
	for (int j=0; j<numExperiments; j++) {
	    mean.add(0.0);
	}

	// For all genes in the cluster, sum the genes' expression values in each experiment
	for (int i=0; i<getSizeOfCluster(); i++) {
	    for (int j=0; j<numExperiments; j++) {
		mean.set(j, mean.get(j) + genesInCluster.get(i).getExpressionInOneExperiment(j));
	    }
	}

	// Divide the sum by the number of genes to obtain the mean (average)
	for (int j=0; j<numExperiments; j++)
	    mean.set(j, mean.get(j) / getSizeOfCluster());

	// Current mean is now accurate
	clusterMeanIsAccurate = true;

	// Create a copy of the curren mean and return it
	Vector<Double> newMean = new Vector<Double>();
	for (int j=0; j<numExperiments; j++)
	    newMean.add(mean.get(j));
	return newMean;
    }

    /**
     * Computes and returns the distance between this <code>Cluster</code> and the specified <code>Cluster</code>.
     * <p>
     * Each <code>Cluster</code> has a vector representing the <code>Cluster</code>'s mean expression
     * across all experiments. The distance between two <code>Clusters</code> is computed as
     * the Euclidean distance between the two <code>Clusters</code>' means.
     *
     * @param   c   a <code>Cluster</code> for which we wish to compute to the distance to this <code>Cluster</code>
     * @return   a <code>double</code> representing the distance between the specified <code>Cluster</code> and this <code>Cluster</code>
     */
    public double getDistanceToCluster(Cluster c) {
	Vector<Double> mean1 = this.getClusterMean();
	Vector<Double> mean2 = c.getClusterMean();

	// Ensure the two cluster means are the same length
	if (mean1.size() != mean2.size()) {
	    System.err.println("Error - the means of two clusters have differing lengths.");
	    return Double.MAX_VALUE;
	}

	// Compute the Euclidean distance between the two vectors (i.e., between the two cluster means)
	double distance = 0.0;
	for (int j=0; j<mean1.size(); j++) {
	    distance += (mean1.get(j) - mean2.get(j)) * (mean1.get(j) - mean2.get(j));
	}
	distance = Math.sqrt(distance);
	return distance;
    }

    /**
     * Adds all genes in the specified <code>Cluster</code> to this <code>Cluster</code>.
     * <p>
     * Thus, this <code>Cluster</code> grows in size by the number of genes
     * in <code>Cluster c</code>. The specified <code>Cluster</code> is not altered.
     *
     * @param   c   a <code>Cluster</code> to be absorbed by this <code>Cluster</code>
     */
    public void absorbCluster(Cluster c) {
	for (int i=0; i<c.getSizeOfCluster(); i++) {
	    this.addGene(c.genesInCluster.get(i));
	}
    }

    /**
     * Returns a <code>String</code> representation of this <code>Cluster</code>.
     * <p>
     * The name of each <code>Gene</code> in this <code>Cluster</code> along with each
     * <code>Gene</code>'s function are included in the <code>String</code> 
     * representation of this <code>Cluster</code>.
     *
     * @return   a <code>String</code> represenation of this <code>Cluster</code>
     */
   public String toString() {
	StringBuilder sb = new StringBuilder();
	for (int i=0; i<getSizeOfCluster(); i++) {
	    sb.append(genesInCluster.get(i).toString());
	}
	sb.append("\n");
	return sb.toString();
    }

}