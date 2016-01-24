import java.util.*;  // Needed for Scanner class and for Vector class
import java.io.*;  // Needed for File class

/******************************************************************
 * An instance of the <code>Clustering</code> class represents a
 * clustering (i.e., grouping or partitioning) of a collection of genes.
 ******************************************************************/
public class Clustering {

    /**************************************************************
     ********************** INSTANCE VARIABLES ********************
     **************************************************************/

    /**
     * A collection of experiment names
     */
    protected Vector<String> experiments;

    /**
     * A collection of genes
     */
    protected Vector<Gene> genes;

    /**
     * A collection of clusters
     */
    protected Vector<Cluster> clusters;



    /**************************************************************
     ********************** CONSTRUCTOR ***************************
     **************************************************************/

    /** 
     * Creates an initially empty <code>Clustering</code>.
     * <p>
     * A set of genes and experiments are determined from the specified <code>String</code> representing the name of a file.
     * Genes and experiments are read-in from the tab-delimited file. Initially, the constructed
     * <code>Clustering</code> is empty.
     *
     * @param   fileName   the name of a tab-delimited text file containing gene and experiment data
     */
    public Clustering(String fileName) {
	experiments = getExperimentNamesFromFile(fileName);
	genes = getGeneInformationFromFile(fileName);
	System.err.println("\nRead in " + getNumExperiments() + " experiments and " + getNumGenes() + " genes.");
	clusters = new Vector<Cluster>();
    }



    /**************************************************************
     ********************** PUBLIC INSTANCE METHODS ***************
     **************************************************************/

    /**
     * Returns the number of experiments in this <code>Clustering</code>.
     *
     * @return   an integer representing the number of experiments
     */
    public int getNumExperiments() {
	return experiments.size();
    }

    /**
     * Returns the total number of genes in this <code>Clustering</code>.
     *
     * @return   an integer representing the number of genes
     */
    public int getNumGenes() {
	return genes.size();
    }

    /**
     * Returns the number of clusters in this <code>Clustering</code>.
     *
     * @return   an integer representing the number of clusters
     */
    public int getNumClusters() {
	return clusters.size();
    }

    /**
     * Returns a collection of experiment names for which there is data in the specified file.
     * <p>
     * The method assumes that the first line of the specified file
     * contains a tab-delimited list of experiment names. The first
     * two tokens on the first line do not correspond to experiment
     * names and are ignored. Thus, if the first line of the file 
     * contains "m+2" tab-delimited tokens, then the line contains
     * "m" experiment names. A <code>Vector</code> containing the
     * names of the experiments is returned.
     *
     * @param   fileName   the name of a tab-delimited text file with gene and experiment data
     * @return   a collection of experiment names as extracted from the specified file
     */
    public Vector<String> getExperimentNamesFromFile(String fileName) {
	Vector<String> experimentNames = new Vector<String>();
	try {
	    Scanner r = new Scanner(new File(fileName));
	    String firstLine = r.nextLine();  // Read first line of file
	    String[] firstLineSplit = firstLine.split("\t", -1);  // Split line based on tabs
	    for (int j=2; j<firstLineSplit.length; j++)  // Ignore first two tokens
		experimentNames.add(firstLineSplit[j]);
	    r.close();
	} catch (FileNotFoundException e) {
	    System.err.println("Error - could not read in file " + fileName);
	}
	return experimentNames;
    }

    /**
     * Returns a collection of <code>Genes</code> based on data in the specified file.
     * <p>
     * The method assumes that the specified file is tab-delimited.
     * The method assumes that the first line of the specified file
     * is a header line that does not contain any gene information.
     * Starting with the second line of the file and continuing until
     * the end of the file, the method assumes that each line is a
     * tab-delimited collection of gene information, including
     * the gene's name, function, and expression values.
     * A <code>Vector</code> of <code>Gene</code> objects is returned.
     *
     * @param   fileName   the name of a tab-delimited text file with gene and experiment data
     * @return   a collection of <code>Genes</code> as extracted from the specified file
     */
    public Vector<Gene> getGeneInformationFromFile(String fileName) {
	genes = new Vector<Gene>();
	try {
	    Scanner r = new Scanner(new File(fileName));
	    String firstLine = r.nextLine();  // Ignore first (header) line of file
	    while (r.hasNextLine()) {  // Continue until we reach end of file
		String line = r.nextLine();  // Read in next line of file
		genes.add(new Gene(line));
	    }
	    r.close();
	} catch (FileNotFoundException e) {
	    System.err.println("Error - could not read in file " + fileName);
	}
	return genes;
    }

    /**
     * Returns a <code>String</code> representation of all <code>Clusters</code> in this <code>Clustering</code>.
     *
     * @return   a <code>String</code> represenation of this <code>Clustering</code>
     */
    public String toString() {
	StringBuilder sb = new StringBuilder();
	for (int k=0; k<getNumClusters(); k++) {
	    sb.append("Cluster # " + k + " containing " + clusters.get(k).getSizeOfCluster() + " genes.\n");
	    sb.append(clusters.get(k).toString());
	}
	return sb.toString();
    }



    /*********************************************************
     ****************** MAIN METHOD **************************
     *********************************************************/

    /**
     * The <code>main</code> method creates an initially empty <code>Clustering</code> from a tab-delimited text file of gene and experiment data.
     * 
     * @param   args   an array of <code>Strings</code> representing any command line arguments
     */
    public static void main(String [] args) {

        if (args.length < 1) {
            System.out.println("\nTo execute this program, a command line argument corresponding\nto a file name is required. For example,\n");
            System.out.println("\tjava Clustering yeast_10.txt\n\n");
            System.exit(0);
        }

	String fileName = args[0];
	Clustering c = new Clustering(fileName);
	
    }

}