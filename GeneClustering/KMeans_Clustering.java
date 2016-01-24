/*
 * author: Pinar Demetci
 */

import java.util.*;

public class KMeans_Clustering extends Clustering{

	int numberClusters; 
    public KMeans_Clustering(String filename, int numClusters){
        super(filename);
        numberClusters=numClusters;
    }
    /*
     * Initializes each cluster so that each cluster contains zero genes
     */

    public void initializeAllClusters(){
    	for(int i = 0; i < numberClusters; i++){
    		Cluster cluster = new Cluster();
    		clusters.add(cluster);
    		cluster.initialize();
    	}
    }
/*
 * Return a collection of the mean (average) expression vectors for all of the clusters.
 */
    public Vector<Vector<Double>> getMeansOfAllClusters()
    {
        Vector<Vector<Double>> means = new Vector<Vector<Double>>();
        for(int i = 0; i < getNumClusters(); i++){
        	Cluster cluster = clusters.get(i);
            means.add(cluster.getClusterMean());
        }
        return means;
    }
/*
 * If any clusters are empty (contain zero genes), then genes are moved from clusters containing multiple genes.
 */
    public void populateEmptyClusters(){
    	List<Cluster> emptyClusters = new ArrayList<Cluster>();
    	List<Cluster> multipleClusters = new ArrayList<Cluster>();
        for(int i = 0; i < getNumClusters(); i++){
        	Cluster cluster = clusters.get(i);
        	if (cluster.getSizeOfCluster() == 0){
        		emptyClusters.add(cluster);
        	}else if (cluster.getSizeOfCluster()>1){
        		multipleClusters.add(cluster);
        	}
        }
        
      	for(int e=0; e<emptyClusters.size(); e++){
      		int r;
      		do{Random random = new Random();
        	r = random.nextInt(multipleClusters.size());
      		}while(multipleClusters.get(r).getSizeOfCluster()<2);
      		
        	Cluster empty = emptyClusters.get(e);
	        Cluster multiple = multipleClusters.get(r);
	        Gene gene = multiple.getGene(0);
	        empty.addGene(gene);
	        multiple.removeGene(gene);
        	}
        	}
    /*
     * Assigns each gene to a random cluster.
     */
    public void randomlyAssignGenesToClusters(){
        initializeAllClusters();
        for(int i = 0; i < getNumGenes(); i++){
        	Gene gene = genes.get(i);
        	Random random = new Random();
            int c = random.nextInt(getNumClusters());
            Cluster cluster = clusters.get(c);
            (cluster).addGene(gene);
        }

        populateEmptyClusters();
    }
    /*
     * Assigns each gene to the cluster whose mean expression vector is closest to the gene.
     */
    public boolean assignGenesToClusters(Vector vector){
    	List<Object> oldClusters = new ArrayList<Object>();
    	List<Object> newClusters = new ArrayList<Object>();
 
    	for(int i=0; i<getNumGenes(); i++){
    		Gene gene = genes.get(i);
    		int clusterNum=-100;
    		Vector<Double> distances = new Vector<Double>();
    		for (int j=0; j<vector.size(); j++){
    			Vector<Double> expression = new Vector<Double>();
    			if (vector.get(j) instanceof Vector<?>){
    				expression = (Vector<Double>)vector.get(j); 
    				//I've been getting compilation errors saying 'unsafe cast from Object to Vector<Double>
    				//but I'm pretty sure getMeansofAllClusters returns Vector<Vector<Double>>. 
    				//I could not figure out a way around it so I just finished writing the code.
    			}
    			double distance = gene.distanceToExpressionVector(expression);
    			distances.add(distance);
    		}
    		double minDistance = Collections.min(distances);
    		vector.indexOf(minDistance);
    		Cluster newCluster = new Cluster();
    		newCluster.addGene(gene);
    		newClusters.add(newCluster);
    		for (int j=0; j<getNumClusters(); j++){
    			if (clusters.get(j).isGeneInCluster(gene)==true){
    				clusterNum = j;
    			}
    		}
    		Cluster oldCluster = clusters.get(clusterNum);
    		oldCluster.removeGene(gene);
    		oldClusters.add(oldCluster);
    		getMeansOfAllClusters();
    	}
    	
    	if (oldClusters.equals(newClusters)){
    		return true;
    	}else{
    		return false;
    	}
    	
    }
    	/*
    	 * Performs k-means clustering.
    	 */
    public void kMeans(){
    	initializeAllClusters();
    	randomlyAssignGenesToClusters();
    	do{
    		if (getMeansOfAllClusters() instanceof Vector<?>){
    		assignGenesToClusters(getMeansOfAllClusters());
    	} 
    	}while (assignGenesToClusters(getMeansOfAllClusters())!=true);
 
    }

    public static void main(String args[]){
        KMeans_Clustering kmeans_clustering = new KMeans_Clustering(args[0], Integer.parseInt(args[1]));
        kmeans_clustering.kMeans();
        System.out.println(kmeans_clustering.toString());
    }
}