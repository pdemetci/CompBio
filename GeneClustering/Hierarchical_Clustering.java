/*
 * author: Pinar Demetci
 */
public class Hierarchical_Clustering extends Clustering{
	//Hierarchical_Clustering inherits from Clustering class
	public Cluster closestCluster1 = new Cluster();
	public Cluster closestCluster2 = new Cluster();
	int numberClusters;
	
	public Hierarchical_Clustering(String filename, int numClusters){
		super (filename);
		numberClusters = numClusters;
	}
	/*
	 * Creates a Clustering based on gene and experiment data from a tab-delimited text file.
	 */
	public static void main(String[] args){
		Hierarchical_Clustering hierarchical_clustering = new Hierarchical_Clustering(args[0], Integer.parseInt(args[1]));
		hierarchical_clustering .hierarchical();
		System.out.println(hierarchical_clustering .toString());
		
	}
	
	/*
	 * Performs centeroid-linkage hierarchical clustering.
	 */
	public void hierarchical(){
		initiallyAssignOneGeneToEachCluster();
		while (getNumClusters()>numberClusters){
			mergeTwoClosestClusters();
		}
	}
	
	/*
	 * Assigns each gene to its own unique cluster.
	 */
	public void initiallyAssignOneGeneToEachCluster(){
		for (int i=0; i<getNumGenes(); i++){
			Cluster cluster = new Cluster();
			cluster.addGene(genes.get(i));
			clusters.add(cluster);
		}
	}
	
	/*
	 * Identifies and merges together the two closest clusters.
	 */
	public void mergeTwoClosestClusters(){
		findClosestClusters();
		closestCluster1.absorbCluster(closestCluster2);
		clusters.remove(closestCluster2);
	}
	/*
	 * Helper Function for mergeTwoClosestClusters
	 * Initially I was storing all distances 
	 * When I searched for it, found this method from:
	 * 
	 */
	private void findClosestClusters(){
		double minDistance = clusters.get(1).getDistanceToCluster(clusters.get(2)); //Randomly initializing the minDistance.
		//will be replaced when we find smaller distances in the following loop:
		for (int i=0; i<getNumClusters()-1; i++){
			for (int j = i+1; j<getNumClusters(); j++){
				Cluster cluster1 = (Cluster)clusters.get(i);
				Cluster cluster2 = (Cluster)clusters.get(j);
				double distance = cluster1.getDistanceToCluster(cluster2);
				if (distance<minDistance){
					minDistance=distance;
					closestCluster1=cluster1;
					closestCluster2=cluster2;
				}
			}
		}
		
}
}
