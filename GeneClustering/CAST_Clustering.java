/*
 * author: Pinar Demetci
 */
public class CAST_Clustering extends Clustering{
	  public double affinityThreshold;
	  public Cluster initialCluster = new Cluster();
	
    public CAST_Clustering(String filename, double d){
        super(filename);
        affinityThreshold = d;
    }
    
    public static void main(String args[]){
        double affinityThreshold= Double.parseDouble(args[1]);
        CAST_Clustering cast_clustering = new CAST_Clustering(args[0], affinityThreshold);
        cast_clustering.cast();
        System.out.println(cast_clustering.toString());
    }
    /*
     * Performs CAST clustering.
     */
   public void cast(){
    	
    	for(int i = 0; i < getNumGenes(); i++){
            initialCluster.addGene(genes.get(i));
    	}
        
        while (initialCluster.getSizeOfCluster()>0){
        	Cluster cluster = new Cluster();
        	addGenesWithHighAffinity(cluster);
        	removeGenesWithLowAffinity(cluster);
        	clusters.add(cluster);
        }
    }
/*
 * Adds to the specified Cluster all unassigned genes that are closer to the genes in the Cluster, on average, than the affinity threshold.
 */
    public int addGenesWithHighAffinity(Cluster cluster){
        for(int i = 0; i < initialCluster.getSizeOfCluster(); i++){
            Gene gene = initialCluster.getGene(i);
            double distance = gene.distanceToExpressionVector(cluster.getClusterMean());
            if(distance < affinityThreshold && cluster.isGeneInCluster(gene)==false){
            	cluster.addGene(gene);
            	initialCluster.removeGene(gene);
            }
        }
        return cluster.getSizeOfCluster();
    }
/*
 * Removes from the specified Cluster any genes that are farther from the genes in the Cluster, on average, than the affinity threshold.
 */
    public int removeGenesWithLowAffinity(Cluster cluster){
        for(int i = 0; i < cluster.getSizeOfCluster(); i++){
        	Gene gene = cluster.getGene(i);
        	double distance = cluster.getGene(i).distanceToExpressionVector(cluster.getClusterMean());
            if( distance > affinityThreshold){
            	cluster.removeGene(gene);
            	initialCluster.addGene(gene);
        }
        }
        return cluster.getSizeOfCluster();
    }
} 