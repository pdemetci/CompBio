/*
 * @author : Pinar Demetci
 * Computational Biology - Project 6
 */
import java.util.*;

/*
 * EM_MotifSearch class inherits from MotifSearch
 * Fields inherited from MotifSearch class are:
 * instanceLocations:
 * matrix:
 * motifLength:
 * numSequences:
 * sequences:
 */

public class EM_MotifSearch extends MotifSearch {

	 public EM_MotifSearch (String fileName, int motifLength) {
	        super(fileName, motifLength); 
	    }
	 
	   public static void main(String args[]) {
	        EM_MotifSearch EM = new EM_MotifSearch(args[0], Integer.parseInt(args[1]));
	        EM.run_EM_multiple_times(Integer.parseInt(args[2]));
	        System.out.println("Motif Matrix:");
	        System.out.println(EM.matrixToString());
	        System.out.println("Consensus sequence: " + EM.getConsensusSequence() + "\n");
	        System.out.println("Information content: " + EM.getInformationContentOfMatrix() + "\n");
	        System.out.println("Motifs in Strings:");
	        System.out.println("\n" + EM.motifInstancesToString());
	    }
	    public void EM() {
	        setRandomLocationsForMotifInstances();
	        double initialInfoC;
	        double recentInfoC;
	        double storeInfoC =-1.0;
	        do{
	        	initialInfoC = storeInfoC;
	        	determineMatrixModel();
	            addPseudocountsToMatrix();
	            recentInfoC = getInformationContentOfMatrix();
	            storeInfoC = recentInfoC;
	            determineMotifInstances();
	        }while(recentInfoC!=initialInfoC);     
	    }
	    
	    public void run_EM_multiple_times(int iterations) {
	    	int i=0;
	    	double InfoContent=-1;
	    	Vector<Integer> candidateMotifIndeces = new Vector <Integer>();
	    	double[][] MotifMatrix = new double[4][motifLength];
	    	do{
	    		EM();
	    		double recentInfoContent = getInformationContentOfMatrix();
	    		if (recentInfoContent>InfoContent){
	    			InfoContent = recentInfoContent;
	    			candidateMotifIndeces = getInstanceLocations();
	    			MotifMatrix = getMatrix();
	    		}
	    		i++;	
	    	}while(i<iterations);
	    }
	 
	 public void setRandomLocationsForMotifInstances(){
	        Random r = new Random();
	        for(int i = 0; i < numSequences; i++) {
	        	String sequence = sequences.get(i);
	            int random = r.nextInt((sequence.length() - motifLength) + 1);
	            instanceLocations.set(i, Integer.valueOf(random));   
	        }
	        
	    }
	 
	 
	 public void determineMatrixModel(){
		 List<Character> nucleotides = Arrays.asList('A', 'C', 'G', 'T');
		 for (int i=0; i<numSequences; i++){
			 String sequence = sequences.get(i);
			 for (int j=0; j<motifLength; j++){
				 char nucleotide = sequence.charAt(instanceLocations.get(i)+j);
				 matrix[nucleotides.indexOf(nucleotide)][j]=(matrix[nucleotides.indexOf(nucleotide)][j]+1.0)/ (numSequences*1.0);
			 }
		 }
	 }
	 
	 public void determineMotifInstances(){
		for (int i=0; i<numSequences; i++){
			Integer index = getMotifForSequence(sequences.get(i));
			if (index == -1){
				System.out.println("There seems to be a problem with getMotifForsequence fuction");
			}else{
				instanceLocations.set(i, index);
			}
		} 
	 }
	 /*
	  * Helper Function for determineMotifInstances
	  */
	 public Integer getMotifForSequence(String s){
		 double maxScore = -1.0;
		 Integer index = -1;
		 for (int i=0; i<(s.length()- motifLength +1); i++){
			 String candidate = s.substring(i, i+motifLength);
			 double score =  motifScore(candidate);
		 if (score > maxScore){
			 maxScore = score;
			 index = i;
		 }
	 } 
		 return index;
}
	 /*
	  * Helper Function for getMotifForSequence
	  */
	 public Double motifScore(String s){
		 double score =1.0;
		 List<Character> nucleotides = Arrays.asList('A', 'C', 'G', 'T');
	     for(int i = 0; i < motifLength;i++){
	    	 char c = s.charAt(i);
	    	 score *= matrix[nucleotides.indexOf(c)][i]; 
	     }
	     return score;
	 }
	 
	
	 public Double getInformationContentOfMatrix(){
		double informationContent = 0.0;
		char[] nucleotide= {'A','C','G','T'};
		for (int i=0; i<motifLength; i++){
			for (int j=0; j<4; j++){
				double score = matrix[j][i]*(Math.log10(matrix[j][i]/ getNucleotideContent(nucleotide[j])) / Math.log10(2.0)) ;
				informationContent+= score;
			}
		} 
		 return informationContent;
	 }
}