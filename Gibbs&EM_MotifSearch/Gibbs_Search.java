/*
 * @author : Pinar Demetci
 * Computational Biology - Project 6
 */

import java.util.*;

public class Gibbs_Search extends EM_MotifSearch{
	
	public Gibbs_Search(String fileName, int motifLength) {
        super(fileName, motifLength);
    }
	
	public static void main(String args[]){
		Gibbs_Search Gibbs = new Gibbs_Search(args[0], Integer.parseInt(args[1]));
		Gibbs.run_EM_multiple_times(Integer.parseInt(args[2]));
		System.out.println("Motif Matrix:");
		System.out.println(Gibbs.matrixToString());
		System.out.println("Consensus sequence: " + Gibbs.getConsensusSequence());
		System.out.println("Information content: " + Gibbs.getInformationContentOfMatrix());
		System.out.println("Motifs in Strings:");
		System.out.println("\n" + Gibbs.motifInstancesToString());
	}
	
	public void determineMotifInstances() {	
		
		for(int i = 0; i < numSequences; i++) {
			Vector<Double> scores = motifScoresForSequence(sequences.get(i));
			int indexOfSampledMotifInstance = getIndexViaSampling(scores);
			instanceLocations.set(i, indexOfSampledMotifInstance);
	}
	}
	
	/*
	 * Helper Function for determineMotifInstances
	 */
	public Vector<Double> motifScoresForSequence(String sequence){
		Vector<Double> scores = new Vector<Double>();
		for (int i=0; i<(sequence.length()- motifLength +1); i++){
			String candidate = sequence.substring(i, i+motifLength);
			double score =  motifScore(candidate);
			scores.add(score); 
		}
		return scores;
	}
	/*
	 * Returns the index of a randomly sampled value in a Vector.
	One value from the Vector is chosen at random and the value's index (not the value itself) is returned. The value is not chosen uniformly at random, but rather via sampling, i.e., higher values are more likely to be chosen and lower values are less likely to be chosen.
	One approach for randomly sampling a collection of values proceeds as follows:
    1. Normalize the values in the collection so that they sum to 1.0. The values now represent a probability distribution.
    2. Convert the values from a probability distribution to a cumulative distribution. In a cumulative distribution, the value at index i represents the sum of all values at indices less than or equal to i in the probability distribution. The final value in a cumulative distribution should be 1.0 since the sum of all values in a probability distribution is 1.0.
    3. Generate a number uniformly at random between 0.0 and 1.0. Return the index of the smallest value in the cumulative distribution that is at least as big as the random number.
	 */
	public int getIndexViaSampling(Vector<Double> scores) {
		//there might be a better way than having countless loops. 
		//but I couldn't figure out one.
		
		//Normalizing the values: 
		double sum=0;
		for (int i=0; i<scores.size(); i++){
			sum+= scores.get(i);
		}
		for (int i=0; i<scores.size();i++){
			scores.set(i, scores.get(i)/sum);
		}
		
		//turning this into a cumulative distribution:
		for (int i=1; i<scores.size(); i++){
			scores.set(i, (scores.get(i)+scores.get(i-1)));
		}
		
		//generate a random value and return score value that is at least as small as the random number:
		int index = scores.size()-1;

		Random r = new Random();
		double rand = r.nextDouble();
		
		for (double c: scores){
			if(rand <= c){
				index= scores.indexOf(c);
			}
		}
		return index;
	}
	
}