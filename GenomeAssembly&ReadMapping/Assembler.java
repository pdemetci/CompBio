import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;

/****************************************************************
 * The <code>Assembler</code> application assembles a genome
 * from a set of short sequencing reads.
 ****************************************************************/

public class Assembler {

    /******************************************************
     ****************** CLASS VARIABLES *******************
     ******************************************************/

  /**
   * Length of k-mers
   */
  private static final int k = 25;
  
  /**
   * Nucleotide characters
   */
  private static final char[] NTs = {'A', 'C', 'G', 'T'};
  




    /*********************************************************
     ****************** INSTANCE VARIABLES *******************
     *********************************************************/

  /**
   * Dictionary of all k-mers found in the sequencing reads.
   * Keys correspond to k-mers and values correspond to the
   * number of occurrences of each k-mer.
   */
  private Hashtable<String, Integer> k_mers;
  
  /**
   * Genome sequence, initially empty, that will be grown and assembled.
   */
  private StringBuilder genome;
  


    /***************************************************
     ****************** CONSTRUCTORS *******************
     ***************************************************/

  /**
   * Creates an <code>Assembler</code> that builds a genome assembly
   * from a file of sequencing reads.
   * 
   * @param   filename   the name of the file containing sequencing reads
   */
  public Assembler(String filename) {
    k_mers = new Hashtable<String, Integer>();
    populate_kmer_dictionary(filename);  // Populate dictionary with k-mers from sequencing reads
    genome = new StringBuilder();
    assemble();
  }
  


  /*******************************************************
    ****************** INSTANCE METHODS *******************
    *******************************************************/

  /**
   * Returns a <code>String</code> representation of the assembly, i.e.,
   * returns the assembled genome sequence.
   * 
   * @return          a <code>String</code> representation of the genome assembly
   */
  public String toString() {
    return genome.toString();
  }
  
  /**
   * Reads in a file of sequencing reads. For each k-mer in a
   * sequencing read, the k-mer is added to the dictionary and 
   * the number of occurrences of the k-mer is updated appropriately.
   * 
   * @param   filename   the name of the file containing sequencing reads
   */
  public void populate_kmer_dictionary(String filename) {
	  String sequence = SequenceOps.sequenceFromFastaFile(filename);
	  for (int i=0; i<sequence.length()-k+1; i=i+1){
		  String kmer= sequence.substring(i, i+k);
		  if (k_mers.containsKey(kmer)){
			  k_mers.put(kmer, k_mers.get(kmer)+1); 
		  }else{
			  k_mers.put(kmer, 1);
		  }
	  }
  }

  /**
   * Assembles a genome sequence based on a dictionary of k-mers.
   * As a starting point, one k-mer from the dictionary is chosen.
   * This k-mer is extended forward as long as possible, one character
   * at a time, into a growing genome sequence.
   * When the genome sequence cannot be extended further forward,
   * then the genome sequence is extended backward as long as possible,
   * one character at a time.
   */
  public void assemble() {
    if (k_mers.isEmpty()) return;  // No k-mers in dictionary
    
    // Initialize genome sequence to be one k-mer from the dictionary
    String starting_kmer = k_mers.keys().nextElement();
    k_mers.remove(starting_kmer); 
    genome.append(starting_kmer);
    extendGenomeSequenceForward();
    extendGenomeSequenceBackward();    
  }
  
  /**
   * Attempts to extend the genome sequence forward (to the right)
   * repeatedly, one character at a time. The genome sequence is
   * extended and a character added to its end if there exists a
   * nucleotide character that can be added to the end of the k-1
   * final characters of the genome sequence to form a k-mer that
   * occurs in the k-mer dictionary. If there are multiple individual
   * nucleotide characters that can be added to the final k-1 characters
   * in the genome sequence to form k-mers in the dictionary, then the
   * character resulting in the k-mer with the largest number of
   * occurrences in the dictionary is chosen. Each time the genome
   * sequence is extended by a character, the corresponding k-mer is 
   * removed from the dictionary so that it will not be used in 
   * future extensions.
   */
  
  public String maxKey(Hashtable<String, Integer> table){
	  int highest_weight = Collections.max(table.values());
	  List<String> candidates = new ArrayList<String>();
	  for (String key: table.keySet()){
		  if (table.get(key) == highest_weight)
			  candidates.add(key);
	  }  
	  if (candidates.size()>1){
		  Random rand = new Random();
			  return candidates.get(rand.nextInt(candidates.size()));
		  }else{
			  return candidates.get(0);
		  }  
  }
  public void updateKmers(Hashtable<String, Integer> candidates){
	  String selected = maxKey(candidates);
	  k_mers.put(selected, k_mers.get(selected)-1);
  		if (k_mers.get(selected)==0){
  			k_mers.remove(selected);
  }
  }
  public String findNext(String starter, boolean forward){
	  Hashtable<String, Integer> candidates = new Hashtable<String, Integer>(); 
	  for (int i=0; i<NTs.length; i++){
		  StringBuilder next = new StringBuilder();
		  if (forward){
		  next.append(starter.substring(1));
		  next.append(NTs[i]);
		  }else{
			  next.append(NTs[i]);
			  next.append(starter.substring(0,3));
		  }
		  if (k_mers.containsKey(next.toString())){
			  candidates.put(next.toString(), k_mers.get(next.toString()));
		  }  
		  }
	  
	  if (candidates.size()<1){
		  return "*";
	  }else{
		  updateKmers(candidates);
	  	  return maxKey(candidates);
	  }	 
  }
 
  public void extendGenomeSequenceForward() {
	  String starter = maxKey(k_mers);
	  genome.append(starter);
	  do{ 
		  String next = findNext(starter, true);
		  genome.append(next.substring(next.length()-1));
		  starter=next;
	}  while (starter!="*");
	}

  
  /**
   * Attempts to extend the genome sequence backward (to the left)
   * repeatedly, one character at a time. The genome sequence is
   * extended and a character added to its front if there exists a
   * nucleotide character that can be added to the front of the k-1
   * first characters of the genome sequence to form a k-mer that
   * occurs in the k-mer dictionary. If there are multiple individual
   * nucleotide characters that can be prepended to the first k-1 characters
   * in the genome sequence to form k-mers in the dictionary, then the
   * character resulting in the k-mer with the largest number of
   * occurrences in the dictionary is chosen. Each time the genome
   * sequence is extended by a character, the corresponding k-mer is 
   * removed from the dictionary so that it will not be used in 
   * future extensions.
   */
  public void extendGenomeSequenceBackward() {
	  String starter = genome.substring(0,4);
	  String next;
      do {
		  next = findNext(starter, false);
		  genome.insert(0, next.substring(0));
		  starter=next;
      }  while (next !="*");
      //To delete the "*" characters:
      genome.deleteCharAt(0); 
      genome.deleteCharAt(genome.length()-1);
      
  }
  
  
  
    /****************************************************
     ****************** CLASS METHODS *******************
     ****************************************************/



    /**
     * The <code>main</code> method executes the <code>Assembler</code>
     * on a file of sequencing reads. The assembled genome is output
     * to standard out.
     *
     * @param   args   a <code>String</code> array of any command line arguments
     */
  public static void main(String[] args) {
    
    if (args.length < 1) {
      System.err.println("When executing this program, as a command line argument please enter the");
      System.err.println("name of a file containing sequencing reads. The application will attempt");
      System.err.println("to assemble the sequencing reads into a genome sequence. Output is to stdout.");
      System.err.println("\n\tjava Assembler data/reads.txt\n");
      return;
    }
    
    Assembler a = new Assembler(args[0]);
    System.out.println(a.toString());  // Print out genome sequence to stdout
    
  }
}
