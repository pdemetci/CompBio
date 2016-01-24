	/**
	 * @author = Pinar Demetci
	 */
	import java.io.*;
	import java.util.*;
	import java.util.Collections;


	
public class SequenceOps {


		public static void main(String args[]){
	}
		
		/**
		 * TASK 1.1
		 * Turns the nucleotide sequence found in the specified FASTA file into string format.
		 * @param filepath The path to the FASTA file.
		 * @return String of the sequence in the FASTA.
		 */
		public static String sequenceFromFastaFile(String filePath){
			
			StringBuilder dnaSequence = new StringBuilder("");
			try {
			    BufferedReader lineReader = new BufferedReader(new FileReader(filePath));
			    String lineText = null;
			 
			    while ((lineText = lineReader.readLine()) != null) {
			    	if (lineText.startsWith(">")==false){ //So that the first line containing info is avoided in the string.
			           dnaSequence.append(lineText);
			           
			    	}
			    }
			    lineReader.close();
			    String dna=dnaSequence.toString();
			    return dna;
			   
			} catch (IOException ex) {
			    System.err.println(ex);
			    return "IOerror";
			}
		}
		
		/**
		 * TASK 1.2
		 * Represents the string of nucleotide sequence in FASTA format.
		 * Turns the string of sequence into a StringBuilder and adds a new line after every 60 nucleotides.
		 * @param s The string of nucleotide sequence
		 * @return String of the sequence in the FASTA format, 60 nucleotides in each line.
		 */
		public static String sequenceToFastaFormat(String s){
			StringBuilder FASTA = new StringBuilder("");
			FASTA.append(s);
			int i=60;
			while (i<s.length()){
				FASTA.insert(i, "\n");
				i=i+61;
			}
			return FASTA.toString();
		}
		
		/**
		 * TASK 1.3
		 * Returns a String corresponding to the reversed version of the specified sequence s.
		 * @param s The string of nucleotide sequence
		 * @return String of the sequence in the FASTA format, 60 nucleotides in each line.
		 */
		public static String reverse(String s){
			String reverse = new StringBuilder(s).reverse().toString();
			return reverse;
		}
		
		/**
		 * TASK 1.4
		 * Returns a String corresponding to the complement of the specified nucleotide sequence s.
		 * @param s The string of nucleotide sequence
		 * @return Complement sequence in String format.
		 */
		public static String complement(String s){

			StringBuilder complement = new StringBuilder("");

			for (int i = 0; i< s.length(); i++){
				if (s.charAt(i)=='A'){
					complement.append('T');
				}
				if (s.charAt(i)=='G'){
					complement.append('C');
				}
				if (s.charAt(i)=='C'){
					complement.append('G');
				}
				if (s.charAt(i)=='T'){
					complement.append('A');
				}
			}
			return complement.toString();
		}
		/**
		 * TASK 1.5
		 * Computes the reverse complement of a specified nucleotide.
		 * @param s The string of nucleotide sequence
		 * @return Complement sequence in String format.
		 */
		public static String reverseComplement(String s){
			String reverse = reverse(s);
			String reverseComplement = complement(reverse);
			return reverseComplement;
		}
		
		/**
		 * TASK 1.6
		 * Computes the GC content of the specified nucleotide sequence.
		 * @param s Nucleotide sequence in String format.
		 * @return GC content ratio to the whole sequence in double format.
		 */
		
		// In the yeast genomic sequence, GC content makes up around 38% (0.38297671574922004) of the sequence.
		// In the E.coli genomic sequence, GC content makes up around 50% (0.5047480343799055) of the sequence.
		// In the human 22 chromosome sequence, GC content makes up  around 48% (0.4791780185638117) of the sequence.
		
		public static double GC_content(String s){
			double GCcount=0;
			for (int i = 0; i < s.length(); i++) {
				if(s.charAt(i) == 'C' || s.charAt(i) == 'G'){
					GCcount++;
				}
			}
			return (GCcount/s.length());
			
		}
		
		/**
		 * TASK 1.7
		 * Computes the GC content of the specified nucleotide sequence.
		 * @param s Nucleotide sequence in String format.
		 * @return GC content ratio to the whole sequence in double format.
		 */
		
		//I was curious to see if there is a built in shuffle function in Java and looked up on StackOverflow.
		//Used what I saw there. This is the link: http://stackoverflow.com/questions/3981420/collections-shuffle
		public static String randomPermutation(String s){
			
		    List<Character> nucleotides = new ArrayList<Character>();
		    for(char c:s.toCharArray()){
	              nucleotides.add(c);
		        }
		    StringBuilder randomSequence = new StringBuilder(s.length());
			Collections.shuffle(nucleotides);
			for(char c:nucleotides){
				randomSequence.append(c);
			}
			return randomSequence.toString();
		    }		
			
		/**
		 * TASK 1.8
		 * Computes the GC content of the specified nucleotide sequence.
		 * @param s Nucleotide sequence in String format.
		 * @return GC content ratio to the whole sequence in double format.
		 */
		public static String randomSequence(int length, double GC_content){

			//number of G nucleotides in the sequence:
			int G_number =  (int) (length*GC_content)/2;
			//number of C nucleotides in the sequence:
			int C_number = (int) (length*GC_content)-G_number;
			//number of A nucleotides in the sequence:
			int A_number =  (int) ((length-(G_number+C_number))/2);
			//number of T nucleotides in the sequence:
			int T_number = (int) (length-(G_number+C_number+A_number));
			
			StringBuilder nucleotides = new StringBuilder("");
			int i = 0;
			for (i=0; i< G_number; i++){
				nucleotides.append('G');
			}
			for (i=0; i< C_number; i++){
				nucleotides.append('C');
			}
			for (i=0; i< A_number; i++){
				nucleotides.append('A');
			}
			for (i=0; i< T_number; i++){
				nucleotides.append('T');
			}
			String s = nucleotides.toString();
			
			return randomPermutation(s);
		}
		
		/**
		 * TASK 1.9
		 * Computes the GC content of the specified nucleotide sequence.
		 * @param s Nucleotide sequence in String format.
		 * @return GC content ratio to the whole sequence in double format.
		 */
		public static String randomSampling(String s){
			int length = s.length();
			double GC_content = GC_content(s);
			return randomSequence(length, GC_content);
		}


		/**
		 * Helper Function for Task 1.10
		 * Turns a tab-delimited .txt file into a HashMap
		 * @param filepath Location of the .txt file in String format.
		 * @return HashMap containing information in the .txt file's first 2 columns.
		 */	
		
	// I looked up the way to do it on StackOverflow and modified it as necessary.*** 
	// The link of the source I used is: 
	// http://stackoverflow.com/questions/8886103/read-from-a-text-file-into-a-hash-map-or-list

	// Note: The translation.txt file I downloaded from the website had spaces instead of tabs between the 
	// second and the third column although it was supposed to have tabs. Therefore, I deleted the last column
	// while using it so it wouldn't affect the HashMap. I am not sure if I was supposed to do it in code but I assumed 
	// it was an unintentional problem. 

	public static Map<String, String> FiletoHashMap(String filepath){
		Map<String, String> map = new HashMap<String, String>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(filepath));
			String line = "";

			while ((line = in.readLine()) != null) {
			    String parts[] = line.split("\t");
			    map.put(parts[0], parts[1]);
			    }
			in.close();
		   return map;
		} catch (IOException ex) {
		    System.err.println(ex); 
		    Map<String, String> mapError = new HashMap <String, String>();
		    mapError.put("error", "error");
		    return mapError;
		}   
	}
		
	/**
	 * TASK 1.10
	 * Depending on the input codon, looks up the corresponding amino acid from the HashMap of codons:aminoacids
	 * @param s String of codon.
	 * @return Amino acid represented by one character.
	 */

	public static char translateCodon(String s){
		Map<String, String> map = FiletoHashMap("/home/pinar/Desktop/translation.txt");
		Object value = map.get(s);
		char aminoacid = value.toString().charAt(0);
		return aminoacid;
		}

	/**
	 * TASK 1.11
	 * Translates the ORF into amino acid sequence
	 * @param s String of nucleotides in ORF.
	 * @return Amino acid sequence in String format.
	 */
	public static String translateORF(String s){
		int i=0;
		char aa;
		StringBuilder aminoacids = new StringBuilder();
		for (i=0; i<s.length()-2; i=i+3){
			String codon= s.substring(i, i+3);
			aa= translateCodon(codon);
			aminoacids.append(aa);
		}
		return aminoacids.toString();
	}
	}





