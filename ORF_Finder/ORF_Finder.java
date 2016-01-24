/**
 * Pinar Demetci
 * Computational Biology - Project #2 - Task 1
 * 24 September 2015
 * Necessary corrections are made to SequenceOps file and uploaded to the server, as well. 
 */

import java.io.BufferedReader; 
import java.io.FileReader;
import java.io.IOException;
import java.util.*; 

public class ORF_Finder {
	public static void main(String args[]){
	System.out.println("Calculating...");
	String sequence = SequenceOps.sequenceFromFastaFile(args[0]);
	String randomSequence = SequenceOps.randomSequence(sequence.length(), SequenceOps.GC_content(sequence));
	List<String> ORF_List = ORF_List(sequence);
	List<String> Aminoacids = ORF_translation(ORF_List);
	List<String> RandomORF_List = ORF_List(randomSequence);
	List<String> Random_Aminoacids = ORF_translation(RandomORF_List);
	System.out.println("**************************************************************");
	System.out.println("***** Genomic Sequence Read-In From Fasta Formatted File *****");
	System.out.println("**************************************************************");
	System.out.println("Number of ORFs is:"+ "\t" + ORF_Number(ORF_List));
	System.out.println("Number of ORFs at least 300nt is:"+"\t" + ORF_300n (ORF_List));
	System.out.println("Number of ORFs with specified amino acid distribution is:"+"\t" + ORF_Dist(Aminoacids));
	System.out.println("Number of ORFs with Kozak sequence is:" + "\t" + Kozak_find(sequence));
	System.out.println("Number of ORFs with likely protein coding codons is:" + "\t" + Coding_ORF((ORF_List)));
	System.out.println("**************************************************************");
	System.out.println("*** Random Sequence With Same Length And GC-Content As Sequence Read-In From File ***");
	System.out.println("**************************************************************");
	System.out.println("Number of ORFs is:"+ "\t" + ORF_Number(RandomORF_List));
	System.out.println("Number of ORFs at least 300nt is:"+"\t" + ORF_300n(RandomORF_List));
	System.out.println("Number of ORFs with specified amino acid distribution is:" + "\t" + ORF_Dist(Random_Aminoacids));
	System.out.println("Number of ORFs with Kozak sequence is:" + "\t" + Kozak_find(randomSequence));
	System.out.println("Number of ORFs with likely protein coding codons is:" + "\t" + Coding_ORF((RandomORF_List)));
	}

/*Helper Function
 * 	Takes in a dna sequence in String format and returns a HashMap of (initial_index(start codon):initial_index(stop codon))
 */
	public static Map<Integer,Integer> ORF_Index(String dna){
		Map<Integer,Integer> ORF = new HashMap<Integer,Integer>();
		for (int i=0; i<dna.length()-2;i++){
			if ((dna.substring(i,i+3)).equals("ATG")){
				for (int n=i; n<dna.length()-2; n=n+3){
					if ((dna.substring(n,n+3)).equals("TGA") | (dna.substring(n,n+3)).equals("TAG") | (dna.substring(n,n+3)).equals("TAA")){
						ORF.put(i,n+3);
						break;
					}
				}
			}
			
		}
		return ORF;
	}
/*Helper Function
 * 	Takes in a string of dna sequence and returns a list of ORFs
 */
	public static List<String> ORF_List(String dna){
		List<String> ORF = new ArrayList<String>();
		for (int i=0; i<dna.length()-2;i++){
			if ((dna.substring(i,i+3)).equals("ATG")){
				for (int n=i; n<dna.length()-2; n=n+3){
					if ((dna.substring(n,n+3)).equals("TGA") | (dna.substring(n,n+3)).equals("TAG") | (dna.substring(n,n+3)).equals("TAA")){
						ORF.add(dna.substring(i,n+3));
						break;
					}
				}
			}
			
		}
		return ORF;
	}
	
/*Task 1
 * Returns the number of ORFs in a given dna sequence
 */
	public static int ORF_Number(List<String> ORF){
		return ORF.size();
	}
	
/*Task 1
 * Takes in a list of ORFs and computes the number of ORFs at least 300 nucleotide.
 * 
 */
	public static int ORF_300n (List<String> ORF){
		int count = 0;
		for (int i=0; i< ORF.size();i++){
			if (ORF.get(i).length()>=300){
				count++;
			}
	}
		return count;
	}
	
	/*
	 * Helper Function	
	 *Translates ORF into aminoacid sequence
	 *Uses the same translation.txt file from Project 1, except U = T.
	 */
	
	public static List<String> ORF_translation(List<String> ORFs){
		List <String> AminoAcids= new ArrayList <String>();
		for (int i=0; i<ORFs.size(); i++){
			AminoAcids.add(SequenceOps.translateORF(ORFs.get(i)));
		}
		return AminoAcids;
	}
	
	/*
	 * Task 1
	 * Takes in list of aminoacids corresponding to ORFs and returns the number of ORFs that code for proteins containing:
	 *less than 4% cysteines and less than 4% tryptophans and more than 6% alanines and more than 6% leucines
	 */
	public static int ORF_Dist(List<String> AminoAcids){
		int ORFs=0;
		for (int i=0; i<AminoAcids.size(); i++){
			List<String> aa = new ArrayList<String>();
			aa.add(AminoAcids.get(i));
			int C_count=0;
			int W_count=0;
			int A_count=0;
			int L_count=0;
			int size = AminoAcids.get(i).length();
			for (String a:aa){	
				if (a=="C"){
					C_count++;
				}
				if (a=="W"){
					W_count++;
				}
				if (a=="A"){
					A_count++;
				}
				if (a=="L"){
					L_count++;
				}
			if ((C_count/size < 0.04) && (W_count/size < 0.04) && (A_count/size > 0.06) && (L_count/size > 0.06)){
				ORFs++;
			}
				
			}
		}
		return ORFs;
	}
	/*Task 1
	 *Takes in a string of dna and returns the number of Kozak sequences. 
	 */
	public static int Kozak_find(String dna){
		int Kozak = 0;
		Map<Integer, Integer> ORFs = ORF_Index(dna);
		for (int key : ORFs.keySet()){
			if ((dna.substring(key-3,key+4).equals("ACCATGG")) | dna.substring(key-3,key+4).equals("GCCATGG")){
				Kozak++;
			}
		}
		return Kozak;
	}

	/* Helper Function
	 * Turns a tab-deliminated file into a HashMap. Takes the filepath and returns the HashMap
	 */
	public static Map<String, String> FiletoHashMap(String filepath, int first, int second){
		Map<String, String> map = new HashMap<String, String>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(filepath));
			String line = "";

			while ((line = in.readLine()) != null) {
			    String parts[] = line.split("\t");
			    map.put(parts[first], parts[second]);
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
	/* Task 1 - Coding ORFs
	 * Takes a list of ORFs and returns the number of ORFs that are likely coding for proteins.
	 */
	public static int Coding_ORF(List<String> ORF){
		int Coding_ORF=0;
		Map<String, String> CodingMap = FiletoHashMap("/home/pinar/Desktop/codonUsage.txt", 0, 1);
		Map<String, String> NonCodingMap = FiletoHashMap("/home/pinar/Desktop/codonUsage.txt", 0, 2);
		for (String a: ORF){
			double coding_score=1;
			double noncoding_score=1;
			for (int i =0; i<a.length(); i=i+3){
				double coding_value = Double.parseDouble(CodingMap.get(a.substring(i, i+3)));
				coding_score *= coding_value;
				double noncoding_value = Double.parseDouble(NonCodingMap.get(a.substring(i, i+3)));
				noncoding_score *= noncoding_value;
			}
			if (coding_score - noncoding_score >0){
				Coding_ORF++;
			}
		}
		return Coding_ORF;
	}
	
}


