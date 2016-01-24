/*
 * author: Pinar Demetci
 */

import java.io.*;
import java.util.*;

public class Mapper {

	/*
	 * Prints out the number of reads in the file.
	 */
	
	// Helper things to store \\
	private static List<String> reads = new ArrayList<String>();
	private static Hashtable<Character, List<Integer>> ranges = new Hashtable<Character, List<Integer>>();
	private static List<Integer> occurences = new ArrayList<Integer>();
	private static List<Integer> storage=new ArrayList<Integer>();
	
	// Things to return in the end \\
	private static int noMaps;
	private static int oneMap;
	private static int multipleMaps;
	
	
	public static void main(String args[]){
		BWT b = new BWT(args[0]);
		getRanges(b);
		System.out.print("Number of total reads processed:                 ");
		System.out.println(numReads(args[1]));

		for (String r:reads){
			storage.clear();
			occurences.clear();
			lookUpMaps(b,r);
			addOcctoStorage();
		// Search the other strand, too:
			String reverseRead = SequenceOps.reverseComplement(r);
			occurences.clear();
			lookUpMaps(b,reverseRead);
			addOcctoStorage();
//		// Calculate the number of maps:
			countMaps();
			
		}		  
		
		System.out.print("Number of unmapped reads:                        ");
		System.out.println(noMaps);
		
		System.out.print("Number of reads mapping to exactly one location: ");
		System.out.println(oneMap);
		
		System.out.print("Number of reads mapping to multiple locations:   ");
		System.out.println(multipleMaps);
	}
	
	/*
	 * Function to find how many times a certain sequence read appears in the reference sequence
	 */
	public static void lookUpMaps (BWT b, String read){
		List<Integer> startRanges = ranges.get(read.charAt(read.length()-1));
		for(int i:startRanges){
			int occurence = b.getNumberOccurrencesPriorToIndex(read.charAt(read.length()-2), i);
			occurences.add(occurence);
		}
		scrapeOcc(occurences);
		int index = 2;
		while (occurences.size()!=0 && index < read.length()){
			List<Integer> rangesOfInterest = lookUpFirstChar(b,read.charAt(read.length()-index), read);
			index ++;
			lookUpLastChar(b, read.charAt(read.length()-index),rangesOfInterest, read);
		}
	}
	
	/*
	 * Looking up the first characters in rotations
	 */
	public static List<Integer> lookUpFirstChar(BWT b, char c, String read ){
		List<Integer> startRanges = ranges.get(c);
		List<Integer> rangesOfInterest = new ArrayList<Integer>();
		for (int o:occurences){
			if (o<startRanges.size()){
			rangesOfInterest.add(startRanges.get(o));
			}
		}
		return rangesOfInterest;
	}
	
	/*
	 * Looking up the last characters in rotations: in BWT transform
	 * In the end, number of occurences will define how many times the read appears
	 */
	public static void lookUpLastChar(BWT b, char c, List<Integer> rangesOfInterest, String read){
		occurences.clear();
		for(int i:rangesOfInterest){
			int occurence = b.getNumberOccurrencesPriorToIndex(c, i+1);
			occurences.add(occurence);
		}
		scrapeOcc(occurences);
	}
	/*
	 * Finally, counting whether a read occurs at all or occurs multiple times
	 */
	public static void countMaps(){
		if (storage.size()==0){
			noMaps ++;
		}else if (storage.size()==1){
			oneMap++;
		}else{
			multipleMaps++;
		}
	}
	
	
	/*
	 * Helper Functions:
	 */
	public static void addOcctoStorage(){
		for (Integer integer : new ArrayList<>(occurences)) {
			storage.add(integer);
		}
	}
	/*
	 * Read in files and output the number of reads
	 */
	public static int numReads(String filePath){
		int readNum=0;
		try {
		    BufferedReader lineReader = new BufferedReader(new FileReader(filePath));
		    String lineText = null;
		 
		    while ((lineText = lineReader.readLine()) != null) {
		    	if (lineText.startsWith(">")==false){ //So that the first line containing info is avoided in the string.
		    		reads.add(lineText);
		    		readNum++;}
		    }lineReader.close();
		    return readNum;
		 } catch (IOException ex) {System.err.println(ex); return 0;}
	}
	
	public static List<Integer> findRanges(BWT b, char c1, char c2){
		List<Integer> output = new ArrayList<Integer>();
		output.add(b.getNumberCharactersLessThan(c1));
		if (c1 == 'T'){
			output.add(b.getLength()-1);
		}else{
		output.add(b.getNumberCharactersLessThan(c2)-1);
		}
		for (int i= output.get(0)+1; i<output.get(1); i++){
			output.add(i);
		}
		Set<Integer> rangeSet = new HashSet<>();
		rangeSet.addAll(output);
		output.clear();
		output.addAll(rangeSet);
		Collections.sort(output);
		return output;
	}
	/*
	 * Storing in the beginning and looking up would be fast
	 */
	public static void getRanges(BWT b){
		//For As:
		List<Integer> As = findRanges(b, 'A', 'C');
		ranges.put('A', As);
		//For Cs:
		List<Integer> Cs = findRanges(b, 'C', 'G');
		ranges.put('C', Cs);
		//For Gs:
		List<Integer> Gs = findRanges(b, 'G', 'T');
		ranges.put('G', Gs);
		//For Ts:
		List<Integer> Ts = findRanges(b, 'T', '$');
		ranges.put('T', Ts);

	}
	/*
	 * Since the length of occurences list will define how many times a reading appears,
	 * we do not want any duplicates in the list.
	 */
	public static List<Integer> scrapeOcc (List<Integer> occurences){
		Set<Integer> occurencesSet = new HashSet<>();
		occurencesSet.addAll(occurences);
		occurences.clear();
		occurences.addAll(occurencesSet);
		for (Integer integer : new ArrayList<>(occurences)) {
		    if (integer == 0) { //0 occurences means none. We should remove to prevent misleading results. 
		        occurences.remove(integer);
		    }
		}
		return occurences;
	}
	}