public class RNA_Basepair {
  
  // Instance variables
  private String sequence;
  private int[][] M;
  
  
  
  // Constructor
  public RNA_Basepair(String fileName) {
    sequence = SequenceOps.sequenceFromFastaFile(fileName).toUpperCase().replace('T', 'U');
    M = new int[sequence.length()][sequence.length()];
  }
  
  
  
  // Public instance methods
  
  // For the RNA sequence, fill in the dynamic programming table, M, to compute
  // the structure that has the largest possible number of basepairings
  public void computeStructureWithMaximumBasepairs() {
    for (int i=sequence.length()-1; i>=0; i--) {
      for (int j=i; j<sequence.length(); j++) {
        
        
        if (j-i < 2)  // If subsequence has length 2 or less then number of basepairs is zero
          M[i][j] = 0;  // Fill in first two diagonals of matrix with zero
        else {  // Compute maximum number of basepairs in subsequence from i to j, inclusive
          int case_A = M[i+1][j-1] + energy_basepair(i,j);  // i and j basepair
          int case_B = M[i+1][j];  // i doesn't pair
          int case_C = M[i][j-1];  // j doesn't pair
          
          int case_D = Integer.MIN_VALUE;  // bifurcation
          // Add your code here to handle case_D, bifurcation
          int bifurcation_MAX = 0;
          
          int k = i;
          int l = i+1;
          while ((k<j) && (l < j)) {
             if (M[i][k]+M[l][j]> bifurcation_MAX){
                bifurcation_MAX = M[i][k]+M[l][j];
              }
             k++;
             l++;
          }
          case_D=bifurcation_MAX;
          
          M[i][j] = max(case_A, case_B, case_C, case_D);  // Fill in dynamic programming table
        }
      }
    }
  }
  
  // Return the number of basepairs in the structure with the maximum number of basepairs
  public int getMaximumNumberBasepairs() {
    return M[0][sequence.length()-1];
  }
  
  // Backtrack through dynamic programming table to create the *structure*
  // with the maximum number of basepairings.
  // Return the structure as a String
  public String getStructure() {
    return getStructure(0, sequence.length()-1);
  }
  
  
  
  // Private instance methods
  
  // Returns 1 if nucleotides at indices i and j basepair (G:C, C:G, A:U, U:A).
  // Returns 0 otherwise.
  private int energy_basepair(int i, int j) {
    if (((sequence.charAt(i) == 'G') && (sequence.charAt(j) == 'C')) ||
        ((sequence.charAt(i) == 'C') && (sequence.charAt(j) == 'G')) ||
        ((sequence.charAt(i) == 'A') && (sequence.charAt(j) == 'U')) ||
        ((sequence.charAt(i) == 'U') && (sequence.charAt(j) == 'A')))
      return 1;
    else return 0;
  }
  
  // Returns the maximum of four integers
  private int max(int a, int b, int c, int d) {
    return Math.max(Math.max(a, b), Math.max(c, d));
  }
  
  // Return a String of the specified length consisting only of "period" characters
  private String periods(int length) {
    StringBuilder sb = new StringBuilder();
    for (int i=0; i<length; i++)
      sb.append('.');
    return sb.toString();
  }
  
  // Backtrack through dynamic programming table to create the *structure*
  // with the maximum number of basepairings.
  // Recursively calculate the optimal structure for the subsequence between i and j, inclusive.
  // Return the structure as a String
  private String getStructure(int i, int j) {
    if (M[i][j] == 0) return periods(j-i+1);
    
    // Determine where current table entry came from (i.e., backtrack)
    int case_A = M[i+1][j-1] + energy_basepair(i,j);
    int case_B = M[i+1][j];
    int case_C = M[i][j-1];
    int case_D = Integer.MIN_VALUE;
    int case_D_index = -1;
    for (int k=i+1; k<j; k++) {
      int temp = M[i][k] + M[k+1][j];
      if (temp > case_D) {
        case_D = temp;
        case_D_index = k;
      }
    }
    
    if (case_A == M[i][j]) {
      if (energy_basepair(i,j) == 1) return "(" + getStructure(i+1, j-1) + ")";
      else return "." + getStructure(i+1, j-1) + ".";
    } else if (case_B == M[i][j])
      return "." + getStructure(i+1, j);
    else if (case_C == M[i][j])
      return getStructure(i, j-1) + ".";
    else  // case_D
      return getStructure(i, case_D_index) + getStructure(case_D_index+1, j);
  }
  
  // Returns a String represenation of a 2D integer array
  private String arrayToString(int[][] a) {
    if (a.length == 0) return "";
    if (a[0].length == 0) return "";
    if (a.length != sequence.length()) {
      System.err.println("Error - when converting 2D array to String, the array length must be the same as the length of the RNA sequence.");
    }
    
    StringBuilder sb = new StringBuilder();
    for (int i=sequence.length()-1; i>=0; i--) {
      sb.append(sequence.charAt(i));
      for (int j=0; j<sequence.length(); j++) {
        if (j < i) sb.append("\t" + "X");
        else sb.append("\t" + a[i][j]);
      }
      sb.append("\n");
    }
    for (int i=0; i<sequence.length(); i++)
      sb.append("\t" + sequence.charAt(i));
    sb.append("\n");
    return sb.toString();
  }
  
  
  
  // Main method
  public static void main(String[] args) {
    
    if (args.length < 1) {
      System.err.println("\nTo execute this program, one command line argument is required corresponding to the name of a FASTA file containing an RNA sequence. For example,\n");
      System.err.println("\tjava RNA_Basepair sequences/tRNA.txt\n");
      return;
    }
    
    RNA_Basepair r = new RNA_Basepair(args[0]);
    r.computeStructureWithMaximumBasepairs();
    //System.out.println(r.arrayToString(r.M));
    System.out.println("\nStructure with maximum number of basepairs:\n");
    System.out.println("\t" + r.getStructure());
    System.out.println("\t" + r.sequence + "\n");
    
    System.out.println("Number of basepairs in structure:\t" + r.getMaximumNumberBasepairs() + "\n");
  }
  
}
