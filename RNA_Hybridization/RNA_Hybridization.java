/*
 * author = Pinar Demetci
 */
public class RNA_Hybridization {
	
	public String seq1;
	public String seq2;
	public String upper = "";
	public String lower = "";
	public String hybridization = "";
	public Energies e;
	public double[][] H;
	public double minimumEnergy = Double.MAX_VALUE;
	public static int minimumColumn = -1; 
	public static int minimumRow = -1;
	
	/*
	 * Constructor for RNA_Hybridization
	 */
	public RNA_Hybridization(String pathToSequence1, String pathToSequence2) {
		seq1 = SequenceOps.sequenceFromFastaFile(pathToSequence1);
		seq2 = SequenceOps.sequenceFromFastaFile(pathToSequence2);
		seq2=SequenceOps.reverse(seq2); //Reversing one of the sequences
		e = new Energies("energyData", seq1, seq2);
		H = new double[seq1.length()][seq2.length()];
	}
	
	public static void main(String args[]){
		RNA_Hybridization r = new RNA_Hybridization(args[0], args[1]);
        r.computeHybridization(); 
        System.out.println(r.upper);
        System.out.println(r.hybridization);
        System.out.println(r.lower);
        System.out.println(r.minimumEnergy);
	}
	
	public void computeHybridization(){
        fillHybridizationMatrix();
        backTracking(minimumRow, minimumColumn);
	}
	
	/*
	 * 
	 */
	private void fillZeros(){
		//Initializing the matrix values
		for(int i=0; i<seq1.length();i++){
			H[i][0]=0.0;
		}
		for(int j=0; j<seq2.length();j++){
			H[0][j]=0.0;
		}
	}
	
	public void fillHybridizationMatrix(){
		fillZeros();
		//Filling out the matrix		
		for (int i =1; i<seq1.length();i++){
			for (int j=1; j< seq2.length(); j++){			
				//for stacking region:
				double case1= e.stacking(i-1, j-1, i,j)+ H[i-1][j-1];
				
				//for bulge loop in the first sequence:
				double case2=Double.MAX_VALUE;
				for (int k=0; k<i-1;k++){
					double bulge1=e.bulge(k, j-1, i, j)+ H[k][j-1];
					if (bulge1<case2){
						case2=bulge1;
					}
				}
				//for bulge loop in the second sequence:
				double case3=Double.MAX_VALUE;
				for (int k=0;k<j-1;k++){
					double bulge2=e.bulge(i-1, k, i, j)+H[i-1][k];
					if (bulge2<case3) {
						case3=bulge2;
					}
				}
				
				//for interior loop:
				double case4=Double.MAX_VALUE;
				for(int k=0; k<i-1;k++){
					for (int l=0;l<j-1;l++){
						double loop = e.interior(k, l, i, j)+ H[k][l];
						if (loop<case4){
							case4 = loop;
						}
					}
				}
				
				H[i][j]=Math.min(case1, Math.min(case2, Math.min(case3, Math.min(case4, 0.0))));
				
		// Update the minimumEnergy & save the columns
		if (H[i][j]<minimumEnergy){
			minimumEnergy=H[i][j];
			minimumColumn=j;
			minimumRow=i;
		}
	}
	}
}
	
	
	/*
	 * A function that generates a String of x number of space character, ' '. For backtracking
	 */
	private String spaces(int x){
		StringBuilder spaces = new StringBuilder();
		for (int i=0; i<x; i++){
			spaces.append(' ');
		}
		return spaces.toString();
	}
	
	/*
	 * A function that generates a String of x number of gap character, '-'. For backtracking
	 */
	private String gaps(int x){
		StringBuilder gaps = new StringBuilder();
		for (int i=0; i<x; i++){
			gaps.append('-');
		}
		return gaps.toString();
	}
	
	/*
	 * Backtracking on the H matrix. Checks to see where the value comes from and 
	 * modifies "upper", "hybridization" and "lower" Strings for outputting the 
	 * hybridization
	 */
	public void backTracking(int i, int j){
		//check to see which situation (i,h) falls into
		
		//Stacking Region
		double case_stacking = e.stacking(i-1, j-1, i,j) + H[i-1][j-1]; 
		
		double case_Bulge1= Double.MAX_VALUE;
		double case_Bulge2= Double.MAX_VALUE;
		double case_interiorLoop = Double.MAX_VALUE;
		
		int Bulge1Index=-1;
		int Bulge2Index=-1;
		int interior1=-1;
		int interior2=-1;

		
		for (int k=0; k<i-1; k++){
			//Bulge in String 1
			double bulge1 = e.bulge(k,j-1,i,j)+H[k][j-1];
			if (bulge1< case_Bulge1){
				case_Bulge1 = bulge1;
				Bulge1Index = k;
			}
			
			double bulge1 = bulge1Engergy(i,j)[0];
			int Bulge1Index = 
			
			for (int l=0; l<j-1; l++){
				//Bulge in String 2
				double bulge2 = e.bulge(i-1,k,i,j)+H[i-1][k];
				if (bulge2< case_Bulge2){
					case_Bulge2 = bulge2;
					Bulge2Index = k;
				}
				//interior Loop
				double intLoop = e.interior(k,l,i,j)+H[k][l];
				if (intLoop< case_interiorLoop){
					case_interiorLoop = intLoop;
					interior1 = k;
					interior2=l;
			}
			}
		}
		
		 if (H[i][j] == 0.0) {
	           int max = Math.max(i, j);
	           StringBuilder zero = new StringBuilder();
	           upper = (zero.append(spaces(max - i)).append(seq1.substring(0, i + 1)).append(upper)).toString();
	           zero.setLength(0);
	           lower = (zero.append(spaces(max - j)).append(seq2.substring(0, j + 1)).append(lower)).toString();
	           zero.setLength(0);
	           hybridization=(zero.append(spaces(max)).append("|").append(hybridization)).toString();
	        }
		
		if(H[i][j]==case_stacking){
			endsStackingRegion(i,j);
		}
		
		else if(H[i][j]==case_Bulge1){
			endsBulge1(i,j, Bulge1Index);
		}
		
		else if(H[i][j]==case_Bulge2){
			endsBulge2(i,j, Bulge2Index);
		}
		
		else if(H[i][j]==case_interiorLoop){
			endsInteriorLoop(i,j, interior1, interior2);
		}
		
	}

	
	/*
	 * Helper Functions for backtracking
	 * Functions are from the backtracking algorithm on http://cs.wellesley.edu/~cs313/projects/project7/project7.html
	 */
	
	private void endsStackingRegion(int i, int j){
		StringBuilder stacking = new StringBuilder();
		upper = (stacking.append(seq1.charAt(i)).append(upper)).toString();
		stacking.setLength(0);
		hybridization= (stacking.append("|").append( hybridization)).toString();
		stacking.setLength(0);
		lower= (stacking.append(seq2.charAt(j)).append(lower)).toString();
		backTracking(i-1 , j-1);
	}
	
	private void endsBulge1(int i, int j, int Bulge1Index){
		StringBuilder bulge1 = new StringBuilder();
		upper= (bulge1.append(seq1.substring(Bulge1Index+1, i+1)).append(upper)).toString();
		bulge1.setLength(0);
		hybridization=(bulge1.append(spaces(i-Bulge1Index-1)).append("|").append( hybridization)).toString();
		bulge1.setLength(0);
		lower=(bulge1.append(gaps(i-Bulge1Index-1)).append(seq2.charAt(j)).append(lower)).toString();
		backTracking(Bulge1Index , j-1);
	}
	
	private void endsBulge2(int i, int j, int Bulge2Index){
		StringBuilder bulge2 = new StringBuilder();
		upper=(bulge2.append(gaps(j-Bulge2Index-1)).append(seq1.charAt(i)).append(upper)).toString();
		bulge2.setLength(0);
		hybridization=(bulge2.append(spaces(j-Bulge2Index-1)).append("|").append( hybridization)).toString();
		bulge2.setLength(0);
		lower=(bulge2.append(seq2.substring(Bulge2Index+1, j+1)).append(lower)).toString();
		backTracking(i-1, Bulge2Index);
	}
	
	private void endsInteriorLoop(int i, int j, int interior1, int interior2){
		 int maxInterior = Math.max(i - interior1 - 1, j - interior2 - 1);
		 StringBuilder loop = new StringBuilder();
         upper = (loop.append(seq1.substring(interior1 + 1, i)).append(gaps(Math.max(maxInterior - (i - interior1 - 1), 0))).append(seq1.charAt(i)).append(upper)).toString();
         loop.setLength(0);
         hybridization=(loop.append(spaces(maxInterior)).append("|").append( hybridization)).toString();
         loop.setLength(0);
         lower=(loop.append(seq2.substring(interior2 + 1, j)).append(gaps(Math.max(maxInterior - (j - interior2 - 1), 0))).append(seq2.charAt(j)).append(lower)).toString();
         backTracking(interior1, interior2);
	}
}
