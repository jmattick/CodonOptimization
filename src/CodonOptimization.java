import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * <h1>Codon Optimization</h1>
 * 
 * <p>Program reduces GC content in a DNA sequence without changing the protein sequence. User can set maximum GC% parameter.
 * Program avoids reducing the frequencies of the codons. 
 * Human codon frequency table adapted from GeneScript: <a href = "https://www.genscript.com/tools/codon-frequency-table">https://www.genscript.com/tools/codon-frequency-table</a> </p>
 * <p>Program assumes codons are sorted from highest frequency/thousand to lowest in CSV input. </p>
 * <p>DNA sequence file should only contain sequence. </p>
 * @author Jessica Mattick <a target = "_blank" href = "https://github.com/jmatty16">https://github.com/jmatty16</a>
 * @version 1.0
 * @since 2018-12-09
 * @see <a href = "Sequence.html">Sequence</a>
 * @see <a href = "AminoAcidKey.html">AminoAcidKey</a>
 * @see <a href = "Codon.html">Codon</a>
 *
 */

public class CodonOptimization {
/**
 * Responsible for IO and calling methods 
 * @param args
 * @throws IOException 
 */
	public static void main(String[] args) throws IOException {
		//Start options
		System.out.println("\n********************************************\n*****     Codon Optimization Tool      *****\n***** Lower GC% and maximize frequency *****\n********************************************");
		String seqName = "";
		System.out.println("Input codon key: \tAminoAcidFrequencies.csv");
		System.out.println("Input sequence file: \tCOinput.txt");
		System.out.println("Codon key species: \tHuman");
		System.out.println("Please enter name of sequence: ");
		Scanner ip = new Scanner(System.in);
		seqName = ip.next();
		System.out.println("Default maximum GC content = 58%. Change? (y/n)");
		String changeGCParam = ip.next();
		double maxGC = 0.58;
		if(changeGCParam.equals("y")) {
			System.out.println("Enter max GC%: ");
			maxGC = ip.nextDouble()/100;
		}
		ip.close();
		String message = ""; //message to send if max GC okay
		
		
		//input and format codon key
		ArrayList<AminoAcidKey> codonKey = formatCodonKeyInput(); // calls formatCodonKeyInput and stores key to codonKey variable
		
		//input and format dna sequence
		String dnaSeq = inputSequence();
		System.out.println("Original Sequence: \n");
		Sequence originalSeq = new Sequence(dnaSeq, codonKey);
		originalSeq.updateSeq(dnaSeq);
		originalSeq.printSeqInfo(); //Display original seq with info in console
		double[] originalFreqArr = originalSeq.returnFreqArr();
		
		
		//Begin Output
		FileOutputStream fileByteStream = null;
		FileOutputStream fileByteStream2 = null;
		PrintWriter outFS = null;
		PrintWriter outFS2 = null;
		fileByteStream = new FileOutputStream(seqName+"_outputSeq.fasta");
		fileByteStream2 = new FileOutputStream(seqName+"_outputInfo.txt");
		outFS = new PrintWriter(fileByteStream);
		outFS2 = new PrintWriter(fileByteStream2);
		outFS.println(">"+seqName+" Original Sequence");
		outFS.println(originalSeq.getDnaSeq());
		outFS2.println("Original Sequence");
		outFS2.println("\tGC content: " + originalSeq.getGcContent());
		outFS2.println("\tLowest Freq: " + originalSeq.getLowestFreq());
		outFS2.println("\tAverage Freq: " + originalSeq.getAvgFreq());
		
		outFS2.println();
		
		//Switch out codons first pass: Will not lower frequency
		Sequence newSeq = originalSeq; //create new Sequence
		ArrayList<Codon> tempArr = newSeq.getSeqArr(); //array of codons in seq
		for (int i = 0; i < newSeq.getSeqArr().size(); i++) { // loop through every codon, replace codon with lower GC and equal or greater frequency codon if available
			Codon replacement = findCodon(replaceCodon(newSeq.getSeqArr().get(i).getDnaSeq(), codonKey), codonKey);
			tempArr.set(i, replacement);
			newSeq.updateSeq(tempArr);
		}
		
		
		
		//Output first optimized sequence
		outFS.println(">"+seqName+" Optimized Sequence");
		outFS.println(newSeq.getDnaSeq());
		outFS2.println("Optimized Sequence");
		outFS2.println("\tGC content: " + newSeq.getGcContent());
		outFS2.println("\tLowest Freq: " + newSeq.getLowestFreq());
		outFS2.println("\tAverage Freq: " + newSeq.getAvgFreq());
		System.out.println("\nOptimized Sequence: \n");
		newSeq.printSeqInfo();
		double[] newFreqArr = newSeq.returnFreqArr();
		
		
		outFS2.println();
		System.out.println();
		
		//Switch out codons second pass: Will lower frequency if also lowering GC content
			Sequence newSeq2 = newSeq;
			ArrayList<Codon> tempArr2 = newSeq2.getSeqArr();
			//make array of difference in frequency between original codon and the next best option
			double[] freqDiffArr = new double[newSeq2.getSeqArr().size()];
			for (int i = 0; i < newSeq2.getSeqArr().size(); i++) {
				freqDiffArr[i] = getFreqDifferences(newSeq2.getSeqArr().get(i).getDnaSeq(), codonKey);
			}
			
			//make codon switches until GC content is low enough
			double condition = newSeq2.getGcContent();
			do {
				int index = findMinDiffIndex(freqDiffArr); //finds index of codon that has the lowest reduction in frequency if codon is switched
				if (index == -1) { //ends do while loop prematurely when out of switches
					message = "Maximum GC% set too high.\n"; //message to send if loop ends before reaching set GC%
					break;
				}
				
				freqDiffArr[index] = 1001.0; //Set difference of switched codons in array to impossible value to avoid repeating
				Codon replacement = findCodon(replaceCodonRepeat(newSeq2.getSeqArr().get(index).getDnaSeq(), codonKey), codonKey); //call 2nd pass replace codon method
				tempArr2.set(index, replacement);
				newSeq2.updateSeq(tempArr2); // update sequence stats each iteration
				condition = newSeq2.getGcContent(); //update condition
				
			} while (condition > maxGC);
			
			//Print 2nd pass sequence to console
			System.out.println("2nd Pass Optimized Sequence: \n\n\tWarning: Sequence may have lower frequency. Determine if frequency is sufficient for application.\n");
			newSeq2.printSeqInfo();
			System.out.println();
			
			double[] new2FreqArr = newSeq2.returnFreqArr();
			
			outFS2.println();
			
			//Print 2nd pass sequence to output
			outFS.println(">"+seqName+" 2nd Pass Sequence Check Frequency Before Using");
			outFS.println(newSeq2.getDnaSeq());
			outFS2.println("2nd Pass Sequence Check Frequency Before Using");
			outFS2.println("\tGC content: " + newSeq2.getGcContent());
			outFS2.println("\tLowest Freq: " + newSeq2.getLowestFreq());
			outFS2.println("\tAverage Freq: " + newSeq2.getAvgFreq());
			
			outFS2.println("Original All Freq: ");
			for (int i = 0; i < originalFreqArr.length; i++) {
				outFS2.println(originalFreqArr[i]);
			}
			outFS2.println("Optimized All Freq: ");
			for (int i = 0; i < newFreqArr.length; i++) {
				outFS2.println(newFreqArr[i]);
			}
			outFS2.println("2nd Pass Optimized All Freq: ");
			for (int i = 0; i < new2FreqArr.length; i++) {
				outFS2.println(new2FreqArr[i]);
			}
			
			//flush and close output
			outFS.flush();
			fileByteStream.close();
			outFS2.flush();
			fileByteStream2.close();
		
			//End message
			System.out.print(message);//if GC max is too low, print message
			System.out.println("Sequences output to \""+seqName + "_outputSeq.fasta\"");
			System.out.println("Sequence info output to \""+seqName + "_outputInfo.txt\"");
			
	}
	/**
	 * Formats input sequence and returns ArrayList of type Codon.
	 * Checks that sequence contains only A, T, G, and C and that the sequence is divisible by 3.
	 * @return string of sequence in all uppercase
	 * @throws IOException
	 */
	public static String inputSequence() throws IOException {
		FileInputStream fileByteStream = null;
		Scanner inFS = null;
		String dnaSeq = "";
//		fileByteStream = new FileInputStream("C:\\Users\\jmatt\\eclipse-workspace\\CodonOptimization\\src\\COinput.txt");
		fileByteStream = new FileInputStream("COinput.txt");

		inFS = new Scanner(fileByteStream);
		dnaSeq = inFS.next().toUpperCase();//inputs dna seq and makes sequence uppercase
		if(!dnaSeq.matches("[ACGT]+")) {
			dnaSeq = "";
			System.out.println("Error: Sequence contains invalid characters");
		}
		if(dnaSeq.length() % 3 != 0) {
			System.out.println("Error: Sequence is not divisible by 3");
		}
		
		
		fileByteStream.close();
		inFS.close();
		return dnaSeq;
	}
	/**
	 * Converts a String containing a DNA sequence into an ArrayList of Codons
	 * @param dnaSeq String containing DNA sequence
	 * @param codonKey Formatted Arraylist<AminoAcidKey>
	 * @return ArrayList<Codon>
	 */
	public static ArrayList<Codon> formatSeqArr(String dnaSeq, ArrayList<AminoAcidKey> codonKey) {
		ArrayList<Codon> seqArr = new ArrayList<Codon>();
		for (int i = 3; i <= dnaSeq.length(); i = i+3) {
			
				seqArr.add(findCodon(dnaSeq.substring(i -3, i), codonKey));
			
		}
		return seqArr;
	}
	
	/**
	 * Parses codon frequency csv file and returns an Arraylist of AminoAcidKey
	 * @return ArrayList<AminoAcidKey>
	 * @throws IOException
	 */
	public static ArrayList<AminoAcidKey> formatCodonKeyInput() throws IOException {
		ArrayList<AminoAcidKey> key = new ArrayList<AminoAcidKey>(); //initialize reference key
		FileInputStream fileByteStream2 = null;
		Scanner inFS2 = null;
		fileByteStream2 = new FileInputStream("AminoAcidFrequencies.csv");

		inFS2 = new Scanner(fileByteStream2);
		inFS2.useDelimiter(",|\\s"); //to parse csv
		inFS2.nextLine();
		
		do {
			
		String codTemp = inFS2.next(); // gets codon sequence
		String aaTemp = inFS2.next(); //gets amino acid name
		double freqTemp = Double.parseDouble(inFS2.next()); //gets amino acid frequency (frequency/thousand)
		double gcTemp = calcGC(codTemp); //calculates CG fraction and adds to to a double variable
		inFS2.nextLine(); 
		Codon tempCodon = new Codon(codTemp, aaTemp, gcTemp, freqTemp); //creates codon variable to add to AminoAcidKey
		
		//Initializes first key item
		if(key.size() == 0) {
			key.add(new AminoAcidKey(aaTemp));
			key.get(0).addCodon(tempCodon);
		} else {
		//adds codons to key items
		for(int i = 0; i < key.size(); i++) {
			//adds codon to existing item
			if (key.get(i).getName().equals(aaTemp)) {
				key.get(i).addCodon(tempCodon);
				break;
			//Creates new AminoAcidKey item and adds codon to it
			} else if (!key.get(i).getName().equals(aaTemp) && i == (key.size() - 1)) {
				key.add(new AminoAcidKey(aaTemp));
				key.get(i+1).addCodon(tempCodon);
				break;
			}
		}
		}
		} while(inFS2.hasNextLine());
		
		fileByteStream2.close();
		inFS2.close();
		return key;
		
	}
	/**
	 * Returns DNA sequence string from ArrayList<Codon>
	 * @param seqArr ArrayList<Codon>
	 * @return String
	 */
	public static String getDNASeq(ArrayList<Codon> seqArr) {
		String res = "";
		for (Codon item: seqArr) {
			res += item.getDnaSeq();
		}
		return res;
	}
	/**
	 * Returns the fraction of g and c in a sequence
	 * @param seq String of DNA sequence
	 * @return double - fraction of g and c in a sequence
	 */
	public static double calcGC(String seq) {
		double gc = 0;
		for (int i = 0; i < seq.length(); i++) {
			if(seq.charAt(i) == 'G' || seq.charAt(i) == 'C') gc++;
		}
		return gc/seq.length();
	}
	/**
	 * Calculate Average Frequency in an ArrayList of Codons
	 * @param list ArrayList<Codon>
	 * @return double - Average Frequency
	 */
	public static double calcAvgFreq(ArrayList<Codon> list) {
		double sum = 0;
		for (Codon item: list) {
			sum += item.getFrequency();
		}
		return sum / list.size();
	}
	/**
	 * Recursive method to return the index at which the max frequency occurs in an ArrayList of Codons
	 * @param list ArrayList<Codon> to search through
	 * @param n int index to start searching
	 * @return int - index of codon with the max frequency
	 */
	public static int findMaxFreqIndex(ArrayList <Codon> list, int n) {
		if (n == list.size() -1) {
			return n;
		} 
		
		int index = findMaxFreqIndex(list, n + 1);
		return list.get(n).getFrequency() > list.get(index).getFrequency() ? n : index;
		
	}
	/**
	 * Recursive method to return the index at which the min frequency occurs in an ArrayList of Codons
	 * @param list this ArrayList of type codon to search through
	 * @param n the index to start searching
	 * @return
	 */
	public static int findMinFreqIndex(ArrayList<Codon> list, int n) {
		if (n == list.size()-1) {
			return n;
		}
		int index = findMinFreqIndex(list, n + 1);
		return list.get(n).getFrequency() < list.get(index).getFrequency() ? n : index;
	}
	/**
	 * Returns index of minimum value in double array
	 * @param arr double[]
	 * @return int - index
	 */
	public static int findMinDiffIndex(double[] arr) {
		double min = 1001.0; //start out with impossible high value
		int res = -1; //if no min return -1
		for (int i = 0; i < arr.length; i++) {
			if(arr[i] < min) {
				min = arr[i];
				res = i;
			}
		}
		return res;
	}
	/**
	 * Returns a String with the sequence of an equivalent codon of lower GC content with equal or greater frequency
	 * @param oldCodon String of original codon sequence
	 * @param codonKey ArrayList<AminoAcidKey>
	 * @return String - sequence of new codon
	 */
	public static String replaceCodon(String oldCodon, ArrayList<AminoAcidKey> codonKey) {
				String newCodon = oldCodon; //initalize newCodon with oldCodon incase better codon not found
				double oldCodonFreq = -1;
				double oldCodonGC = -1;
				int aaIndex = -1;
				//finds codon in key
				for (int i = 0; i < codonKey.size(); i++) {
					ArrayList<Codon> tempCodons = codonKey.get(i).getCodons();
					for (int j = 0; j < tempCodons.size(); j++) {
						if(tempCodons.get(j).getDnaSeq().equals(oldCodon)) {
							aaIndex = i;
							oldCodonFreq = tempCodons.get(j).getFrequency();
							oldCodonGC = tempCodons.get(j).getGc();
						} 
					}
				}
				
				//finds lower gc codon with max frequency in AminoAcidKey
				if (aaIndex != -1) {
					AminoAcidKey tempKey = codonKey.get(aaIndex);
					int maxIndex = -1;
					newCodon = oldCodon;
					for (int i = 0; i < tempKey.getCodons().size(); i++) {
						maxIndex = findMaxFreqIndex(tempKey.getCodons(), i);
						if (tempKey.getCodons().get(maxIndex).getGc() <= oldCodonGC) {
							if (tempKey.getCodons().get(maxIndex).getFrequency() >= oldCodonFreq) {
								newCodon = tempKey.getCodons().get(maxIndex).getDnaSeq();
								return newCodon;
							} 
							
						} 
					}
				}
				return newCodon;
		
	}
	/**
	 * Returns the difference in frequency between the next possible switch that lowers GC
	 * @param oldCodon String with codon sequence
	 * @param codonKey ArrayList<AminoAcidKey>
	 * @return double - difference of frequencies of codon and next best codon
	 */
	public static double getFreqDifferences(String oldCodon, ArrayList<AminoAcidKey> codonKey) {
		double diff = 1001.0; //set to 1001 since frequency is per thousand so this value is impossible
		double oldCodonFreq = -1;
		double oldCodonGC = -1;
		int aaIndex = -1;
		
		//finds codon in key
		for (int i = 0; i < codonKey.size(); i++) {
			ArrayList<Codon> tempCodons = codonKey.get(i).getCodons();
			for (int j = 0; j < tempCodons.size(); j++) {
				if(tempCodons.get(j).getDnaSeq().equals(oldCodon)) {
					aaIndex = i;
					oldCodonFreq = tempCodons.get(j).getFrequency();
					oldCodonGC = tempCodons.get(j).getGc();
					
				} 
			}
		}
		
		//finds lower gc codon with minimum difference in frequency
		if (aaIndex != -1) {
			AminoAcidKey tempKey = codonKey.get(aaIndex);
			for (int i = 0; i < tempKey.getCodons().size(); i++) {
				if (tempKey.getCodons().get(i).getGc() < oldCodonGC) { // set to < instead of <= like in first pass since we only want to change if it will lower the GC content
					double temp = tempKey.getCodons().get(i).getFrequency();
					if (oldCodonFreq - temp < diff ) {
						diff = oldCodonFreq - temp;
						return diff;	
					}
				}
			}
		}
		return diff;
	}
	/**
	 * Returns a String with the sequence of an equivalent codon of lower GC content with no regard to frequency. Used to switch codons in the second pass
	 * @param oldCodon String of original codon sequence
	 * @param codonKey ArrayList<AminoAcidKey>
	 * @return String - sequence of new codon 
	 */
	public static String replaceCodonRepeat(String oldCodon, ArrayList<AminoAcidKey> codonKey) {
		String newCodon = "";
//		double oldCodonFreq = -1;
		double oldCodonGC = -1;
		int aaIndex = -1;
		
		//finds codon in key
		for (int i = 0; i < codonKey.size(); i++) {
			ArrayList<Codon> tempCodons = codonKey.get(i).getCodons();
			for (int j = 0; j < tempCodons.size(); j++) {
				if(tempCodons.get(j).getDnaSeq().equals(oldCodon)) {
					aaIndex = i;
					oldCodonGC = tempCodons.get(j).getGc();
				} 
			}
		}
		
		//finds lower gc codon with minimum difference in frequency
		if (aaIndex != -1) {
			AminoAcidKey tempKey = codonKey.get(aaIndex);
			newCodon = oldCodon;
			for (int i = 0; i < tempKey.getCodons().size(); i++) {
				if (tempKey.getCodons().get(i).getGc() < oldCodonGC) { //less than not less than and equal to
						newCodon = tempKey.getCodons().get(i).getDnaSeq();
						return newCodon;
				} 
			}
		}
		return newCodon;
	}
	/**
	 * Finds Codon in key given string containing dna sequence of codon
	 * @param dna String of DNA sequence
	 * @param codonKey ArrayList<AminoAcidKey>
	 * @return Codon
	 */
	public static Codon findCodon(String dna, ArrayList<AminoAcidKey> codonKey) {
		Codon res = null;
		for (int i = 0; i < codonKey.size(); i++) {
			ArrayList<Codon> tempCodons = codonKey.get(i).getCodons();
			for (int j = 0; j < tempCodons.size(); j++) {
				if(tempCodons.get(j).getDnaSeq().equals(dna)) {
					res = tempCodons.get(j);
					break;
				} 
			}
		}
		if(res == null) {
			System.out.println("Error: Codon \""+dna+"\" not found");
		}
		return res;
	}
}