import java.util.ArrayList;
/**
 * <h1>Sequence</h1>
 * 
 * @author Jessica Mattick <a target = "_blank" href = "https://github.com/jmatty16">https://github.com/jmatty16</a>
 * @version 1.0
 * @since 2018-12-09
 * @see <a href = "AminoAcidKey.html">AminoAcidKey</a>
 * @see <a href = "Codon.html">Codon</a>
 */
public class Sequence {
	private String dnaSeq;
	private ArrayList<AminoAcidKey> codonKey = new ArrayList<AminoAcidKey>();
	private ArrayList<Codon> seqArr = new ArrayList<Codon>();
	private double gcContent;
	private double lowestFreq;
	private double avgFreq;
	
	//Constructor
	Sequence(String dnaSeq, ArrayList<AminoAcidKey> codonKey, ArrayList<Codon> seqArr, double gcContent, double lowestFreq, double avgFreq) {
		this.dnaSeq = dnaSeq;
		this.codonKey = codonKey;
		this.seqArr = seqArr;
		this.gcContent = gcContent;
		this.lowestFreq = lowestFreq;
		this.avgFreq = avgFreq;
	}
	
	Sequence(String dnaSeq, ArrayList<AminoAcidKey> codonKey) {
		this.dnaSeq = dnaSeq;
		this.codonKey = codonKey;
	}
	
	Sequence() {
		
	}

	//Set/Get

	public String getDnaSeq() {
		return dnaSeq;
	}

	public void setDnaSeq(String dnaSeq) {
		this.dnaSeq = dnaSeq;
	}

	public ArrayList<AminoAcidKey> getCodonKey() {
		return codonKey;
	}

	public void setCodonKey(ArrayList<AminoAcidKey> codonKey) {
		this.codonKey = codonKey;
	}

	public ArrayList<Codon> getSeqArr() {
		return seqArr;
	}

	public void setSeqArr(ArrayList<Codon> seqArr) {
		this.seqArr = seqArr;
	}
	public double getGcContent() {
		return gcContent;
	}

	public void setGcContent(double gcContent) {
		this.gcContent = gcContent;
	}

	public double getLowestFreq() {
		return lowestFreq;
	}

	public void setLowestFreq(double lowestFreq) {
		this.lowestFreq = lowestFreq;
	}

	public double getAvgFreq() {
		return avgFreq;
	}

	public void setAvgFreq(double avgFreq) {
		this.avgFreq = avgFreq;
	}
	
	//Methods
	/**
	 * Updates Sequence data based on dnaSeq
	 * @param dnaSeq String of DNA sequence
	 */
	public void updateSeq(String dnaSeq) {
		this.seqArr = CodonOptimization.formatSeqArr(dnaSeq, this.codonKey);
		this.gcContent = CodonOptimization.calcGC(this.dnaSeq);
		this.lowestFreq = seqArr.get(CodonOptimization.findMinFreqIndex(seqArr, 0)).getFrequency();
		this.avgFreq = CodonOptimization.calcAvgFreq(seqArr);
		
	}
	
	/**
	 * Updates Sequence data based on seqArr
	 * @param seqArr ArrayList<Codon>
	 */
	public void updateSeq(ArrayList<Codon> seqArr) {
		this.dnaSeq = CodonOptimization.getDNASeq(seqArr);
		this.gcContent = CodonOptimization.calcGC(this.dnaSeq);
		this.lowestFreq = seqArr.get(CodonOptimization.findMinFreqIndex(seqArr, 0)).getFrequency();
		this.avgFreq = CodonOptimization.calcAvgFreq(seqArr);
		
	}
	/**
	 * Returns an array of the frequencies in a sequence
	 * @return double[]
	 */
	public double[] returnFreqArr() {
		double[] res = new double[this.seqArr.size()];
		for (int i = 0; i < this.seqArr.size(); i++) {
			res[i] = this.seqArr.get(i).getFrequency();
		}
		return res;
	}
	/**
	 * Prints sequence, GC content, Lowest Frequency, and Average Frequency 
	 */
	public void printSeqInfo() {
		System.out.println("\t"+ this.dnaSeq);
		System.out.println("\tGC content: \t"+ this.gcContent);
		System.out.println("\tLowest Freq: \t" + this.lowestFreq);
		System.out.println("\tAverage Freq: \t" + this.avgFreq);
	}
	
	
	
	
	

}
