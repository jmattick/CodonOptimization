/**
 * <h1>Codon</h1>
 * 
 * @author Jessica Mattick <a target = "_blank" href = "https://github.com/jmatty16">https://github.com/jmatty16</a>
 * @version 1.0
 * @since 2018-12-09
 *
 */
public class Codon {
	private String dnaSeq;
	private String aminoAcid;
	private double gc;
	private double frequency;
	
	//Constructor
	public Codon(String dnaSeq, String aminoAcid, double gc, double frequency) {
		this.dnaSeq = dnaSeq;
		this.aminoAcid = aminoAcid;
		this.gc = gc;
		this.frequency = frequency;	
	}
	
	//Set/Get
	public String getDnaSeq() {
		return dnaSeq;
	}

	public void setDnaSeq(String dnaSeq) {
		this.dnaSeq = dnaSeq;
	}

	public String getAminoAcid() {
		return aminoAcid;
	}

	public void setAminoAcid(String aminoAcid) {
		this.aminoAcid = aminoAcid;
	}

	public double getGc() {
		return gc;
	}

	public void setGc(double gc) {
		this.gc = gc;
	}

	public double getFrequency() {
		return frequency;
	}

	public void setFrequency(double frequency) {
		this.frequency = frequency;
	}
	
	//Methods
	public void display() {
		System.out.println("Codon: " + dnaSeq + "\nAminoAcid: " + aminoAcid + "\nGC%: " + gc + "\nFrequency: " + frequency);
	}
}
