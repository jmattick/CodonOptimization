import java.util.ArrayList;
/**
 * <h1>AminoAcidKey</h1>
 * 
 * @author Jessica Mattick <a target = "_blank" href = "https://github.com/jmatty16">https://github.com/jmatty16</a>
 * @version 1.0
 * @since 2018-12-09
 * @see <a href = "Codon.html">Codon</a>
 */
public class AminoAcidKey {
	private String name;
	private ArrayList<Codon> codons = new ArrayList<Codon>();
	
	//Constructor
	public AminoAcidKey(String name) {
		this.name = name;
		
	}
	
	//Set/Get
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public ArrayList<Codon> getCodons() {
		return codons;
	}
	public void setCodons(ArrayList<Codon> codons) {
		this.codons = codons;
	}
	
	//Methods
	public void addCodon(Codon newCodon) {
		codons.add(newCodon);
	}
	public void display() {
		System.out.println(name);
		for(int i = 0; i < codons.size(); i++) {
			System.out.println("\t" + codons.get(i).getDnaSeq() + "\t" + codons.get(i).getFrequency() + "\t" + codons.get(i).getGc());
		}
	}
	

}
