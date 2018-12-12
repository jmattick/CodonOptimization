# CodonOptimization

Optimizes DNA sequences for cloning by lowering the GC content through codon optimization, while avoiding reducing the codon frequencies. 

DNA sequences go through either one or two passes depending on the input sequence and the parameters. The first pass optimization only changes codons if there is an alternate codon that has a lower GC% and an equal or higher frequency. If it is impossible to lower the GC content of the sequence without lowering the frequency, then a second pass optimization will occur. The second pass optimization determines the reduction in frequency of the next best codon at each index of the sequence and replaces the codons starting at the smallest frequency difference until the sequence meets the GC parameter. 



See docs: https://jmatty16.github.io/CodonOptimization/doc/
