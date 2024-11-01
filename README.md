# ProjectGTK
repository for supplementary material and code 

**Skripts**
**1. #1(bindingsites).ipynb**
Identifies specific nucleotide patterns in a genomic sequence. The script scans for pattern_1 and pattern_2 within a defined distance and calculates the nucleotide distance between them.

Inputs: genome_Seq: genomic sequence to analyze.
        pattern_1: A specific nucleotide sequence ("GATG").
        pattern_2: A list of potential matching nucleotide patterns ("AGATAG", "TGATAG", "TGATAA", "AGATAA").

Outputs: Prints the occurrences of pattern_1 and pattern_2 in genome_Seq, with each position and distance of the patterns.
         Displays the nucleotide distance between matched patterns if any are found.

**2. #2(boyer-moore).py**
Searches for occurrences of pattern_1 and pattern_2 within a genome sequence using the Boyer-Moore algorithm. It efficiently identifies pattern_1 sites and locates pattern_2 patterns at a specified gap distance from the pattern_1 occurrences.

Inputs:gff3_file: HARgen.gff3
       fasta_file: HARgenome.fasta
       pattern_1: The initial nucleotide pattern to locate ("GATG").
       pattern_2: List of secondary nucleotide patterns ("AGATAG", "TGATAG", "TGATAA", "AGATAA").
       exact_gap: 4 nucleotides separating pattern_1 and pattern_2.

Outputs: Results are saved to findings1.txt and include occurrences of pattern_1 and pattern_2, their positions, and distances between them.

**3. #3(geneids).py**
This script parses the GFF3 file containing the genome of Harmonia axyridis to search for Pannier binding site regions and outputs the gene IDs with corresponding chromosome locations in 258genes(wchr).txt.

Inputs:txt_file: findings(4dist)reduced.txt listing genomic regions identified as potential Pannier binding sites.
       gff3_file: HARgen.gff3 is the genome annotation file 

Outputs: Writes results to 258genes(wchr).txt, listing chromosome IDs alongside genomic regions where Pannier binding sites overlap gene features.

**4. #4(geneidspromoter).py**
This script identifies all Pannier binding sites that are within 2000 base pairs (bp) upstream of the start of each gene among the previously identified 258 genes in the Harmonia axyridis genome. The output is saved in finalgeneids(110).txt.

Inputs:findings_file: findings(4dist)reduced.txt, containing genomic coordinates of Pannier binding sites.
       gff_file: The GFF3 annotation file with gene information for Harmonia axyridis (HARgen.gff3).

Outputs:Writes results to finalgeneids(110).txt, listing gene IDs with binding site coordinates in the upstream promoter region.

**Supplementary material**
1. Figure 2-5
   Figure 2: Crystal structures of Pannier, 3VD6 and superimposed structures
   Figure 3: multiple sequence alignment of 3VD6 and Pannier
   Figure 4: Protein-Protein interface of Pannier
   Figure 5: Elytral patterns at different diffusion rates
