# This python script makes use of the Boyer-Moore algorithm to find the Pannier binding sites in the entire
# genome of Harmonia axyridis.
from Bio import SeqIO

def preprocess_bad_char(pattern):
    bad_char_shift = {}
    length = len(pattern)
    for i in range(length - 1):
        bad_char_shift[pattern[i]] = length - i - 1
    return bad_char_shift

# Boyer-Moore 
def boyer_moore_search(text, pattern):
    bad_char_shift = preprocess_bad_char(pattern)
    n = len(text)
    m = len(pattern)
    i = 0
    matches = []
    
    while i <= n - m:
        j = m - 1 
        while j >= 0 and pattern[j] == text[i + j]:
            j -= 1
        if j < 0:
            matches.append(i)
            i += (m if i + m < n else 1)  
        else:
            i += max(1, bad_char_shift.get(text[i + j], m))
    
    return matches

def find_patterns_bm(genome_seq, pattern_1, pattern_2, output_file, exact_gap=4):
    distances = []
    
    pattern_1_matches = boyer_moore_search(genome_seq, pattern_1)
    
    for i in pattern_1_matches:
        j = i + len(pattern_1) + exact_gap  
        
        if j + 6 <= len(genome_seq):  
            current_subseq = genome_seq[j:j+6]  
            if current_subseq in pattern_2:
                
                distance = j - (i + len(pattern_1)) 
                distances.append((i, j, distance, current_subseq))
    
    with open(output_file, 'w') as f:
        if not distances:
            f.write("No occurrences of the patterns were found.\n")
        else:
            for idx, (start_p1, start_p2, distance, found_pattern) in enumerate(distances):
                f.write(f"Find {idx + 1}:\n")
                f.write(f"Pattern_1 ({pattern_1}) found at index {start_p1}\n")
                f.write(f"Pattern_2 ({found_pattern}) found at index {start_p2}\n")
                f.write(f"Number of nucleotides between the patterns: {distance}\n\n")

def read_genome_from_gff3(gff3_file, fasta_file):
    genome_seq = ""
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_seq += str(record.seq)
    
    sequences = []
    with open(gff3_file, 'r') as gff3:
        for line in gff3:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                if len(columns) > 8 and columns[2] == "CDS": 
                    start = int(columns[3])
                    end = int(columns[4])
                    seq = genome_seq[start-1:end]
                    sequences.append(seq)
    
    return "".join(sequences)

gff3_file = "HARgen.gff3"  
fasta_file = "HARgenome.fasta" 
genome_seq = read_genome_from_gff3(gff3_file, fasta_file)

pattern_1 = "GATG"  
pattern_2 = ["AGATAG", "TGATAG", "TGATAA", "AGATAA"]  
output_file = "findings1.txt"  

find_patterns_bm(genome_seq, pattern_1, pattern_2, output_file, exact_gap=4)
