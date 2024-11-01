# This python script is used to parse the gff3 file containing the complete genome of Harmonia axyridis, this allowed to 
# search for all Pannier binding sites 2000bp upstream of the start region in all 258 identified genes and print out the
# gene ids in a output file finalgeneids(110).txt. 
from Bio import SeqIO
from BCBio import GFF

def parse_findings(file_path):

    findings = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Find"):
                # start and end 
                start, end = map(int, line.split(":")[1].strip().split("-"))
                findings.append((start, end))
    return findings

def parse_gff(gff_file):
    genes = []
    with open(gff_file, 'r') as f:
        for rec in GFF.parse(f):
            for feature in rec.features:
                if feature.type == "gene":
                    genes.append(feature)
    return genes

def find_genes_with_upstream_binding(genes, findings, upstream_distance=2000):
    # if 2000bp upstream 
    result_genes = []
    for gene in genes:
        gene_start = gene.location.start
        gene_end = gene.location.end
        gene_strand = gene.strand

        # Calculate the upstream region based on the gene's strand
        if gene_strand == 1:  # +strand
            upstream_start = max(0, gene_start - upstream_distance)
            upstream_end = gene_start
        else:  # -strand
            upstream_start = gene_end
            upstream_end = gene_end + upstream_distance

        # Check if any finding overlaps with the upstream region
        for start, end in findings:
            if upstream_start <= start <= upstream_end or upstream_start <= end <= upstream_end:
                result_genes.append((gene, (start, end)))
                break

    return result_genes

def write_output(result_genes, output_file):
    with open(output_file, 'w') as f:
        for gene, binding_site in result_genes:
            gene_id = gene.qualifiers.get("ID", ["Unknown"])[0]
            binding_start, binding_end = binding_site
            f.write(f"{gene_id}\t{binding_start}-{binding_end}\n")

def main():
    # File paths
    findings_file = "findings(4dist)reduced.txt"
    gff_file = "HARgen.gff3"  
    output_file = "finalgeneids(110).txt"


    findings = parse_findings(findings_file)
    genes = parse_gff(gff_file)

    result_genes = find_genes_with_upstream_binding(genes, findings)

    write_output(result_genes, output_file)

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main() 