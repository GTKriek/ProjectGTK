# This python script is used to parse the gff3 file containing the complete genome of Harmonia axyridis, this allowed to 
# search for all Pannier binding sites in the entire genome and print the gene ids in a output file 258genes(wchr).txt. 
from BCBio import GFF
import re

def parse_regions(txt_file):
    """Parse the regions from the .txt file."""
    regions = []
    with open(txt_file, 'r') as f:
        for line in f:
            match = re.match(r'Find \d+: (\d+)-(\d+)', line)
            if match:
                start = int(match.group(1))
                end = int(match.group(2))
                regions.append((start, end))
    return regions

def search_gff3(gff3_file, regions):
    """Search for the regions in the .gff3 file using BCBio.GFF and match them to chromosomes."""
    matches = []
    with open(gff3_file, 'r') as f:
        # Parse the HARgen.gff3
        for rec in GFF.parse(f):
            chrom = rec.id
            for feature in rec.features:
                gff_start = int(feature.location.start) + 1  
                gff_end = int(feature.location.end)

                # Check region overlaps with GFF entry
                for start, end in regions:
                    if gff_start <= end and gff_end >= start:
                        matches.append((chrom, f"{start}-{end}"))
                        break
    return matches

def write_output(output_file, matches):
    """Write the matched regions and chromosomes to an output file."""
    with open(output_file, 'w') as f:
        for chrom, region in matches:
            f.write(f"{chrom}\t{region}\n")

def main():
    txt_file = 'findings(4dist)reduced.txt'  
    gff3_file = 'HARgen.gff3' 
    output_file = '258genes(wchr).txt'  

    # Parse the regions from text file
    regions = parse_regions(txt_file)

    # Search for regions in the HARgen.gff3 file
    matches = search_gff3(gff3_file, regions)

    write_output(output_file, matches)

    print(f"Output written to {output_file}")

if __name__ == '__main__':
    main()
