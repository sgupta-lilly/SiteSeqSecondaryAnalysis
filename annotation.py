import pandas as pd
import re
import sys
import os

def load_gtf(gtf_file):
    """Load the GTF file into a DataFrame."""
    gtf_columns = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    return pd.read_csv(gtf_file, sep='\t', comment='#', names=gtf_columns)

def load_info(info_file):
    """Load the info file into a DataFrame."""
    return pd.read_csv(info_file, sep='\t', comment='#')

def load_cosmic_data(cosmic_file):
    """Load the COSMIC gene information file into a DataFrame."""
    cosmic_columns = [
        "Gene Symbol", "Name", "Entrez GeneId", "Genome Location", "Tier", 
        "Hallmark", "Chr Band", "Somatic", "Germline", "Tumour Types(Somatic)", 
        "Tumour Types(Germline)", "Cancer Syndrome", "Tissue Type", 
        "Molecular Genetics", "Role in Cancer", "Mutation Types", 
        "Translocation Partner", "Other Germline Mut", "Other Syndrome", "Synonyms"
    ]
    return pd.read_csv(cosmic_file, sep=',', header=None, names=cosmic_columns)

def extract_gene_name(attribute):
    """Extract gene name from GTF attribute."""
    match = re.search(r'gene_name "([^"]+)"', attribute)
    return match.group(1) if match else None

def annotate_bed(row, exons, genes):
    """Annotate genomic coordinates based on overlap with GTF regions."""
    chrom = row['Chr']
    start = row['OT_Start']
    end = row['OT_End']
    
    # Check for exon overlap
    exon_overlap = exons[(exons['chr'] == chrom) & (exons['start'] <= end) & (exons['end'] >= start)]
    if not exon_overlap.empty:
        gene_name = exon_overlap['gene_name'].values[0]
        return gene_name, "exon", 1
    
    # Check for gene overlap (introns)
    gene_overlap = genes[(genes['chr'] == chrom) & (genes['start'] <= end) & (genes['end'] >= start)]
    if not gene_overlap.empty:
        gene_name = gene_overlap['gene_name'].values[0]
        return gene_name, "intron", 2
    
    return None, "intergenic", 3

def find_cosmic_info(gene, cosmic_data):
    """Find the Tier and Role in Cancer based on gene symbol."""
    cosmic_entry = cosmic_data[cosmic_data['Gene Symbol'] == gene]
    if not cosmic_entry.empty:
        tier = cosmic_entry['Tier'].values[0]
        role_in_cancer = cosmic_entry['Role in Cancer'].values[0]
        return tier, role_in_cancer
    return None, None

def main(gtf_file, info_file, cosmic_file):
    """Main function to process files and annotate genomic coordinates."""
    if not os.path.exists(gtf_file):
        print(f"Error: GTF file '{gtf_file}' does not exist.")
        sys.exit(1)
    
    if not os.path.exists(info_file):
        print(f"Error: Info file '{info_file}' does not exist.")
        sys.exit(1)
    
    if not os.path.exists(cosmic_file):
        print(f"Error: COSMIC gene information file '{cosmic_file}' does not exist.")
        sys.exit(1)

    gtf = load_gtf(gtf_file)

    # Add gene name column to the GTF DataFrame
    gtf['gene_name'] = gtf['attribute'].apply(extract_gene_name)

    # Filter GTF for relevant features
    exons = gtf[gtf['feature'] == 'exon']
    genes = gtf[gtf['feature'] == 'gene']

    # Load info file
    info = load_info(info_file)
    info['OT_Start'] = pd.to_numeric(info['OT_Start'], errors='coerce')
    info['OT_End'] = pd.to_numeric(info['OT_End'], errors='coerce')

    # Apply annotation function
    annotations = info.apply(lambda row: annotate_bed(row, exons, genes), axis=1)
    annotations_df = pd.DataFrame(annotations.tolist(), columns=['Gene', 'Annotation', 'Rank_Annotation'])
    
    # Load COSMIC data
    cosmic_data = load_cosmic_data(cosmic_file)

    # Find COSMIC Tier and Role in Cancer based on gene symbol
    annotations_df[['Tier', 'Role in Cancer']] = annotations_df['Gene'].apply(
        lambda gene: find_cosmic_info(gene, cosmic_data)
    ).apply(pd.Series)

    # Concatenate the original DataFrame with the new annotations
    final_result = pd.concat([info, annotations_df], axis=1)

    # Save annotated info.txt
    output_file = info_file + '_annotated.txt'
    final_result.to_csv(output_file, sep='\t', index=False)
    print(f"Annotated file saved as {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python annotate_genomes.py <path_to_info.txt> <path_to_cosmic_file.csv>")
        sys.exit(1)

    gtf_file = "gencode.v46.annotation.gtf"  # Update the path if necessary
    info_file = sys.argv[1]
    cosmic_file = sys.argv[2]  # New argument for COSMIC file

    main(gtf_file, info_file, cosmic_file)
