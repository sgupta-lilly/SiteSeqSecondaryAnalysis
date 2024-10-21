import sys
import pandas as pd
from pybedtools import BedTool
import os
import glob
import time
import subprocess
import shutil
####Usage
#secondary analysis v1 
#1. create result_secondary_pipeline
#2. pathtosamplesheet = path to samplesheet
#3. pathtobam = path to bam files. These should be in user_data/sorted_bam_bai/ 
#3a. check name of bam file _REP1.mLb.clN.sorted.bam; if different change line 83 and 122
#4. path to results_offtarget where the named *OT_sum.tsv per sample. The expected name of files = "sample" column from samplesheet (used to run nf-lly-siteseq) + "_cameron_macs2_combined_OT_sum.tsv". If the suffix is different then 
#change line 35
#/lrlhps/users/c195933/mambaforge/bin/python3 secondary_analysis.py pathtosamplesheet pathtobam pathresults_offtarget

def copy_file_to_current(src_file_path):
    """
    Copies a file from a given source path to the current directory.
    
    Args:
    src_file_path (str): Path to the source file
    
    Returns:
    None
    """
    # Get the base name of the file (i.e., just the file name without directories)
    file_name = os.path.basename(src_file_path)
    
    # Define the destination as the current directory
    dest_file_path = os.path.join(os.getcwd(), file_name)
    
    # Copy the file
    shutil.copy(src_file_path, dest_file_path)
    
    print(f"File '{file_name}' copied to the current directory.")



#####PART 1
#input samplesheet
#df_metafile = pd.read_csv("/lrlhps/genomics/prod/lgm/dna_editing/BN24-11550_APOC3_Top8_SiteSeq/BN24-11550_SITESeq_APOC3Top8_nf_08-22-24.csv", sep = ",")
df_metafile = pd.read_csv(sys.argv[1], sep = ",")
#bam_path = "/lrlhps/genomics/prod/lgm/dna_editing/BN24-11550_APOC3_Top8_SiteSeq/user_data/sorted_bam_bai/"
bam_path = sys.argv[2]
results_offtarget = sys.argv[3]

results_secondary_analysis = results_offtarget+"result_secondary_pipeline/"

#sample_list = (df_metafile['sample'].tolist())
sample_list = (df_metafile.iloc[:, 0].tolist())
for eachsample in sample_list:
    print(eachsample)
# Define the input file from macs and cameron_OT_Sum.csv
#input_file = sys.argv[1]
    input_file = results_offtarget+eachsample+"_cameron_macs2_combined_OT_sum.tsv"
    copy_file_to_current(input_file)

# Load the data into a DataFrame and select only the necessary columns
    df = pd.read_csv(input_file, sep='\t', usecols=['Chr', 'OT_Start', 'OT_End'])

# Remove duplicate rows
    df = df.drop_duplicates()
    df = df.sort_values(by=['Chr', 'OT_Start', 'OT_End'])

# Save the cleaned data to a BED file
#bed_file = 'output.bed'
    bed_file = eachsample+"_output.bed"
    df.to_csv(bed_file, sep='\t', header=False, index=False)

# Use pybedtools to remove overlaps
    bed = BedTool(bed_file)
    bed_non_overlapping = bed.sort()

# Save the non-overlapping BED file
    final_out = eachsample+"_non_overlapping_output.bed"
    bed_non_overlapping.saveas(final_out)
    remove_cmd = (f"rm {bed_file}")
    os.system(remove_cmd)
    print('BED file with non-overlapping regions created: non_overlapping_output.bed')
#### file for merge 

# Load the data into a DataFrame and select only the necessary columns
    df_merge = pd.read_csv(input_file, sep='\t', usecols=[
    'OT_ID', 'Chr', 'DSB_pos', 'OT_Start', 'OT_End', 'Strand',
    'OT_Mismatch', 'OT_length', 'OT_pam_pos', 'PAM_Mismatch', 
    'RNA_bulge', 'DNA_bulge', 'Cleavage_pos', 'IGV_region'])    
    df_merge['ON_Target'] = df_merge.apply(lambda row: '1' if (row['OT_Mismatch'] <= 6 and row['PAM_Mismatch'] == 0 and row['RNA_bulge'] == 0 and row['DNA_bulge'] == 0) else '0', axis=1)
    if 'Cleavage_pos' not in df_merge.columns:
        df_merge['Cleavage_pos'] = ''  # Or provide a default value

# Group by 'Chr', 'start', and 'end', then concatenate values
    grouped_df = df_merge.groupby([
    'OT_ID', 'Chr', 'OT_Start', 'OT_End', 'Strand', 
    'OT_Mismatch', 'OT_length', 'OT_pam_pos', 'PAM_Mismatch', 
    'RNA_bulge', 'DNA_bulge'], as_index=False).agg({
    'DSB_pos': lambda x: ','.join(map(str, x)),
    'Cleavage_pos': lambda x: ','.join(map(str, x))
    })




# run multiBamsummary
    bamfile = bam_path+eachsample+"_REP1.mLb.clN.sorted.bam"
    count_npz = eachsample+"_counts.npz"
    count_tab = eachsample+"_sites_read_counts.tab"
    multiBamSummary_cmd = (
    f"/lrlhps/users/c195933/mambaforge/bin/multiBamSummary BED-file "
    f" --BED {final_out} --ignoreDuplicates --minMappingQuality 30 "
    f" --minFragmentLength 60 --bamfiles {bamfile} "
    f" -out {count_npz} --outRawCounts {count_tab} "
    )
    # Execute the multiBamSummary command
    print("Running multiBamSummary with deepTools...")
    os.system(multiBamSummary_cmd)


# Wait if the coverage file is not ready
    while not os.path.exists(count_tab):
        print(f"Waiting for {count_tab} to be available...")
        time.sleep(20)


## combine coverage
    df_second = pd.read_csv(count_tab, sep='\t', header=None, names=['Chr', 'OT_Start', 'OT_End', 'Cov'], skiprows=1)

    df_processed = grouped_df

    df_processed['OT_Start'] = df_processed['OT_Start'].astype(int)
    df_processed['OT_End'] = df_processed['OT_End'].astype(int)
    df_second['OT_Start'] = df_second['OT_Start'].astype(int)
    df_second['OT_End'] = df_second['OT_End'].astype(int)

# Merge the DataFrames on 'Chr', 'OT_Start', and 'OT_End'
    merged_df = pd.merge(df_processed, df_second, on=['Chr', 'OT_Start', 'OT_End'], how='inner')
    merged_output_file = eachsample + '_merged_output.bed'
    merged_df.to_csv(merged_output_file, sep='\t', header=True, index=False)




# Run multiBamSummary on each chunk
updated_list = [bam_path+element + "_REP1.mLb.clN.sorted.bam" for element in sample_list]
space_separated_bams = ' '.join(updated_list)

#####PART 2
# create coverage across all the sites
# Get a list of all 'sample*.merged.bed' files in the current folder
bed_files = [f for f in os.listdir('.') if f.endswith('_non_overlapping_output.bed')]

# Initialize an empty dataframe to store concatenated data
df_concat = pd.DataFrame()

# Iterate through each bed file and concatenate them
for bed_file in bed_files:
    print(f"Processing {bed_file}")
    df = pd.read_csv(bed_file, sep='\t', header=None, names=['Chr', 'Start', 'End'])
    df_concat = pd.concat([df_concat, df], ignore_index=True)

# Sort the concatenated dataframe by 'Chr', 'Start', and 'End' columns
df_concat_sorted = df_concat.sort_values(by=['Chr', 'Start', 'End'])

# Drop duplicates
df_final = df_concat_sorted.drop_duplicates()

# Number of lines per file
lines_per_file = 6000

# Split the dataframe into chunks of 1000 rows
num_chunks = (len(df_final) // lines_per_file) + 1  # Calculate how many files needed
for i in range(num_chunks):
    start_idx = i * lines_per_file
    end_idx = (i + 1) * lines_per_file
    chunk = df_final.iloc[start_idx:end_idx]
    
    # Define the output filename for each chunk
    output_file = f'all_site_part_{i+1}.bed'
    
    # Save each chunk to a separate file
    chunk.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"Saved {output_file} with {len(chunk)} lines.")


    #bam_file = bam_path+'*sorted.bam'  # Replace with actual BAM file path
    output_npz = f'all_site_part_{i+1}_counts.npz'
    output_raw_counts = f'all_site_part_{i+1}_read_counts.tab'

    multiBamSummary_cmd = (
    f"/lrlhps/users/c195933/mambaforge/bin/multiBamSummary BED-file "
    f" --BED {results_secondary_analysis+output_file} --ignoreDuplicates --minMappingQuality 30 "
    f" --minFragmentLength 60 --bamfiles {space_separated_bams} "
    f" -out {results_secondary_analysis+output_npz} --outRawCounts {results_secondary_analysis+output_raw_counts} "
    )



   # print(multiBamSummary_cmd)
   # print("Running multiBamSummary with deepTools...")
   # os.system(multiBamSummary_cmd)

    # Create a qsub script for the job
    qsub_script = f"""
    #!/bin/bash
    #PBS -N multiBamSummary_job_{i+1}
    #PBS -l walltime=24:00:00
    #PBS -l nodes=1:ppn=4
    #PBS -j oe

        cd $PBS_O_WORKDIR
    echo "Running multiBamSummary for {output_file}"
    {multiBamSummary_cmd}
    """


    # Save the qsub script to a file
    qsub_script_filename = f"multiBamSummary_job_{i+1}.sh"
    with open(qsub_script_filename, 'w') as f:
        f.write(qsub_script)

    # Submit the job using qsub
    print(f"Submitting job for {output_file}")
    os.system(f"qsub {qsub_script_filename}")
