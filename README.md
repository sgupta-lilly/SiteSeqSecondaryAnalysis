# SiteSeqSecondaryAnalysis
SiteSeqSecondaryAnalysis
This repo is to be run after SiteSeq nf-lly-siteseq workflow

Usage
secondary analysis v1 
1. create result_secondary_pipeline
2. pathtosamplesheet = path to samplesheet
3. pathtobam = path to bam files. These should be in user_data/sorted_bam_bai/<br />
3a. check name of bam files _REP1.mLb.clN.sorted.bam; if different change line 83 and 122<br />
Usage<br />
**/lrlhps/users/c195933/mambaforge/bin/python3 secondary_analysis.py pathtosamplesheet pathtobam pathresults_offtarget**


5. path to results_offtarget where the named *OT_sum.tsv per sample. The expected name of files = "sample" column from samplesheet (used to run nf-lly-siteseq) + "_cameron_macs2_combined_OT_sum.tsv". If the suffix is different then 
change line 35

