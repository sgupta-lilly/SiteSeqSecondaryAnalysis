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


<br />
/lrlhps/users/c195933/mambaforge/bin/python3 secondary_analysis.py /lrlhps/genomics/prod/lgm/dna_editing/BN24-11550_APOC3_Top8_SiteSeq/BN24-11550_SITESeq_APOC3Top8_nf_08-22-24.csv /lrlhps/genomics/prod/lgm/dna_editing/BN24-11550_APOC3_Top8_SiteSeq/user_data/sorted_bam_bai/ /lrlhps/genomics/prod/lgm/dna_editing/BN24-11550_APOC3_Top8_SiteSeq/results_DSB_Cameron_Macs2_combined_08-25-2024/

<br />
4. path to results_offtarget where the named *OT_sum.tsv per sample. <br />
5. The expected name of files = "sample" column from samplesheet (used to run nf-lly-siteseq) + "_cameron_macs2_combined_OT_sum.tsv". If the suffix is different then 
change line 35

