# SiteSeqSecondaryAnalysis

## This repo is to be run after SiteSeq nf-lly-siteseq workflow

## Usage
### secondary analysis v1 
1. create result_secondary_pipeline
2. pathtosamplesheet = path to samplesheet
3. pathtobam = path to bam files. These should be in user_data/sorted_bam_bai/<br />
4. pathto2secondary analysis
3a. check name of bam files _REP1.mLb.clN.sorted.bam; if different change line 83 and 122<br />

### command <br />
**/lrlhps/users/c195933/mambaforge/bin/python3 secondary_analysis.py pathtosamplesheet pathtobam pathto2secondary analysis **

<br />

4. path to  pathto2secondary analysis where the named *OT_sum.tsv per sample. <br />
5. The expected name of files = "sample" column from samplesheet (used to run nf-lly-siteseq) + "_cameron_macs2_combined_OT_sum.tsv". If the suffix is different then 
change line 35


<br />
* Specific example

<br />
** python3 secondary_analysis.py BN24-11550_SITESeq_APOC3Top8_nf_08-22-24.csv ../user_data/sorted_bam_bai/ ../results_offtarget/ ../pathto2secondaryanalysis/
<br />
<br />
<br />
After the secondary_analysis.py 

* run Rmarkdown
* guidename as in the samplesheet 
* results_folder = user_data
* samplesheet = samplesheet
* result_secondary_pipeline = folder create to run secondary_analysis.py
<br />
Rscript -e "library(rmarkdown);rmarkdown::render('SiteSeqGuideReport.rmd',  params = list(guidename = '17rev', results_folder='/user_data/', samplesheet='BN24-11550_SITESeq_APOC3Top8_nf_08-22-24.csv', result_secondary_pipeline='/result_secondary_pipeline/'), output_file = '17rev.html')"
<br />

Finally annotated with gencode.v46.annotation.gtf from https://www.gencodegenes.org/human/release_21.html and cancer_gene_census.csv
<br />
python3 annotation.py 17rev.OT_rank_out cancer_gene_census.csv

17rev.OT_rank_out is a table with all offtarget sites for 17rev from the rmarkdown. Annotation.py annotates the sites with gencode.v46.annotation.gtf that should be in the directory and cosmic census (licenced file) attached to this repo


