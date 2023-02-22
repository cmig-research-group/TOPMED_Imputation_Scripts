This collection of scripts and functions is used for pre and post processing of [Michigan](https://imputation.biodatacatalyst.nhlbi.nih.gov/) or [TOPMED](https://imputation.biodatacatalyst.nhlbi.nih.gov/#!) imputation servers. 



# Dependancies
Use of these scripts requires the following to be installed:

* [plink](https://www.cog-genomics.org/plink/)
* [plink2](https://www.cog-genomics.org/plink/2.0/)
* [hstlib](https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856#file-install-samtools-bcftools-and-htslib-md)


You can use the links above to install these libraries, alternatively should also be able to install these using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#installing-conda-on-a-system-that-has-other-python-installations-or-packages):
```
conda install -c bioconda plink plink2 hstlib
```


# Pre-Imputation

Before imputation `vcf_convert_gzip.py` takes a plink file (merged across chromosomes), seperates them into chromosomes, converts to vcf and gzips for upload to the servers. Usage is shown below:

```
python vcf_convert_gzip.py -bed_file /path/to/merged_bed_file_prefix -out_file /path/to_output/vcf/files_to_upload_chr
```

# Post-Imputation

Imputation result from server comes as vcf.gz files without rsids. After all steps are complete it should procduce plink files (.bed, .bim, .fam) at a best guess threshold of prob>0.9, with rsids instead of chr:bp:a1:a2 names that come out from the imputation servers. **Note:** these scripts will also remove duplicate IDs - this will remove triallelic variants. 

`prep_RSID_files.py` downloads snp RSIDs (hg38 for TOPMED and hg19 for HRC) to a given directory for renaming SNP ids. This process only needs to be done once and does not need to be performed for each new imputation result.

`plink_conversion.py` takes the unziped files from the imputation server and converts to plink, removes duplicate ids (triallelic variants) and renames to availible RSIDs that are processed from `prep_RSID_files.py`

Example of how these scripts are used is shown below.


```
snp_dir='/path/to/hg38'
post_imput_dir='/path/to/post_impute_dir'
# Password emailed from Imputation Server
password='PASSWORD' 

mkdir -p $snp_dir
# Download RSID files for renaming SNPS (hg38 for TOPMED)
python prep_RSID_files.py -snp_dir $snp_dir -build hg38
# Unzip files (may take time)
for i in ${post_imput_dir}/*zip; do unzip -P $password "$i" -d ${post_imput_dir} & done
# Convert to plink, rename snps and remove duplicates
python plink_conversion.py -post_imput_dir $post_imput_dir -snp_dir $snp_dir
```