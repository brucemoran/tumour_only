# Tumour_Only DNAseq Analysis
## Call SNV somatic sample data
### How to Setup
#### Dependencies:
[NextFlow](https://www.nextflow.io/index.html#GetStarted), [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
#### Reference Generation:
Run the included `download-references.nf`
```
nextflow run brucemoran/somatic_n-of-1/download-references.nf \
  -profile        [str]     singularity,refs
  --assembly      [str]     GRCh37 or GRCh38
  --exomeTag      [str]     naming for exome/panel outputs when supplied; tag is then used in somatic_n-of-1 and batch_somatic main.nf pipelines to select relevant exome reference data

  one of:
  --exomeBedURL   [str]     URL to exome bed file for intervals
  --exomeBedFile  [str]     locally downloaded exome bed file for intervals

General Optional Arguments:

  --exomeAssembly [str]     what assembly is the exome BED based on? (default: GRCH37)
  --cosmicUser    [str]     COSMIC login credentials, user email
  --cosmicPass    [str]     COSMIC login credentials, password
```
### Tumour_Only Pipeline
#### About the pipeline:
This pipeline was developed to analyse and report on clinical cancer research data. To this end we have tried to make a useful and clinician-readable output. This has been achieved largely on the back of [PCGR](https://github.com/sigven/pcgr)/[CPSR](https://github.com/sigven/cpsr) which provide really excellent HTML reports and annotation from multiple clinically relevant sources.

Variant calling uses MuTect2.

#### To run the pipeline:
```
nextflow run brucemoran/tumour_only -profile standard,singularity
```

#### Arguments:
```
nextflow run brucemoran/tumour_only --help
```
#### Formats of --sampleCsv, --sampleCat
To input specific fastq files per samples:
```
sampleID,meta,read1,read2
soma1,primary_tumour,/full/path/to/soma1.R1.fastq.gz,/full/path/to/soma1.R2.fastq.gz
soma2,metastasis,/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz
```
To input directories containing multiple fastqs named {sampleID}*{ext}
```
sampleID,meta,dir,ext
germ1,whole_blood,/full/path/to/dir/with/germ1/fqs,R1.fastq.gz;R2.fastq.gz
soma1,primary_tumour,/full/path/to/dir/with/soma1/fqs,R1.fastq.gz;R2.fastq.gz
soma2,metastasis,/full/path/to/dir/with/soma2/fqs,R1.fastq.gz;R2.fastq.gz
```
The `meta` column is used for reporting in PCGR, CPSR where `sampleID` may include clinical/personal data, for example.

Headers of `sample.csv` file must match above exactly, and you should have only one germline/normal sample per run.

Both germline CPSR results via HaplotypeCaller, and somatic PCGR results via MuTect2 calling, are returned for the sample.

####Singularity Hub and containers
Unfortunately shub where the containers were housed is gone, so you may need to build the containers. We are working on a quay.io repo. Please pull our github Singularity repo (brucemoran/Singularity) and use our recipe_builder.sh script:
```
git clone https://github.com/brucemoran/Singularity
cd Singularity
sudo sh ./utilities/recipe_builder.sh <recipe_name> <$NXF_SINGULARITY_CACHE> shub
```
This builds a correctly named container into the $NXF_SINGULARITY_CACHE so the pipeline will run.
