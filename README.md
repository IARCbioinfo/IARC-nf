# IARC bioinformatics nextflow pipelines

## Nextflow pipelines list

### Raw NGS data processing
| Name      | Description     |	Tools used	|
|-----------|-----------------|-----------------|
| [alignment-nf](https://github.com/IARCbioinfo/alignment-nf)    | Performs BAM realignment or fastq alignment, with/without local indel realignment and base quality score recalibration |[bwa](https://github.com/lh3/bwa), [samblaster](https://github.com/GregoryFaust/samblaster), [sambamba](https://github.com/lomereiter/sambamba) & in option: [samtools](http://samtools.sourceforge.net/), [AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval), [GATK](www.broadinstitute.org/gatk/download), [k8 javascript execution shell](https://sourceforge.net/projects/bio-bwa/files/bwakit/), [bwa-postalt.js](https://github.com/lh3/bwa/tree/master/bwakit) |
| [BQSR-nf](https://github.com/IARCbioinfo/BQSR-nf)   | Performs base quality score recalibration of bam files using GATK |[samtools](http://samtools.sourceforge.net/), [samblaster](https://github.com/GregoryFaust/samblaster), [sambamba](https://github.com/lomereiter/sambamba), [GATK](www.broadinstitute.org/gatk/download)|
| [abra-nf](https://github.com/IARCbioinfo/abra-nf)   | Runs ABRA (Assembly Based ReAligner) |[ABRA](https://github.com/mozack/abra), [bedtools](http://bedtools.readthedocs.io/en/latest/), [bwa](http://bio-bwa.sourceforge.net), [sambamba](http://lomereiter.github.io/sambamba/), [samtools](http://www.htslib.org/) |
| [RNAseq-nf](https://github.com/IARCbioinfo/RNAseq-nf)   | Performs RNAseq mapping, quality control, and reads counting - See also [RNAseq_analysis_scripts](https://github.com/IARCbioinfo/RNAseq_analysis_scripts) for post-processing  |[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt), [cutadapt](http://cutadapt.readthedocs.io/en/stable/installation.html), Python version > 2.7, [trim_galore](https://github.com/FelixKrueger/TrimGalore), [RESeQC](http://rseqc.sourceforge.net/), [multiQC](http://multiqc.info/docs/), [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), [htseq](http://www-huber.embl.de/HTSeq/doc/install.html#install) & in option: [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml), [GATK](www.broadinstitute.org/gatk/download), [samtools](http://samtools.sourceforge.net/)|
| [GATK-Alignment-nf](https://github.com/IARCbioinfo/GATK-Alignment-nf)   | Performs bwa alignment and pre-processing (realignment and recalibration) following GATK best practices (less performant than [alignment-nf](https://github.com/IARCbioinfo/alignment-nf) ) |[bwa](https://github.com/lh3/bwa), picard, [GATK](www.broadinstitute.org/gatk/download)|

### QC
| Name      | Description     |	Tools used	|
|-----------|-----------------|-----------------|
| [conpair-nf](https://github.com/IARCbioinfo/conpair-nf)   | Runs conpair (concordance and contamination estimator) |[conpair](https://github.com/nygenome/Conpair), [Python 2.7](www.python.org), [numpy 1.7.0 or higher](www.numpy.org), [scipy 0.14.0 or higher](www.scipy.org), [GATK 2.3 or higher](www.broadinstitute.org/gatk/download)|
| [damage-estimator-nf](https://github.com/IARCbioinfo/damage-estimator-nf)   | Runs "Damage Estimator" |[Damage Estimator](https://github.com/Ettwiller/Damage-estimator), [samtools](http://samtools.sourceforge.net/), [R](https://www.r-project.org) with GGPLOT2 package|
| [bamsurgeon-nf](https://github.com/IARCbioinfo/bamsurgeon-nf)   | Runs bamsurgeon with step of variant simulation |[Python 2.7](www.python.org), [bamsurgeon](http://github.com/adamewing/bamsurgeon/), [R software](https://www.r-project.org/) (tested with R version 3.2.3)|

### Variant calling
| Name      | Description     |	Tools used	|
|-----------|-----------------|-----------------|
| [platypus-nf](https://github.com/IARCbioinfo/platypus-nf)   | Runs Platypus (germline variant caller) |[Platypus](https://github.com/andyrimmer/Platypus)|
| [GATK-Calling-GVCF-nf](https://github.com/IARCbioinfo/GATK-Calling-GVCF-nf)   | Runs variant calling in GVCF mode on bam files, joint genotyping and variant recalibration (SNPs and indels) following GATK best practices - still in development!!|[GATK](www.broadinstitute.org/gatk/download)|
| [CODEX-nf](https://github.com/IARCbioinfo/CODEX-nf)   | Performs copy number variant calling from whole exome sequencing data using CODEX |[R](https://www.r-project.org) with package Codex, Rscript |
| [needlestack](https://github.com/IARCbioinfo/needlestack)   | Performs multi-sample somatic variant calling |[perl](https://www.perl.org),  [bedtools](http://bedtools.readthedocs.org/en/latest/), [samtools](http://samtools.sourceforge.net/) and Rscript |
| [mutect-nf](https://github.com/IARCbioinfo/mutect-nf)   | Runs Mutect on tumor-matched normal bam pairs |[Mutect](https://github.com/broadinstitute/mutect) and its dependencies (Java 1.7 and Maven 3.0+), [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)|
| [strelka-nf](https://github.com/IARCbioinfo/strelka-nf)   | Runs Strelka |[Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/home/strelka-workflow-installation)|

### Other
| Name      | Description     |	Tools used	|
|-----------|-----------------|-----------------|
| [addreplacerg-nf](https://github.com/IARCbioinfo/addreplacerg-nf)   | Adds and replaces read group tags in BAM files |[samtools](http://samtools.sourceforge.net/)|

## Nextflow 

### Installation : 

1. Install [java](https://java.com/download/) JRE if you don't already have it (7 or higher).

2. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```
  
### Configuration file

### Pipelines updates :

You can update the nextflow sofware and the pipeline itself simply using:
```bash
nextflow -self-update
nextflow pull iarcbioinfo/pipeline_name
```

You can also automatically update the pipeline when you run it by adding the option `-latest` in the `nextflow run` command. Doing so you will always run the latest version from Github.

### Display help :

```bash
nextflow run iarcbioinfo/pipeline_name --help
```

## Docker

Install [docker](https://www.docker.com).
	
This is very system specific (but quite easy in most cases), follow  [docker documentation](https://docs.docker.com/installation/). Also follow the optional configuration step called `Create a Docker group` in their documentation.

To run nextflow pipeline with Docker, simply add the `-with-docker` option.

