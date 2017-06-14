# IARC bioinformatics nextflow pipelines

## Nextflow pipelines list

### Raw NGS data processing
| Name      | Description     |
|-----------|-----------------| 
| [alignment-nf](https://github.com/IARCbioinfo/alignment-nf)    | Performs BAM realignment or fastq alignment, with/without local indel realignment and base quality score recalibration |
| [BQSR-nf](https://github.com/IARCbioinfo/BQSR-nf)   | Performs base quality score recalibration of bam files using GATK |
| [abra-nf](https://github.com/IARCbioinfo/abra-nf)   | Runs ABRA (Assembly Based ReAligner) |
| [RNAseq-nf](https://github.com/IARCbioinfo/RNAseq-nf)   | Performs RNAseq mapping, quality control, and reads counting - See also [RNAseq_analysis_scripts](https://github.com/IARCbioinfo/RNAseq_analysis_scripts) for post-processing  |
| [GATK-Alignment-nf](https://github.com/IARCbioinfo/GATK-Alignment-nf)   | Performs bwa alignment and pre-processing (realignment and recalibration) following GATK best practices (less performant than [alignment-nf](https://github.com/IARCbioinfo/alignment-nf) )

 |

### QC
| Name      | Description     |
|-----------|-----------------| 
| [conpair-nf](https://github.com/IARCbioinfo/conpair-nf)   | Runs conpair (concordance and contamination estimator) |
| [damage-estimator-nf](https://github.com/IARCbioinfo/damage-estimator-nf)   | Runs "Damage Estimator" |
| [bamsurgeon-nf](https://github.com/IARCbioinfo/bamsurgeon-nf)   | Runs bamsurgeon with step of variant simulation |

### Variant calling
| Name      | Description     |
|-----------|-----------------| 
| [platypus-nf](https://github.com/IARCbioinfo/platypus-nf)   | Runs Platypus (germline variant caller) |
| [GATK-Calling-GVCF-nf](https://github.com/IARCbioinfo/GATK-Calling-GVCF-nf)   | Runs variant calling in GVCF mode on bam files, joint genotyping and variant recalibration (SNPs and indels) following GATK best practices - still in development!!|
| [CODEX-nf](https://github.com/IARCbioinfo/CODEX-nf)   | Performs copy number variant calling from whole exome sequencing data using CODEX |
| [needlestack](https://github.com/IARCbioinfo/needlestack)   | Performs multi-sample somatic variant calling |
| [mutect-nf](https://github.com/IARCbioinfo/mutect-nf)   | Runs Mutect on tumor-matched normal bam pairs |
| [strelka-nf](https://github.com/IARCbioinfo/strelka-nf)   | Runs Strelka |

### Other
| Name      | Description     |
|-----------|-----------------| 
| [addreplacerg-nf](https://github.com/IARCbioinfo/addreplacerg-nf)   | ...... |

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

