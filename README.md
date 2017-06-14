# IARC bioinformatics nextflow pipelines

## Nextflow pipelines list

### Raw NGS data processing
| Name      | Description     |
|-----------|-----------------| 
| [alignment-nf](https://github.com/IARCbioinfo/alignment-nf)    | Performs BAM realignment or fastq alignment, with/without local indel realignment and base quality score recalibration |
| [BQSR-nf](https://github.com/IARCbioinfo/BQSR-nf)   | ...... |
| [abra-nf](https://github.com/IARCbioinfo/abra-nf)   | ...... |
| [RNAseq-nf](https://github.com/IARCbioinfo/RNAseq-nf)   | RNAseq mapping, quality control, and reads counting nextflow pipeline - See also [RNAseq_analysis_scripts](https://github.com/IARCbioinfo/RNAseq_analysis_scripts) for post-processing  |
| [GATK-Alignment-nf](https://github.com/IARCbioinfo/GATK-Alignment-nf)   | Performs bwa alignment and pre-processing (realignment and recalibration) following GATK best practices (less performant than [alignment-nf](https://github.com/IARCbioinfo/alignment-nf) )

 |

### QC
| Name      | Description     |
|-----------|-----------------| 
| [conpair-nf](https://github.com/IARCbioinfo/conpair-nf)   | ...... |
| [damage-estimator-nf](https://github.com/IARCbioinfo/damage-estimator-nf)   | Nextflow pipeline to run "Damage Estimator" |
| [bamsurgeon-nf](https://github.com/IARCbioinfo/bamsurgeon-nf)   | ...... |

### Variant calling
| Name      | Description     |
|-----------|-----------------| 
| [platypus-nf](https://github.com/IARCbioinfo/platypus-nf)   | Platypus germline variant calling |
| [GATK-Alignment-nf](https://github.com/IARCbioinfo/GATK-Alignment-nf)   | Runs variant calling in GVCF mode on bam files, joint genotyping and variant recalibration (SNPs and indels) following GATK best practices |
| [CODEX-nf](https://github.com/IARCbioinfo/CODEX-nf)   | ...... |
| [needlestack](https://github.com/IARCbioinfo/needlestack)   | A multi-sample somatic variant caller |
| [mutect-nf](https://github.com/IARCbioinfo/mutect-nf)   | Mutect pipeline on tumor-matched normal bam folder |
| [strelka-nf](https://github.com/IARCbioinfo/strelka-nf)   | Strelka pipeline with Nextflow |

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

