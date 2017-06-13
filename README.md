# IARC bioinformatics nextflow pipelines

## Nextflow pipelines list

### Raw NGS data processing
* [alignment-nf](https://github.com/IARCbioinfo/alignment-nf)
* [BQSR-nf](https://github.com/IARCbioinfo/BQSR-nf)
* [abra-nf](https://github.com/IARCbioinfo/abra-nf)
* [RNAseq-nf](https://github.com/IARCbioinfo/RNAseq-nf)

### QC
* [conpair-nf](https://github.com/IARCbioinfo/conpair-nf)
* [damage-estimator-nf](https://github.com/IARCbioinfo/damage-estimator-nf)
* [bamsurgeon-nf](https://github.com/IARCbioinfo/bamsurgeon-nf)

### Variant calling
* [platypus-nf](https://github.com/IARCbioinfo/platypus-nf)
* [GATK-Alignment-nf](https://github.com/IARCbioinfo/GATK-Alignment-nf)
* [CODEX-nf](https://github.com/IARCbioinfo/CODEX-nf)
* [needlestack](https://github.com/IARCbioinfo/needlestack)
* [mutect-nf](https://github.com/IARCbioinfo/mutect-nf)
* [strelka-nf](https://github.com/IARCbioinfo/strelka-nf)

### Other
* [addreplacerg-nf](https://github.com/IARCbioinfo/addreplacerg-nf)

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

