# IARC bioinformatics nextflow pipelines (updated on 05/11/2020)

This page lists all the pipelines developed at IARC (mostly nextflow pipelines which are suffixed with -nf) and explains how to use them (at the bottom of the page)

## IARC pipelines list

### Raw NGS data processing
| Name      |Latest version|Maintained | Description     |	Tools used	|
|-----------|--------------|-----------|-----------------|----------------------|
| [alignment-nf](https://github.com/IARCbioinfo/alignment-nf)    |v1.2 - July 2020|Yes| Performs BAM realignment or fastq alignment, with/without local indel realignment and base quality score recalibration |[bwa](https://github.com/lh3/bwa), [samblaster](https://github.com/GregoryFaust/samblaster), [sambamba](https://github.com/lomereiter/sambamba), [samtools](http://samtools.sourceforge.net/), [AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval), [GATK](www.broadinstitute.org/gatk/download), [k8 javascript execution shell](https://sourceforge.net/projects/bio-bwa/files/bwakit/), [bwa-postalt.js](https://github.com/lh3/bwa/tree/master/bwakit) |
| [gatk4-DataPreProcessing-nf](https://github.com/IARCbioinfo/gatk4-DataPreProcessing-nf)   |Nov 2018|?| Performs bwa alignment and pre-processing (mark duplicates and recalibration) following GATK4 best practices - compatible with hg38 |[bwa](https://github.com/lh3/bwa), [picard](https://broadinstitute.github.io/picard/), [GATK4](https://software.broadinstitute.org/gatk/download/), [sambamba](https://github.com/lomereiter/sambamba), [qualimap](http://qualimap.bioinfo.cipf.es/)|
| [GATK-Alignment-nf](https://github.com/IARCbioinfo/GATK-Alignment-nf)   |June 2017|No| Performs bwa alignment and pre-processing (realignment and recalibration) following first version of GATK best practices (less performant than [alignment-nf](https://github.com/IARCbioinfo/alignment-nf) ) |[bwa](https://github.com/lh3/bwa), [picard](https://broadinstitute.github.io/picard/), [GATK](www.broadinstitute.org/gatk/download)|
| [BQSR-nf](https://github.com/IARCbioinfo/BQSR-nf)   |v1.1 - Apr 2020|Yes| Performs base quality score recalibration of bam files using GATK |[samtools](http://samtools.sourceforge.net/), [samblaster](https://github.com/GregoryFaust/samblaster), [sambamba](https://github.com/lomereiter/sambamba), [GATK](www.broadinstitute.org/gatk/download)|
| [abra-nf](https://github.com/IARCbioinfo/abra-nf)   |v3.0 - Apr 2020|Yes| Runs ABRA (Assembly Based ReAligner) |[ABRA](https://github.com/mozack/abra), [bedtools](http://bedtools.readthedocs.io/en/latest/), [bwa](http://bio-bwa.sourceforge.net), [sambamba](http://lomereiter.github.io/sambamba/), [samtools](http://www.htslib.org/) |
| [PostAlignment-nf](https://github.com/IARCbioinfo/PostAlignment-nf)   |Aug 2018|?| Perform post alignment on bam files | [samtools](http://samtools.sourceforge.net/), [sambamba](https://github.com/lomereiter/sambamba), [bwa-postalt.js](https://github.com/lh3/bwa/tree/master/bwakit)|
|*************** |||
| [marathon-wgs](https://github.com/IARCbioinfo/marathon-wgs)   |June 2018|?| Studies intratumor heterogeneity with Canopy|[bwa](https://github.com/lh3/bwa), [platypus](https://github.com/andyrimmer/Platypus), [strelka2](https://github.com/Illumina/strelka), [vt](https://github.com/atks/vt), [annovar](http://annovar.openbioinformatics.org/en/latest/), [R](https://www.r-project.org), [Falcon](https://cran.r-project.org/web/packages/falcon/index.html), [Canopy](https://github.com/yuchaojiang/Canopy)|
| [ITH-nf](https://github.com/IARCbioinfo/ITH-nf)   |Sept 2018|?| Perform intra-tumoral heterogeneity (ITH) analysis |[Strelka2](https://github.com/Illumina/strelka) , [Platypus](https://www.well.ox.ac.uk/platypus), [Bcftools](https://samtools.github.io/bcftools/bcftools.html), [Tabix](http://www.htslib.org/doc/tabix.html), [Falcon](https://omictools.com/falcon-3-tool), [Canopy](https://github.com/yuchaojiang/Canopy)|
| [ITH_pipeline](https://github.com/IARCbioinfo/ITH_pipeline)   |August 2019|?| Study intra-tumoral heterogeneity (ITH) through subclonality reconstruction |[HATCHet](https://github.com/raphael-group/hatchet) , [DeCiFer](https://github.com/raphael-group/decifer), [ClonEvol](https://github.com/hdng/clonevol)|

### RNA Seq
| Name      |Latest version|Maintained | Description     |	Tools used	|
|-----------|--------------|-----------|-----------------|----------------------|
| [RNAseq-nf](https://github.com/IARCbioinfo/RNAseq-nf)   |v2.3 - May 2020|Yes| Performs RNAseq mapping, quality control, and reads counting - See also [RNAseq_analysis_scripts](https://github.com/IARCbioinfo/RNAseq_analysis_scripts) for post-processing  |[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt), [RESeQC](http://rseqc.sourceforge.net/), [multiQC](http://multiqc.info/docs/), [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), [htseq](http://www-huber.embl.de/HTSeq/doc/install.html#install), [cutadapt](http://cutadapt.readthedocs.io/en/stable/installation.html), Python version > 2.7, [trim_galore](https://github.com/FelixKrueger/TrimGalore), [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml), [GATK](www.broadinstitute.org/gatk/download), [samtools](http://samtools.sourceforge.net/)|
| [RNAseq-transcript-nf](https://github.com/IARCbioinfo/RNAseq-transcript-nf)   |v2.2 - June 2020|Yes| Performs transcript identification and quantification from a series of BAM files |[StringTie](https://github.com/gpertea/stringtie)|
| [RNAseq-fusion-nf](https://github.com/IARCbioinfo/RNAseq-fusion-nf)   |v1.1 - Aug 2020|Yes| Perform fusion-genes discovery from RNAseq data|[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki)|
| [quantiseq-nf](https://github.com/IARCbioinfo/quantiseq-nf)   |v1.1 - July 2020|Yes| Quantify immune cell content from RNA-seq data|[quanTIseq](https://icbi.i-med.ac.at/software/quantiseq/doc/) |

### QC
| Name      |Latest version|Maintained | Description     |	Tools used	|
|-----------|--------------|-----------|-----------------|----------------------|
| [NGSCheckMate](https://github.com/IARCbioinfo/NGSCheckMate)   |v1.1 - Aug 2020|Yes| Runs NGSCheckMate on BAM files to identify data files from a same indidual (i.e. check N/T pairs) |[NGSCheckMate](https://github.com/parklab/NGSCheckMate)|
| [conpair-nf](https://github.com/IARCbioinfo/conpair-nf)   |June 2018|?| Runs conpair (concordance and contamination estimator) |[conpair](https://github.com/nygenome/Conpair), [Python 2.7](www.python.org), [numpy 1.7.0 or higher](www.numpy.org), [scipy 0.14.0 or higher](www.scipy.org), [GATK 2.3 or higher](www.broadinstitute.org/gatk/download)|
| [damage-estimator-nf](https://github.com/IARCbioinfo/damage-estimator-nf)   |June 2017|?| Runs "Damage Estimator" |[Damage Estimator](https://github.com/Ettwiller/Damage-estimator), [samtools](http://samtools.sourceforge.net/), [R](https://www.r-project.org) with GGPLOT2 package|
| [QC3](https://github.com/IARCbioinfo/QC3)   |May 2016|No| Runs QC on DNA seq data (raw data, aligned data and variant calls - forked from [slzhao](https://github.com/slzhao/QC3) |[samtools](http://samtools.sourceforge.net/)|
| [fastqc-nf](https://github.com/IARCbioinfo/fastqc-nf)   |v1.1 - July 2020|Yes| Runs fastqc and multiqc on DNA seq data (fastq data) |[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Multiqc](https://multiqc.info/)|
| [qualimap-nf](https://github.com/IARCbioinfo/qualimap-nf)   |v1.1 - Nov 2019|Yes| Performs quality control on bam files (WES, WGS and target alignment data) |[samtools](http://samtools.sourceforge.net/), [Qualimap](http://qualimap.bioinfo.cipf.es/), [Multiqc](https://multiqc.info/)|
| [mpileup-nf](https://github.com/IARCbioinfo/mpileup-nf)   |Jan 2018|?| Computes bam coverage with samtools mpileup (bed parallelization) |[samtools](http://samtools.sourceforge.net/),[annovar](http://annovar.openbioinformatics.org/en/latest/)|
| [bamsurgeon-nf](https://github.com/IARCbioinfo/bamsurgeon-nf)   |Mar 2019|?| Runs bamsurgeon (tool to add mutations to bam files) with step of variant simulation |[Python 2.7](www.python.org), [bamsurgeon](http://github.com/adamewing/bamsurgeon/), [R software](https://www.r-project.org/) (tested with R version 3.2.3)|

### Variant calling
| Name      |Latest version|Maintained | Description     |	Tools used	|
|-----------|--------------|-----------|-----------------|----------------------|
| [needlestack](https://github.com/IARCbioinfo/needlestack)   |v1.1 - May 2019|Yes| Performs multi-sample somatic variant calling |[perl](https://www.perl.org),  [bedtools](http://bedtools.readthedocs.org/en/latest/), [samtools](http://samtools.sourceforge.net/) and [R software](https://www.r-project.org/) |
| [target-seq](https://github.com/IARCbioinfo/target-seq)   |Aug 2019|?| Whole pipeline to perform multi-sample somatic variant calling using Needlestack on targeted sequencing data|[abra2](https://github.com/IARCbioinfo/abra-nf),[QC3](https://github.com/slzhao/QC3) ,[needlestack](https://github.com/IARCbioinfo/needlestack), [annovar](http://annovar.openbioinformatics.org/en/latest/) and [R software](https://www.r-project.org/) |
| [strelka2-nf](https://github.com/IARCbioinfo/strelka2-nf)   |v1.2 - Oct 2020|Yes| Runs Strelka 2 (germline and somatic variant caller)|[Strelka2](https://github.com/Illumina/strelka)|
| [strelka-nf](https://github.com/IARCbioinfo/strelka-nf)   |Jun 2017|No| Runs Strelka (germline and somatic variant caller)|[Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/home/strelka-workflow-installation)|
| [mutect-nf](https://github.com/IARCbioinfo/mutect-nf)   |v2.2b - Oct 2020|Yes| Runs Mutect on tumor-matched normal bam pairs |[Mutect](https://github.com/broadinstitute/mutect) and its dependencies (Java 1.7 and Maven 3.0+), [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)|
| [gatk4-HaplotypeCaller-nf](https://github.com/IARCbioinfo/gatk4-HaplotypeCaller-nf)   |Dec 2019|?| Runs variant calling in GVCF mode on bam files following GATK best practices|[GATK](https://software.broadinstitute.org/gatk/download/)|
| [gatk4-GenotypeGVCFs-nf](https://github.com/IARCbioinfo/gatk4-GenotypeGVCFs-nf)   |Apr 2019|?| Runs joint genotyping on gvcf files following GATK best practices|[GATK](https://software.broadinstitute.org/gatk/download/)|
| [GVCF_pipeline-nf](https://github.com/IARCbioinfo/GVCF_pipeline-nf)   |Nov 2016|?| Performs bam realignment and recalibration + variant calling in GVCF mode following GATK best practices|[bwa](https://github.com/lh3/bwa), [samblaster](https://github.com/GregoryFaust/samblaster), [sambamba](https://github.com/lomereiter/sambamba), [GATK](https://software.broadinstitute.org/gatk/download/)|
| [platypus-nf](https://github.com/IARCbioinfo/platypus-nf)   |v1.0 - Apr 2018|?| Runs Platypus (germline variant caller) |[Platypus](https://github.com/andyrimmer/Platypus)|
| [TCGA_platypus-nf](https://github.com/IARCbioinfo/TCGA_platypus-nf)   |Aug 2018|?| Converts TCGA Platypus vcf in format for annotation with annovar |[vt](https://github.com/atks/vt),[VCFTools](https://github.com/vcftools/vcftools)|
| [vcf_normalization-nf](https://github.com/IARCbioinfo/vcf_normalization-nf)   |v1.1 - May 2020|Yes| Decomposes and normalizes variant calls (vcf files) |[bcftools](https://github.com/samtools/bcftools),[samtools/htslib](http://www.htslib.org/)|
| [TCGA_germline-nf](https://github.com/IARCbioinfo/TCGA_germline-nf) |May 2017|?| Extract germline variants from TCGA data for annotation with annovar (vcf files) |[R software](https://www.r-project.org/)|
| [gama_annot-nf](https://github.com/IARCbioinfo/gama_annot-nf) |Aug 2020|Yes| Filter and annotate batch of vcf files (annovar + strand + context) |[annovar](http://annovar.openbioinformatics.org/en/latest/), [R](https://www.r-project.org)|
| [table_annovar-nf](https://github.com/IARCbioinfo/table_annovar-nf) |v1.1 - May 2020|Yes| Annotate variants with annovar (vcf files) |[annovar](http://annovar.openbioinformatics.org/en/latest/)|
|*************** |||
| [MutSpec](https://github.com/IARCbioinfo/mutspec)   |v2.0 - May 2017|?| Suite of tools for analyzing and interpreting mutational signatures |[annovar](http://annovar.openbioinformatics.org/en/latest/)|
|*************** |||
| [CODEX-nf](https://github.com/IARCbioinfo/CODEX-nf)   |Mar 2017|?| Performs copy number variant calling from whole exome sequencing data using CODEX |[R](https://www.r-project.org) with package Codex, Rscript |
| [facets-nf](https://github.com/IARCbioinfo/facets-nf)   |v2.0 - Oct 2020|Yes| Performs fraction and copy number estimate from tumor/normal sequencing data using facets |[facets](https://github.com/mskcc/facets) , [R](https://www.r-project.org) |
| [svaba-nf](https://github.com/IARCbioinfo/svaba-nf)   |v1.0 - August 2020|Yes| Performs structural variant calling using SvABA |[SvABA](https://github.com/walaj/svaba) , [R](https://www.r-project.org) |


### Other tools/pipelines
| Name      |Latest version|Maintained | Description     |	Tools used	|
|-----------|--------------|-----------|-----------------|----------------------|
| [template-nf](https://github.com/IARCbioinfo/template-nf)   |May 2020|Yes | Empty template for nextflow pipelines |NA|
| [data_test](https://github.com/IARCbioinfo/data_test)   |Aug 2020|Yes| Small data files to test IARC nextflow pipelines |NA|
||||
| [scanMyWorkDir](https://github.com/IARCbioinfo/scanMyWorkDir)   |May 2018|?| Non-destructive and informative scan of a nextflow work folder |NA|
| [addreplacerg-nf](https://github.com/IARCbioinfo/addreplacerg-nf)   |Jan 2017|?| Adds and replaces read group tags in BAM files |[samtools](http://samtools.sourceforge.net/)|
| [bametrics-nf](https://github.com/IARCbioinfo/bametrics-nf)   |Mar 2017|?| Computes average metrics from reads that overlap a given set of positions |NA|
| [Gviz_multiAlignments](https://github.com/IARCbioinfo/Gviz_multiAlignments)   |Aug 2017|?| Generates multiple BAM alignments views using Gviz bioconductor package|[Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html)||
| [nf_coverage_demo](https://github.com/IARCbioinfo/nf_coverage_demo)   |v2.3 - July 2020|Yes| Plots mean coverage over a series of BAM files |[bedtools](http://bedtools.readthedocs.io/en/latest/), [R software](https://www.r-project.org/)|
| [LiftOver-nf](https://github.com/IARCbioinfo/LiftOver-nf) |Nov 2017|?| Converts BED/VCF between hg19 and hg38 |[picard](https://broadinstitute.github.io/picard/)|
| [PVAmpliconFinder](https://github.com/IARCbioinfo/PVAmpliconFinder)   |Aug 2020|Yes| Identify and classify known and potentially new papilliomaviridae sequences from amplicon deep-sequencing with degenerated papillomavirus primers.|Python and Perl + FastQC, MultiQC, Trim Galore, VSEARCH, Blast, RaxML-EPA, PaPaRa, CAP3, KRONA)|
| [integration_analysis_scripts](https://github.com/IARCbioinfo/integration_analysis_scripts)   |Mar 2020|Yes| Performs unsupervised analyses (clustering) from transformed expression data (e.g., log fpkm) and methylation beta values |[R software](https://www.r-project.org/) with iClusterPlus, gplots and lattice R packages|
| [mpileup2readcounts](https://github.com/IARCbioinfo/mpileup2readcounts)|Apr 2018|?| Get the readcounts at a locus by piping samtools mpileup output - forked from [gatoravi](https://github.com/gatoravi/mpileup2readcounts) |[samtools](http://samtools.sourceforge.net/)|
| [Methylation_analysis_scripts](https://github.com/IARCbioinfo/Methylation_analysis_scripts)   |v1.0 - June 2020|Yes| Perform Illumina EPIC 850K array pre-processing and QC from idat files|[R software](https://www.r-project.org)| 
| [DRMetrics](https://github.com/IARCbioinfo/DRMetrics)   |Oct 2020|Yes| Evaluate the quality of projections obtained after using dimensionality reduction techniques|[R software](https://www.r-project.org/)|
| [acnviewer-singularity](https://github.com/IARCbioinfo/acnviewer-singularity)   |Jul 2019|?| Build a singularity image of aCNViewer (tool for visualization of absolute copy number and copy neutral variations) (|[Singularity](https://sylabs.io/singularity/)|
| [polysolver-singularity](https://github.com/IARCbioinfo/polysolver-singularity)   |Dec 2019|?| Build a singularity image of Polysolver (tool for HLA typing based on whole exome seq)|[Singularity](https://sylabs.io/singularity/)|

### NEW !!
| Name      |Latest version|Maintained | Description     |	Tools used	|
|-----------|--------------|-----------|-----------------|----------------------|
| [MinION_pipes](https://github.com/IARCbioinfo/MinION_pipes)   |Jan 2020|?| Analyze MinION sequencing data for the reconstruction of viral genomes |Guppy V3.1.5+, Porechop V0.2.4, Nanofilt V2.2.0, Filtlong V0.2.0, SPAdes V3.10.1, CAP3 02/10/15, BLAST V2.9.0+, MUSCLE V3.8.1551, Nanopolish V0.11.0, Minimap2 V2.15, Samtools version 1.9|
| [DraftPolisher](https://github.com/IARCbioinfo/DraftPolisher)   |Jan 2020|?| Fast polishing of draft sequences (draft genome assembly) |[MUSCLE](http://drive5.com/muscle/downloads.htm), [Python3](https://www.python.org/downloads/release/python-360/) |
| [Imputation-nf](https://github.com/IARCbioinfo/Imputation-nf)   |Aug 2020 - v1.0|Yes| Pipeline to perform dataset genotyping imputation |LiftOver, Plink, Admixture, Perl, Term::ReadKey, Becftools, Eagle, Minimac4 and samtools|

### Courses
| Name      | Description     |	Tools used	|
|-----------|-----------------|-----------------|
| [nextflow-course-2018](https://github.com/IARCbioinfo/nextflow-course-2018)   | Nextflow course |NA|
| [SBG-CGC_course2018](https://github.com/IARCbioinfo/SBG-CGC_course2018)   | Analyzing TCGA data in SBG-CGC |NA|
| [Medical-genomics-course2020](https://github.com/IARCbioinfo/medical_genomics_course)   | Medical Genomics course held at the INSA Lyon - Fall 2020|NA|

### Tricks
| Name      |Latest version|Maintained | Description     |	Tools used	|
|-----------|--------------|-----------|-----------------|----------------------|
| [BAM-tricks](https://github.com/IARCbioinfo/BAM-tricks)   ||| Tips and tricks for BAM files |[samtools](http://samtools.sourceforge.net/), freebayes, [bedtools](http://bedtools.readthedocs.io/en/latest/), biobambam2, [Picard](http://broadinstitute.github.io/picard/), [rbamtools](https://cran.r-project.org/web/packages/rbamtools/index.html)|
| [VCF-tricks](https://github.com/IARCbioinfo/VCF-tricks) ||| Tips and tricks for VCF files |[samtools](http://samtools.sourceforge.net/),[bcftools](https://github.com/samtools/bcftools), [vcflib](https://github.com/vcflib/vcflib), [vcftools](https://github.com/vcftools/vcftools), R scripts|
| [R-tricks](https://github.com/IARCbioinfo/R-tricks)||| Tips and tricks for R |NA|
| [EGA-tricks](https://github.com/IARCbioinfo/EGA-tricks)||| Tips and tricks to use the European Genome-Phenome Archive from the European Bioinformatics Institute |[EGA client](https://www.ebi.ac.uk/ega/sites/ebi.ac.uk.ega/files/documents/EGA_download_client_2.2.2.zip)|
| [GDC-tricks](https://github.com/IARCbioinfo/GDC-tricks)||| Tips and tricks to use the [GDC data portal](https://gdc-portal.nci.nih.gov/) |NA|
| [awesomeTCGA](https://github.com/IARCbioinfo/awesome-TCGA) ||| Curated list of resources to access TCGA data |NA|
| [LSF-Tricks](https://github.com/IARCbioinfo/LSF-tricks)   ||| Tips and tricks for LSF HPC scheduler |NA|

### Coming soon... (only dev branches yet)
| Name      | Description     |	Tools used	|
|-----------|-----------------|----------------------|
| [Nextflow_DSL2](https://github.com/IARCbioinfo/Nextflow_DSL2)   | Repository with modules for nextflow DSL2 |NA|
| [methylkey](https://github.com/IARCbioinfo/methylkey)   | Pipeline for 450k and 850k array analysis (bisulfite data analysis using Minfi, Methylumi, Comet, Bumphunter and DMRcate packages)|[R software](https://www.r-project.org/)|
| [variantflag](https://github.com/IARCbioinfo/variantflag)   | Merge and annotate variants from different callers ||
| [PRSice](https://github.com/IARCbioinfo/PRSice-nf)   | Pipeline to compute polygenic risk scores ||
| [EPIDRIVER2020](https://github.com/IARCbioinfo/EPIDRIVER2020)   | Scripts for EPIDRIVER Project ||

## Installation 

### Nextflow

1. Install [java](https://java.com/download/) JRE if you don't already have it (7 or higher).

2. Install [nextflow](http://www.nextflow.io/).

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	And move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```

### Docker

To avoid having to installing all dependencies each time you use a pipeline, you can instead install [docker](https://www.docker.com) and let nextflow dealing with it. Installing docker is system specific (but quite easy in most cases), follow  [docker documentation](https://docs.docker.com/installation/) (docker CE is sufficient). Also follow the post-installation step to manage Docker as a non-root user ([here](https://docs.docker.com/engine/installation/linux/linux-postinstall/) for Linux), otherwise you will need to change the `sudo` option in nextflow `docker` config scope as described in the nextflow documentation [here](https://www.nextflow.io/docs/latest/config.html#scope-docker).

To run nextflow pipeline with Docker, simply add the `-with-docker` option in the `nextflow run` command.

### Singularity

To avoid having to installing all dependencies each time you use a pipeline, you can also install [singularity](https://www.docker.com) and let nextflow dealing with it. 

See documentation [here](https://sylabs.io/guides/3.5/user-guide/).

In case you want to use the same singularity container - with the exactly same versions of pipeline and tools - on several data over time you may want to pull the container and archive it somewhere :
```bash
singularity pull shub://IARCbioinfo/pipeline-nf:v2.2
```
This will create a singularity container file: pipeline-nf_v2.2.sif that you can then use by specifying it in the nextflow command (see usage)

### Configuration file

## Usage

```bash
nextflow run iarcbioinfo/pipeline_name -r X --input_folder xxx --output_folder xxx -params-file xxx.yml -w /scratch/work
```
OR USING SINGULARITY
```bash
nextflow run iarcbioinfo/pipeline_name -r X -profile singularity --input_folder xxx --output_folder xxx -params-file xxx.yml -w /scratch/work
```
OR USING SINGULARITY WITH SPECIFIC CONTAINER
```bash
nextflow run iarcbioinfo/pipeline_name -r X -with-singularity XXX.sif --input_folder xxx --output_folder xxx -params-file xxx.yml -w /scratch/work
```

### Pipelines updates

You can update the nextflow sofware and the pipeline itself simply using:
```bash
nextflow -self-update
nextflow pull iarcbioinfo/pipeline_name
```

You can also automatically update the pipeline when you run it by adding the option `-latest` in the `nextflow run` command. Doing so you will always run the latest version from Github.

### Display help

```bash
nextflow run iarcbioinfo/pipeline_name --help
```

