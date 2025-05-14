<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-ontgeno_logo_dark.png">
    <img alt="nf-core/ontgeno" src="docs/images/nf-core-ontgeno_logo_light.png">
  </picture>
</h1>[![GitHub Actions CI Status](https://github.com/nf-core/ontgeno/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/ontgeno/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/ontgeno/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/ontgeno/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/ontgeno/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/ontgeno)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23ontgeno-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/ontgeno)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/ontgeno** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->
1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
3. Map reads to reference genome with minimap2
4. Perform variant calling with Clair3, DeepVariant, and GATK4 HaplotypeCaller
5. Joint genotype calling with GLNexus
6. Phasing with Eagle2 and Whatshap
7. Variant annotation with Ensembl VEP
8. Structural variant calling with Sniffles2 and cuteSV

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate): -->

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,runid,library,fastq
patient01,run1,A,/path_to_data/patient01.fastq.gz
```

Each row represents a fastq file. Files with the same sample id are merged during the analysis. The combination of sample, runid and library should be unique for each file.

Next, prepare a YAML file for the pipeline parameters, referring to the samplesheet with the `input` parameter:

`params.yaml`:
```yaml
# General settings:
#-----------------------
project: "p2024-0032"
outdir: "results_sup"

# Sample details:
#-----------------------
input: "test/samplelist.sup.csv"

# References and intervals:
#-----------------------
fasta: "/data/references/Homo_sapiens/Ensembl/GRCh38/Sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
bed: "test/region.bed"
#str_file: "assets/human_GRCh38_no_alt_analysis_set.trf.bed"
panel: "assets/1kGP_high_coverage_Illumina.7_117470000_117680000.filtered.SNV_INDEL_SV_phased_panel.bcf"

# Workflow settings:
#-----------------------
basecalling_model: 'sup'
trim_length: 0
store_trimmed: false
min_coverage: 10
genotype_model: "assets/r1041_e82_400bps_sup_v410"
glnexus_config: "assets/clair3.yml"

# Annotation settings:
#-----------------------
genome: 'GRCh38'
species: 'homo_sapiens'
vep_cache_version: 112
vep_cache: "/data/databases/vep"
```

### Running on the IBU cluster

To run the pipeline, start an interactive session on the IBU cluster:

```bash
srun --partition pibu_el8 --account p2024-0032 --cpus-per-task=1 --mem=8000 --time=144:00:00 --pty bash
module load Java
export NXF_SINGULARITY_CACHEDIR=/data/projects/p2024-0032_comprehensive_cftr_gene_sequencing_using_long_read_nanopore_technology/pipelines/singularity_cache
export NXF_TEMP=$SCRATCH
```

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run main.nf \
  -profile unibe_ibu \
  -params-file test/params.sup.yaml
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/ontgeno/usage) and the [parameter documentation](https://nf-co.re/ontgeno/parameters).


### Local execution

To run the pipeline locally on your computer, use the provided `local.config` file:

```bash
nextflow run main.nf \
  -c local.config \
  -params-file test/params.sup.yaml
```

Make sure to set the environment variable `NXF_SINGULARITY_CACHEDIR` to avoid having to download Singularity containers repeatedly:
```bash
echo NXF_SINGULARITY_CACHEDIR
```

If not defined, set it before running the pipeline to a directory with sufficient free space:
```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_cache
```

### Clean-up

To completely clean up all previous pipeline runs, the following files and folders need to be deleted:

```bash
rm -r work            # The working cache of the pipeline
rm -r .nextflow       # The Nextflow database containing information for nextflow log
rm -r <outdir>        # The folder with published results as defined in --outdir
rm .nextflow.log*     # The log files of previous runs
```
After this, it is no longer possible to resume a previous run, so be careful what you delete. To just delete a specific pipeline run, use `nextflow log` to get the ids of all runs and remove the run with `nextflow clean <RUN_NAME> -f`.

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/ontgeno/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/ontgeno/output).

## Credits

nf-core/ontgeno was originally written by Alexander Nater.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#ontgeno` channel](https://nfcore.slack.com/channels/ontgeno) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/ontgeno for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) --><!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
