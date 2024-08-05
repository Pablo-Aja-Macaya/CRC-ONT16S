## Summary
This repository contains the explanations/code to reproduce the bioinformatic analyses done in [**xxxxxxxx**](https://doi.org/xxxxxxxxx), where Oxford Nanopore Technologies (ONT) reads are used to study the microbiome of colorectal cancer patients (CRC) through full 16S rRNA sequencing (V1V9). Additionally, three basecalling models (fast, hac and sup), two databases (rrnDB+NCBI-16S-RefSeq and SILVA) and two different approaches (ONT-V1V9 and Illumina-V3V4) are compared. **Main analysis is shown in HTML report inside `analysis_20240705.zip`**.

> [!IMPORTANT]
> Nanopore POD5 signals are classified and analyzed through the scripts shown here, but the Illumina reads used for comparison are classified using the code available in [**this repository**](https://github.com/Pablo-Aja-Macaya/CRC-16S-study).

<!-- ## Index
- [Environments](#environments)
- [Input and metadata](#input-and-metadata)
- [Post processing](#post-processing)
- [Other analyses](#other-analyses) -->

## Environments
Each step has different requirements, but a base environment is used to run the Snakemake and Python scripts. A Conda/Mamba installation is required (with mamba being highly recommended).

```sh
# Main env (scripts are called from here)
mamba create -n snakemake_env
conda activate snakemake_env
mamba install snakemake==7.18.2 python==3.8.10 pandas==1.5.0 numpy==1.23.1 colorama matplotlib seaborn
```

The rest of the environments are available at `envs/*.yml`. They are all used by `classification_pipeline.smk` at variable definition `CONDA_ENVS` (beginning of script) except `ancombc2.yml`, which is used by `analysis.R`. Additionally, an installation of dorado and its models is required for `classification_pipeline.smk`, and is also defined at the beginning of the script.

## Input and metadata

Basecalled Oxford Nanopore and Illumina reads are available in bioproject **PRJNA911189**. Metadata can be downloaded from the bioproject (Using both biosample and SRA metadata), but it is also available in this repository at `data/metadata_table.tsv`. POD5 files are not available in NCBI at the moment, which means `classification_pipeline.smk` can not be executed, but its output is provided (see [Post processing](#post-processing)).

The expected metadata for the following scripts consists in a .tsv file with the columns sample-id, subject, origin, sample-type, sample-nature and sequencing-run. Where origin is the type of sample (subgingival-fluid, saliva, faeces...), sample-type is the subject classification (crc or non-crc). ONT samples will additionally have a barcode column.

## Running classification on ONT reads

In order to run `classification_pipeline.smk` POD5s are needed (not available for download at the moment but the pipeline output is provided, see [Post processing](#post-processing)) and the variables in the `Config variables` section must be configured.

Then, the script can be ran like this:

```sh
# Activate main environment
conda activate snakemake_env 
threads="60"

# Check input is OK
snakemake -s classification_pipeline.smk -c $threads --use-conda --keep-going -n

# Run
snakemake -s classification_pipeline.smk -c $threads --use-conda --keep-going
```


## Post processing
A script has been developed to format and clean data and to compare alpha-diversity, beta-diversity and differential abundance. This script is available at `analysis.Rmd` and its HTML result is contained in `analysis_20240705.zip`. Additionally, input for `analysis.Rmd` is provided in `data/`, which is the result of processing ONT reads through `classification_pipeline.smk` and merging the output with Illumina's results:
- `metadata_table.tsv` 
- `taxonomy_table.tsv`
- `feature_table.json`

<!-- <p align="center">
  <img src="post-processing/examples.png">
</p> -->

## Other analyses
The explanation for other analyses done are available at the Methodology section of the paper under "Bioinformatic analysis". 
