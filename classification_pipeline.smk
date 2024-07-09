import sys
sys.path.append(".")

"""
Steps
- [x] Convert to pod
- [x] Split pod5s into channels
- [x] Basecall each run with dorado duplex
- [x] Split reads, searching middle strand adaptors/barcodes
- [x] Separate into simplex and duplex reads
    - Reasons: https://github.com/nanoporetech/dorado/issues/600#issuecomment-1943011355, https://github.com/nanoporetech/dorado/issues/679#issuecomment-1988583542
    - Demultiplex simplex reads (with trimming)
    - Assign barcode to duplex reads according to the parent reads barcode. There can be duplex reads with parents from different barcodes, remove them
    - Merge duplex and simplex into one file per barcode
- [x] Switch barcode names to sample-id names
- [x] Filter by quality, min length, max length. Optional: Cut from front and back?
- [x] Emu
"""

###############################################
# ------- Libraries/Functions set up -------- #
###############################################
import pandas as pd
from pathlib import Path
import glob
import numpy
import shutil
from subprocess import run
from tqdm import tqdm
from Bio import SeqIO
import re
from colorama import Fore, Style
from snakemake.logging import logger
import plotly.express as px
from functools import reduce

# Environments
CONDA_ENVS = {
    "pod5": "pod5",
    "chopper": "chopper",
    "nanoplot": "bact_ont",
    "emu": "emu",
    "qiime2": "qiime2-2022.8",
    "duplex_tools": "bact_ont",
    "pychopper": "pychopper",
}

###############################################
# ------------ Config variables ------------- #
###############################################
# -- Inputs --
METADATA = "/home/usuario/Proyectos/Results/ONT16S/data/metadata.tsv"
SEQUENCING_RUN_IDS = {
    "20230621_1356_MC-112596_FAW78429_3ad37fea": {
        "pod5_folder": "/media/usuario/15d7f455-39ea-40e7-add2-de57c58767eb/Archive/Nanopore/ONT16S/16SFULLCCR/no_sample/20230621_1356_MC-112596_FAW78429_3ad37fea/pod5" #"/home/usuario/Proyectos/Results/tests/pod5_0" 
    },
    "20240312_1555_MC-112596_FAY65338_e65cc192": {
        "pod5_folder": "/media/usuario/15d7f455-39ea-40e7-add2-de57c58767eb/Archive/Nanopore/ONT16S/full16S_120424/no_sample/20240312_1555_MC-112596_FAY65338_e65cc192/pod5" #"/home/usuario/Proyectos/Results/tests/pod5" # 
    }, 
}

# -- Params --
# Basecalling
DORADO = "/home/usuario/Proyectos/programas/dorado/dorado-0.5.3-linux-x64/bin/dorado"
BARCODE_KITS = "SQK-NBD114-96"

dorado_models_folder = "/home/usuario/Proyectos/programas/dorado/dorado-0.5.3-linux-x64/models"
DORADO_MODELS = {
    "dna_r10.4.1_e8.2_400bps_fast@v4.1.0": {
        "model": f"{dorado_models_folder}/dna_r10.4.1_e8.2_400bps_fast@v4.1.0",
        "short_name": "fast",
        "minimum_quality": 7 
    },
    "dna_r10.4.1_e8.2_400bps_hac@v4.1.0": {
        "model": f"{dorado_models_folder}/dna_r10.4.1_e8.2_400bps_hac@v4.1.0",
        "short_name": "hac",
        "minimum_quality": 12
    },
    "dna_r10.4.1_e8.2_400bps_sup@v4.1.0": {
        "model": f"{dorado_models_folder}/dna_r10.4.1_e8.2_400bps_sup@v4.1.0",
        "short_name": "sup",
        "minimum_quality": 12
    },
}
PRIMERS_CONFIG = "/home/usuario/Proyectos/Results/ONT16S/data/primer_config.txt"
PRIMERS_FASTA = "/home/usuario/Proyectos/Results/ONT16S/data/primers.fasta"

# Emu
EMU_DATABASES = {
    "silva": "/home/usuario/Proyectos/Innova/data/db/emu/emu_silva_db",
    "default": "/home/usuario/Proyectos/Innova/data/db/emu/emu_default_db",
}

# -- Illumina --
# These are tables calculated outside of this script using Illumina V3V4 for the same samples
illu_folder = "/media/usuario/15d7f455-39ea-40e7-add2-de57c58767eb/analyses/Kelly/KellyCCR/analisis_ccr_kelly_2022-10-25_noreads/qiime/export_results"
ILLUMINA_DATA_INPUT = {
    "feature_table": f"{illu_folder}/feature-table.txt",
    "taxonomy_table": f"{illu_folder}/taxonomy.tsv",
    "metadata_table": f"/home/usuario/Proyectos/CRC-16S-study/data/non-ffpe_metadata.tsv",
}



###############################################
# -------------- Folder set up -------------- #
###############################################
# Main
BASE_FOLDER = Path("/home/usuario/Proyectos/Results/ONT16S/pipeline")

# Output paths
SPLIT_POD5_OUTPUT = f"{BASE_FOLDER}/01_split_pod5s"
DORADO_OUTPUT = f"{BASE_FOLDER}/02_dorado"
NANOPLOT_OUTPUT = f"{BASE_FOLDER}/02.2_nanoplot"
SEPARATE_READ_TYPES_OUTPUT = f"{BASE_FOLDER}/03_separate_read_types"
DEMUX_SIMPLEX_OUTPUT = f"{BASE_FOLDER}/04_demux_simplex"
FILTER_DUPLEX_OUTPUT = f"{BASE_FOLDER}/05_filter_duplex"
GET_DEMUX_READS_OUTPUT = f"{BASE_FOLDER}/06_demuxed_reads"
CHOPPER_OUTPUT = f"{BASE_FOLDER}/07_chopper"
SPLIT_ON_ADAPTER_OUTPUT = f"{BASE_FOLDER}/08_split_on_adapter"
DETECT_PRIMERS_OUTPUT = f"{BASE_FOLDER}/08_primer_detection"
FILTER_16S_CONCATAMERS = f"{BASE_FOLDER}/09_filter_16s_concatamers"
EMU_OUTPUT = f"{BASE_FOLDER}/10_emu"
MERGE_EMU_OUTPUT = f"{BASE_FOLDER}/11_merge_emu"
MERGE_ALL_EMU_TABLES_PER_GROUP = f"{BASE_FOLDER}/12_merge_all_emu_tables_per_group"
MERGE_ALL_EMU_TABLES = f"{BASE_FOLDER}/13_merge_all_emu_tables"
MERGE_ONT_AND_ILLUMINA = f"{BASE_FOLDER}/14_merge_ont_and_illumina"
BIOM_CONVERSION = f"{BASE_FOLDER}/15_biom_conversion"

###############################################
# ----------------- Input ------------------- #
###############################################

def check_sequencing_run_ids_in_metadata(metadata_file=METADATA, sequencing_run_ids=list(SEQUENCING_RUN_IDS.keys()), seq_run_col="sequencing-run"):
    df = pd.read_csv(METADATA, sep="\t")
    
    assert len(df[seq_run_col].isin(sequencing_run_ids))>=1, f"No rows in column '{seq_run_col}' metadata with these sequencing run ids: {sequencing_run_ids}"

    temp = df[~df[seq_run_col].isin(sequencing_run_ids)]
    if len(temp) != 0:
        logger.info(f"WARNING: There are {len(temp)} rows in metadata with unexpected sequencing run ids:\n{Fore.LIGHTBLACK_EX}{temp}{Style.RESET_ALL}")

    df = df[df[seq_run_col].isin(sequencing_run_ids)]

    return df


def check_no_repeated_samples(metadata_df):
    # Duplicated sample-ids
    duplicated = metadata_df[metadata_df["sample-id"].duplicated()]
    assert len(duplicated) == 0, f"There are duplicated sample-ids: {duplicated}"

    # Group by sequencing-run and ensure no barcode repetition inside sequencing-run
    sampleids_per_bc = metadata_df.groupby(["sequencing-run","barcode"])[["sample-id"]].count().rename(columns={"sample-id":"sample_id_count"})
    sampleids_per_bc = sampleids_per_bc[sampleids_per_bc["sample_id_count"]!=1]
    assert len(sampleids_per_bc)==0, f"There are barcodes in a sequencing-run which refer to multiple sample-ids: {sampleids_per_bc}" 

    return metadata_df

METADATA_DF = check_sequencing_run_ids_in_metadata()
METADATA_DF = check_no_repeated_samples(METADATA_DF)


###############################################
# ---------------- Output ------------------- #
###############################################

# def get_output(wc, metadata_df=METADATA_DF):
#     d = {}
#     for record in metadata_df.to_dict("records"):
#         sequencing_run_id = record["sequencing-run"]
#         barcode = record["barcode"]
#         fs = [
#             Path(GET_DEMUX_READS_OUTPUT, sequencing_run_id, basecalling_model_id, barcode, f"{barcode}.fastq.gz"),
#             Path(EMU_OUTPUT, sequencing_run_id, basecalling_model_id, barcode, f"{barcode}_read-assignment-distributions.tsv"),
#             Path(NANOPLOT_OUTPUT, sequencing_run_id, basecalling_model_id, barcode),
#             Path(MERGE_EMU_OUTPUT, sequencing_run_id, basecalling_model_id, "raw_feature_table.tsv"),
#         ]
#         if d.get(sequencing_run_id):
#             d[sequencing_run_id] += [*fs]
#         else:
#             d[sequencing_run_id] = [*fs]

#     fs = [j for i in d.values() for j in i]
#     return fs
        

rule all:
    input:
        # lambda wc: get_output(wc),
        expand(Path(DORADO_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{sequencing_run_id}__{basecalling_model_id}.fastq.gz"), sequencing_run_id=list(SEQUENCING_RUN_IDS.keys()), basecalling_model_id=list(DORADO_MODELS.keys())),
        # expand(Path(NANOPLOT_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}"), sequencing_run_id=list(SEQUENCING_RUN_IDS.keys()), basecalling_model_id=list(DORADO_MODELS.keys())),
        expand(Path(MERGE_EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "raw_feature_table.tsv"), sequencing_run_id=list(SEQUENCING_RUN_IDS.keys()), basecalling_model_id=list(DORADO_MODELS.keys()), emu_database_id=list(EMU_DATABASES.keys())),
        expand(Path(MERGE_ALL_EMU_TABLES_PER_GROUP, "{basecalling_model_id}", "{emu_database_id}", "feature_table.tsv"), basecalling_model_id=list(DORADO_MODELS.keys()), emu_database_id=list(EMU_DATABASES.keys())),
        Path(MERGE_ALL_EMU_TABLES, "feature_table.tsv"),
        Path(MERGE_ONT_AND_ILLUMINA, "feature_table.tsv"),
        Path(BIOM_CONVERSION, "feature_table.biom.json"),


# ---- Convert ----
# "pod5 convert fast5 [INPUT_FOLDER] -t 60 -o [OUTPUT_FOLDER] --one-to-one"
# "pod5 convert fast5 [INPUT_FOLDER] -o [OUTPUT_FOLDER] --one-to-one [INPUT_FOLDER]"


rule split_pod5s_per_channel:
    """
    Split folder of pod5s into a pod5 per channel
    This makes basecalling faster and allows for separate basecalling of each file maintaining
    duplex consistency (if not done like this duplex reads can be lost due to being in different pod5s)
    """
    input:
        pod5_folder = lambda wc: SEQUENCING_RUN_IDS[wc.sequencing_run_id]["pod5_folder"],
    output:
        folder = directory(Path(SPLIT_POD5_OUTPUT, "{sequencing_run_id}")),
        split_pod5_folder = temp(directory(Path(SPLIT_POD5_OUTPUT, "{sequencing_run_id}", "split_pod5s"))),
        channel_summary = Path(SPLIT_POD5_OUTPUT, "{sequencing_run_id}", "summary.tsv"),
    conda: CONDA_ENVS["pod5"]
    threads: 60
    shell:
        """
        pod5 view {input.pod5_folder} --include "read_id, channel" --output {output.channel_summary} --threads {threads}
        pod5 subset {input.pod5_folder} --summary {output.channel_summary} --columns channel --output {output.split_pod5_folder} --threads {threads}
        """

rule dorado:
    """
    Basecalling with dorado
    """
    input:
        pod5_folder = rules.split_pod5s_per_channel.output.split_pod5_folder,
    output:
        folder = directory(Path(DORADO_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}")),
        reads = protected(Path(DORADO_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{sequencing_run_id}__{basecalling_model_id}.fastq.gz")),
    params:
        dorado = DORADO,
        model = lambda wc: DORADO_MODELS[wc.basecalling_model_id]["model"],
        batchsize = 64, # 64 for low ram gpus
    priority: 999
    threads: 60
    shell:
        """
        {params.dorado} duplex --emit-fastq -b 64 {params.model} {input.pod5_folder} | pigz -p {threads} > {output.reads}
        """

rule nanoplot:
    input:
        reads = rules.dorado.output.reads,
    output:
        folder = directory(Path(NANOPLOT_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}")),
    conda: CONDA_ENVS["nanoplot"]
    threads: 12
    shell:
        """
        NanoPlot --fastq {input.reads} -o {output.folder} -t {threads}
        """


rule separate_read_types:
    """
    Separate read types into duplex and simplex
    Duplex reads IDs are formed by two simplex parental IDs, separated by ";"
    """
    input:
        reads = rules.dorado.output.reads,
    output:
        folder = directory(Path(SEPARATE_READ_TYPES_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}")),
        simplex_reads = temp(Path(SEPARATE_READ_TYPES_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "simplex.fastq")),
        simplex_reads_ids = Path(SEPARATE_READ_TYPES_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "simplex_ids.txt"),
        duplex_reads = temp(Path(SEPARATE_READ_TYPES_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "duplex.fastq")),
        duplex_reads_ids = Path(SEPARATE_READ_TYPES_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "duplex_ids.txt"),
    threads: 60
    shell:
        """
        mkdir -p {output.folder}

        # Separate into duplex and simplex
        seqkit seq -j {threads} -n {input.reads} | awk '/;/ {{print $0 > "{output.duplex_reads_ids}"; next}}{{print $0 > "{output.simplex_reads_ids}"}}'

        # Select
        seqkit grep -j {threads} --pattern-file {output.simplex_reads_ids} {input.reads} > {output.simplex_reads}
        seqkit grep -j {threads} --pattern-file {output.duplex_reads_ids} {input.reads} > {output.duplex_reads}
        """

rule demux_simplex:
    """
    Demultiplex simplex reads
    """
    input:
        reads = rules.separate_read_types.output.simplex_reads,
    output:
        folder = directory(Path(DEMUX_SIMPLEX_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}")),
    threads: 60
    params:
        dorado = DORADO,
        barcoding_kit = BARCODE_KITS,
    shell:
        """
        # Demux
        {params.dorado} demux --kit-name {params.barcoding_kit} --output-dir {output.folder} --emit-fastq --threads {threads} {input.reads}

        # Rename output reads to barcodeNN.fastq instead of [params.barcoding_kit]_barcodeNN.fastq
        rename 's/{params.barcoding_kit}_//' {output.folder}/*.fastq
        """

rule filter_duplex:
    """
    Remove duplex reads where parental simplex reads are discordant in barcode or length
    """
    input:
        simplex_reads_demux_folder = rules.demux_simplex.output.folder,
        duplex_reads = rules.separate_read_types.output.duplex_reads,
        duplex_reads_ids = rules.separate_read_types.output.duplex_reads_ids,
    output:
        folder = directory(Path(FILTER_DUPLEX_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}")),
        duplex_qc = Path(FILTER_DUPLEX_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "duplex_qc.tsv"),
        accepted_duplex_ids_per_bc_folder = directory(Path(FILTER_DUPLEX_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "accepted_duplex_per_bc")),
    threads: 60
    run:
        simplex_demux_folder = input.simplex_reads_demux_folder
        duplex_reads_ids = input.duplex_reads_ids
        duplex_reads_fastq = input.duplex_reads

        # -- Get assignation of each simplex read --
        l = []
        barcodes_seen = []
        for f in glob.glob(f"{simplex_demux_folder}/*.fastq"):
            bc = Path(f).name.replace(".fastq","")

            for record in SeqIO.parse(f, "fastq"):
                l.append([record.id, bc, len(record.seq)])

            if bc not in barcodes_seen:
                barcodes_seen.append(bc)

        simplex_reads_demux = pd.DataFrame(l)
        simplex_reads_demux.columns = ["read_id", "barcode", "length"]
        simplex_reads_demux

        # -- Now get duplex reads and their parents --
        duplex_reads = pd.read_csv(duplex_reads_ids, header=None)
        duplex_reads.columns = ["read_id"]
        duplex_reads[["duplex_parent_1", "duplex_parent_2"]] = duplex_reads[
            "read_id"
        ].str.split(";", expand=True)

        # -- Check if parents barcode are the same --
        df = pd.merge(
            duplex_reads, simplex_reads_demux, left_on="duplex_parent_1", right_on="read_id"
        )
        df = pd.merge(df, simplex_reads_demux, left_on="duplex_parent_2", right_on="read_id")
        df = df[
            [
                "read_id_x",
                "duplex_parent_1",
                "duplex_parent_2",
                "barcode_x",
                "barcode_y",
                "length_x",
                "length_y",
            ]
        ]
        df = df.rename(
            columns={
                "read_id_x": "read_id",
                "barcode_x": "duplex_parent_1_barcode",
                "barcode_y": "duplex_parent_2_barcode",
                "length_x": "duplex_parent_1_length",
                "length_y": "duplex_parent_2_length",
            }
        )

        # -- Check quality --
        # Parents with equal barcodes
        df["consensus_barcode"] = df.apply(
            lambda row: row["duplex_parent_1_barcode"]
            if row["duplex_parent_1_barcode"] == row["duplex_parent_2_barcode"]
            else "BAD",
            axis=1,
        )
        # Lengths not differing more than 50%
        df["consensus_length"] = df.apply(
            lambda row: None
            if abs(100 - 100 * row["duplex_parent_1_length"] / row["duplex_parent_2_length"])
            <= 50
            else "big_difference",
            axis=1,
        )

        # -- Save stats per duplex read --
        df.to_csv(
            output.duplex_qc,
            sep="\t",
            index=None,
        )

        # Apply filters
        filtered_duplex = df[
            (df["consensus_barcode"] != "BAD")
            & (df["consensus_length"] != "big_difference")
        ]
        # filtered_duplex["read_id"].to_csv(
        #     output.accepted_duplex_ids,
        #     sep="\t",
        #     index=None,
        #     header=None
        # )

        # -- Save all passing duplex reads into one file per barcode --
        shell("mkdir -p {output.accepted_duplex_ids_per_bc_folder}")
        for bc in barcodes_seen:
            temp = filtered_duplex[filtered_duplex["consensus_barcode"] == bc]
            simplex = f"{simplex_demux_folder}/{bc}.fastq"
            temp["read_id"].to_csv(f"{output.accepted_duplex_ids_per_bc_folder}/{bc}.txt", sep="\t", index=None, header=None)

        # passing_reads = pd.concat(
        #     [
        #         # Simplex reads
        #         simplex_reads_demux[["read_id", "barcode"]],
        #         # Duplex reads passing QC
        #         filtered_duplex[["read_id", "consensus_barcode"]].rename(
        #             columns={"consensus_barcode": "barcode"}
        #         ),
        #     ]
        # )


rule get_demux_reads:
    input:
        simplex_reads_demux_folder = rules.demux_simplex.output.folder,
        duplex_reads = rules.separate_read_types.output.duplex_reads,
        duplex_reads_ids_per_bc_folder = rules.filter_duplex.output.accepted_duplex_ids_per_bc_folder,
    output:
        folder = directory(Path(GET_DEMUX_READS_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}")),
        reads = Path(GET_DEMUX_READS_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}", "{barcode}.fastq.gz"),
    threads: 12
    shell:
        """
        demuxed_simplex_fastq="{input.simplex_reads_demux_folder}/{wildcards.barcode}.fastq"
        duplex_reads_ids="{input.duplex_reads_ids_per_bc_folder}/{wildcards.barcode}.txt"

        # Add simplex reads to file
        cat $demuxed_simplex_fastq > {output.reads}.temp

        # Extract duplex reads
        seqkit grep --pattern-file $duplex_reads_ids {input.duplex_reads} >> {output.reads}.temp

        # Compress
        pigz --stdout --processes {threads} {output.reads}.temp > {output.reads}
        rm {output.reads}.temp
        """


rule chopper:
    input:
        reads = rules.get_demux_reads.output.reads,
    output:
        folder = directory(Path(CHOPPER_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}")),
        reads = temp(Path(CHOPPER_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}", "{barcode}.fastq.gz")),
    params:
        quality = lambda wc: DORADO_MODELS[wc.basecalling_model_id]["minimum_quality"],
        headcrop = 20,
        tailcrop = 20,
        minlength = 1200,
        maxlength = 1900,
    conda: CONDA_ENVS["chopper"]
    threads: 6
    shell:
        """
        mkdir -p {output.folder}
        pigz -p {threads} --decompress --stdout {input.reads} | \
            chopper -q {params.quality} --headcrop {params.headcrop} --tailcrop {params.tailcrop} --minlength {params.minlength} --maxlength {params.maxlength} -t {threads} 2> {output.folder}/chopper.log | \
            pigz -p {threads} > {output.reads}
        echo "quality:{params.quality}; headcrop:{params.headcrop}; tailcrop:{params.tailcrop}; minlength:{params.minlength}; maxlength:{params.maxlength}" >> {output.folder}/chopper.log
        """

rule split_on_adapter:
    input:
        folder = rules.chopper.output.folder,
        reads = rules.chopper.output.reads,
    output:
        folder=directory(Path(SPLIT_ON_ADAPTER_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}")),
        reads=temp(Path(SPLIT_ON_ADAPTER_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}", "{barcode}.fastq.gz")),
        split_reads=Path(SPLIT_ON_ADAPTER_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}", "split_reads.txt"),
    params:
        sample_type="Native",
    threads: 6
    conda:
        CONDA_ENVS["duplex_tools"]
    shell:
        """
        # Split on adapter
        duplex_tools split_on_adapter {input.folder} {output.folder}/split {params.sample_type} --threads {threads} --allow_multiple_splits > {output.folder}/duplex_tools.log
        
        # Grab read ids that end in _1, _2 (these are the ones with adapters in the midlle originally), deduplicate the ids and save into file
        seqkit grep -r -i -p "_\d+$" {output.folder}/split/*.fastq.gz | seqkit replace -p "(.*)_\d+" -r '$1' | seqkit rmdup | seqkit seq -n -i > {output.split_reads}

        # Grab all reads that did not contain a adapter
        seqkit grep --invert-match --pattern-file {output.split_reads} {output.folder}/split/*.fastq.gz > {output.reads}
        
        # Clean
        rm -r {output.folder}/split
        """

# rule detect_primers:
#     """
#     Detect specific primers
#     pychopper is designed for simple primer pairs. If given a lot of primers the tool will detect them and write it in output.bed
#     but it will think that anything outside of reads with two primers is wrong. We only care about the bed file where positions for each primer are given
#     """
#     input:
#         reads = rules.chopper.output.reads,
#     output:
#         folder = directory(Path(DETECT_PRIMERS_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}")),
#         bed = Path(DETECT_PRIMERS_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}", "scores.bed"),
#     params:
#         primers_config = PRIMERS_CONFIG,
#         primers_fasta = PRIMERS_FASTA,
#         minimum_segment_length = 500,
#     threads: 10
#     conda: CONDA_ENVS["pychopper"]
#     shell:
#         """
#         mkdir -p {output.folder}

#         pychopper -m edlib -p -z {params.minimum_segment_length} -t {threads} \
#             -b {params.primers_fasta} -c {params.primers_config} \
#             -r {output.folder}/report.pdf -S {output.folder}/statistics.tsv -A {output.bed} \
#             -u {output.folder}/unclassified.fastq -w {output.folder}/rescued.fastq \
#             {input.reads} {output.folder}/output.fastq &> {output.folder}/pychopper.log

#         rm {output.folder}/*.fastq {output.folder}/*.pdf
#         """

# rule filter_16s_concatamers:
#     input:
#         reads = rules.chopper.output.reads,
#         primer_bed = rules.detect_primers.output.bed,
#     output:
#         folder = directory(Path(FILTER_16S_CONCATAMERS, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}")),
#         reads = temp(Path(FILTER_16S_CONCATAMERS, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}", "{barcode}.fastq.gz")),
#         concatamers_ids = Path(FILTER_16S_CONCATAMERS, "{sequencing_run_id}", "{basecalling_model_id}", "{barcode}", "concatamers_ids.txt"),
#     threads: 6
#     params:
#         counted_primers_minimum=5, # minimum number of primers to be detected, reads with less hits are removed
#         repeated_primers_maximum=4, # maximum number of primers repeated in read
#         repeat_search_window_separation=100, # separation between repeated primer hits to be considered a repeat (there can be multiple hits in the same position, but having hits very separated can indicate concatamers)
#     run:
#         shell("mkdir -p {output.folder}")

#         from scripts.check_concatamers import check_concatamers
#         _ = check_concatamers(input.primer_bed, output.concatamers_ids, params.counted_primers_minimum, params.repeated_primers_maximum, params.repeat_search_window_separation)

#         shell("seqkit grep -j {threads} --invert-match --pattern-file {output.concatamers_ids} {input.reads} > {output.reads}")


rule emu:
    input:
        reads = rules.split_on_adapter.output.reads,
    output:
        folder = directory(Path(EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "{barcode}")),
        relative_abundance = Path(EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "{barcode}", "{barcode}_rel-abundance.tsv"),
        read_assignment = Path(EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "{barcode}", "{barcode}_read-assignment-distributions.tsv"),
    params:
        db= lambda wc : EMU_DATABASES[wc.emu_database_id],
        read_type="map-ont",
    conda: CONDA_ENVS["emu"]
    threads: 12
    shell:
        """
        emu abundance --db {params.db} --type {params.read_type} \
            --keep-counts --keep-read-assignments --output-dir {output.folder} \
            --threads {threads} --output-basename {wildcards.barcode} {input.reads}
        """


# Propagate unidentified levels
# "g__Klebsiella";"s__" --> "g__Klebsiella;s__Klebsiella NA"
def propagate_levels(row):
    previous_tax = ""
    for k, v in row.items():
        if not v:
            v = f"{previous_tax}"
            if not v.endswith("NA"):
                v += " NA"
        row[k] = v
        previous_tax = v
    return row

rule merge_emu:
    input:
        emu_res = lambda wc: expand(
            rules.emu.output.relative_abundance, 
            sequencing_run_id=wc.sequencing_run_id, 
            basecalling_model_id=wc.basecalling_model_id, 
            emu_database_id=wc.emu_database_id, 
            barcode=METADATA_DF[METADATA_DF["sequencing-run"]==wc.sequencing_run_id]["barcode"].to_list()
        )
    output:
        folder = directory(Path(MERGE_EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}")),
        relative_feature_table = Path(MERGE_EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "relative_feature_table.tsv"),
        raw_feature_table = Path(MERGE_EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "raw_feature_table.tsv"),
        taxonomy_table = Path(MERGE_EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "taxonomy_table.tsv"),
        relative_abundance_plot = Path(MERGE_EMU_OUTPUT, "{sequencing_run_id}", "{basecalling_model_id}", "{emu_database_id}", "relative_abundance.html"),
    run:
        emu_files = input.emu_res

        # -- Read and merge emu results --
        def merge_emu_results(emu_files, measure, db=wildcards.emu_database_id):
            # ---- Merge emu results ----
            l = []
            for i in emu_files:
                # Find barcode id
                bc = Path(i).name.replace("_rel-abundance.tsv", "")

                # Read
                df = pd.read_csv(i, sep="\t")
                
                # Rename columns
                df = df.rename(columns={
                    "tax_id": "feature_id",
                    "abundance": "relative_abundance",
                    "estimated counts": "raw_abundance",
                })

                # If database is default compress taxonomy into column lineage
                # (Silva comes by default in column "lineage")
                if db=="default":
                    cols = [
                        'superkingdom', 
                        'phylum',
                        'class', 
                        'order', 
                        'family',
                        'genus', 
                        'species', 
                        'subspecies',
                    ]
                    df["lineage"] = df[cols].fillna("").astype(str).T.apply(lambda x: x.str.cat(sep=";"))
                    df = df.drop(cols, axis=1)

                # Transform results
                df["relative_abundance"] = 100*df["raw_abundance"]/df["raw_abundance"].sum()
                df["raw_abundance"] = df["raw_abundance"].astype(int)
                df.loc[df["feature_id"]=="unassigned","lineage"] = "[Unassigned]"
                
                # Rename selected measure with the name of the barcode
                df = df.rename(columns={measure: bc})

                # Select columns and store
                df = df[["feature_id","lineage", bc]]
                l.append(df)

            main_df = l[0]
            for i in l[1:]:
                main_df = pd.merge(main_df, i, how="outer", on=["feature_id","lineage"])
            main_df = main_df.fillna(0)

            # -- Separate taxonomy and feature table --
            tax_df = main_df[["feature_id","lineage"]]
            main_df = main_df.drop("lineage", axis=1)

            # -- Edit taxonomy table --
            levels = [
                "superkingdom",
                "phylum", 
                "class", 
                "order", 
                "family", 
                "genus", 
                "species", 
            ]
            tax_df[levels + ["extra"]] = tax_df["lineage"].str.split(";", expand=True)
            tax_df = tax_df.fillna(value='')
            tax_df = tax_df.drop(["lineage","extra"], axis=1)

            # Add genera to species
            if db=="silva":
                tax_df["species"] = tax_df.apply(lambda r: f"{r['genus']} {r['species']}" if r["genus"] and r["species"] else "", axis=1)

            tax_df[levels] = tax_df[levels].apply(propagate_levels, axis=1)

            return main_df, tax_df

        
        raw_feature_table, tax_df = merge_emu_results(emu_files, "raw_abundance")
        relative_feature_table, tax_df = merge_emu_results(emu_files, "relative_abundance")

        # -- Save tables --
        raw_feature_table.to_csv(output.raw_feature_table, sep="\t", index=None)
        relative_feature_table.to_csv(output.relative_feature_table, sep="\t", index=None)
        tax_df.to_csv(output.taxonomy_table, sep="\t", index=None)

        # -- Transform main to long and plot --
        temp = pd.merge(relative_feature_table, tax_df[["feature_id","species"]], on="feature_id")
        df = pd.melt(temp, id_vars=['feature_id', "species"], value_vars=[i for i in temp.columns if i.startswith("barcode")], value_name="abundance", var_name="bc")
        df = df.sort_values(["bc","species"], ascending=[True, False])
        fig = px.bar(
            df, x="bc", y="abundance", color="species", 
        )
        fig.update_layout(legend_traceorder="reversed")
        fig.update_layout(yaxis_range=(0, 100))
        fig.update_layout(
            xaxis_title="", yaxis_title="Relative abundance (%)"
        )
        with open(output.relative_abundance_plot, "wt") as handle:
            handle.write(fig.to_html())




rule merge_all_emu_tables_per_group:
    input:
        unpack(lambda wc: {
            k: [
                expand(rules.merge_emu.output.raw_feature_table, sequencing_run_id=k, basecalling_model_id=wc.basecalling_model_id, emu_database_id=wc.emu_database_id)[0],
                expand(rules.merge_emu.output.taxonomy_table, sequencing_run_id=k, basecalling_model_id=wc.basecalling_model_id, emu_database_id=wc.emu_database_id)[0],
            ] for k in SEQUENCING_RUN_IDS.keys()
        })
    output:
        folder = directory(Path(MERGE_ALL_EMU_TABLES_PER_GROUP, "{basecalling_model_id}", "{emu_database_id}")),
        feature_table = Path(MERGE_ALL_EMU_TABLES_PER_GROUP, "{basecalling_model_id}", "{emu_database_id}", "feature_table.tsv"),
        taxonomy_table = Path(MERGE_ALL_EMU_TABLES_PER_GROUP, "{basecalling_model_id}", "{emu_database_id}", "taxonomy_table.tsv"),
        metadata_table = Path(MERGE_ALL_EMU_TABLES_PER_GROUP, "{basecalling_model_id}", "{emu_database_id}", "metadata_table.tsv"),
    params:
        metadata_df = METADATA_DF,
        sequencing_run_ids=SEQUENCING_RUN_IDS,
        basecalling_model_ids=DORADO_MODELS,
    run:
        # input_dict = {
        #     "20230621_1356_MC-112596_FAW78429_3ad37fea":{
        #         "feature_table": "/home/usuario/Proyectos/Results/ONT16S/pipeline/10_merge_emu/20230621_1356_MC-112596_FAW78429_3ad37fea/raw_feature_table.tsv",
        #         "taxonomy_table": "/home/usuario/Proyectos/Results/ONT16S/pipeline/10_merge_emu/20230621_1356_MC-112596_FAW78429_3ad37fea/taxonomy_table.tsv"
        #     },
        #     "20240312_1555_MC-112596_FAY65338_e65cc192":{
        #         "feature_table": "/home/usuario/Proyectos/Results/ONT16S/pipeline/10_merge_emu/20230621_1356_MC-112596_FAW78429_3ad37fea/raw_feature_table.tsv",
        #         "taxonomy_table": "/home/usuario/Proyectos/Results/ONT16S/pipeline/10_merge_emu/20230621_1356_MC-112596_FAW78429_3ad37fea/taxonomy_table.tsv"        
        #     }
        # }

        input_dict = {id: {"feature_table": input[id][0], "taxonomy_table": input[id][1]} for id in params.sequencing_run_ids.keys()}
        metadata_df = params.metadata_df

        # Read tables from each run, transform into long and add sequencing-run column
        tax_table_l = []
        feature_table_l = []
        for run_id, results in input_dict.items():
            feature_table = pd.read_csv(results["feature_table"], sep="\t")
            taxonomy_table = pd.read_csv(results["taxonomy_table"], sep="\t")

            feature_table = pd.melt(
                feature_table, 
                id_vars='feature_id', 
                value_vars=[i for i in feature_table.columns if i!="feature_id"],
                var_name="barcode",
                value_name="abundance"
            )
            feature_table["sequencing-run"] = run_id

            feature_table_l.append(feature_table)
            tax_table_l.append(taxonomy_table)

        # Concat tables
        taxonomy_table = pd.concat(tax_table_l).drop_duplicates()
        feature_table = pd.concat(feature_table_l)

        # Pivot table to wide
        feature_table = pd.pivot(feature_table, index=["sequencing-run","barcode"], columns="feature_id").fillna(0)

        # Flatten multilevel index
        feature_table.columns = [tup[1] for tup in feature_table.columns.to_flat_index()]
        feature_table.reset_index()

        # Merge with metadata, replacing sequencing-run+barcode by sample-id
        feature_table = pd.merge(metadata_df[["sample-id","sequencing-run","barcode"]], feature_table, on=["sequencing-run","barcode"])
        feature_table = feature_table.drop(columns=['sequencing-run', 'barcode'])

        # Set sample-id as index, transpose and rename index
        feature_table = feature_table.set_index("sample-id")
        feature_table = feature_table.T
        feature_table.index.names = ['feature_id']
        feature_table = feature_table.reset_index()

        # Save
        feature_table.to_csv(output.feature_table, sep="\t", index=None)
        taxonomy_table.to_csv(output.taxonomy_table, sep="\t", index=None)
        metadata_df.to_csv(output.metadata_table, sep="\t", index=None)

rule merge_all_emu_tables:
    input:
        unpack(lambda wc: {
            f"{model_id}__{db_id}": [
                expand(rules.merge_all_emu_tables_per_group.output.feature_table, basecalling_model_id=model_id, emu_database_id=db_id)[0],
                expand(rules.merge_all_emu_tables_per_group.output.taxonomy_table, basecalling_model_id=model_id, emu_database_id=db_id)[0],
                expand(rules.merge_all_emu_tables_per_group.output.metadata_table, basecalling_model_id=model_id, emu_database_id=db_id)[0],
            ] for model_id in DORADO_MODELS.keys() for db_id in EMU_DATABASES.keys() 
        })
    output:
        folder = directory(Path(MERGE_ALL_EMU_TABLES)),
        feature_table = Path(MERGE_ALL_EMU_TABLES, "feature_table.tsv"),
        taxonomy_table = Path(MERGE_ALL_EMU_TABLES, "taxonomy_table.tsv"),
        metadata_table = Path(MERGE_ALL_EMU_TABLES, "metadata_table.tsv"),
    run:
        # Snakemake input structure
        # input = {
        #     "dna_r10.4.1_e8.2_400bps_sup@v4.1.0__default": [
        #         "/home/usuario/Proyectos/Results/ONT16S/pipeline/12_merge_all_emu_tables/dna_r10.4.1_e8.2_400bps_sup@v4.1.0/default/feature_table.tsv",
        #         "/home/usuario/Proyectos/Results/ONT16S/pipeline/12_merge_all_emu_tables/dna_r10.4.1_e8.2_400bps_sup@v4.1.0/default/taxonomy_table.tsv",
        #         "/home/usuario/Proyectos/Results/ONT16S/pipeline/12_merge_all_emu_tables/dna_r10.4.1_e8.2_400bps_sup@v4.1.0/default/metadata_table.tsv",
        #     ],
        #     "dna_r10.4.1_e8.2_400bps_sup@v4.1.0__silva": [
        #         "/home/usuario/Proyectos/Results/ONT16S/pipeline/12_merge_all_emu_tables/dna_r10.4.1_e8.2_400bps_sup@v4.1.0/silva/feature_table.tsv",
        #         "/home/usuario/Proyectos/Results/ONT16S/pipeline/12_merge_all_emu_tables/dna_r10.4.1_e8.2_400bps_sup@v4.1.0/silva/taxonomy_table.tsv",
        #         "/home/usuario/Proyectos/Results/ONT16S/pipeline/12_merge_all_emu_tables/dna_r10.4.1_e8.2_400bps_sup@v4.1.0/silva/metadata_table.tsv",
        #     ],
        # }

        # Prepare input into dict
        input_dict = {}
        for m in DORADO_MODELS.keys():
            for db in EMU_DATABASES.keys():
                i = f"{m}__{db}"
                input_dict[i] = {
                    "feature_table": input[i][0],
                    "taxonomy_table": input[i][1],
                    "metadata_table": input[i][2],
                    "model": m,
                    "database": db,
                }

        # Format for merge
        tax_tables = []
        md_tables = []
        feat_tables = []
        for k, v in input_dict.items():
            print(k)

            # -- Read --
            feature_table = pd.read_csv(v["feature_table"], sep="\t")
            taxonomy_table = pd.read_csv(v["taxonomy_table"], sep="\t")
            metadata_table = pd.read_csv(v["metadata_table"], sep="\t")
            model = v["model"]
            model_short_name = DORADO_MODELS[model]["short_name"]
            database = v["database"]

            # -- Add iddentifiers to sample-ids --
            sample_ids = list(metadata_table["sample-id"])
            suffix = f"_{model_short_name}_{database}"

            # Suffixes in metadata and new columns
            metadata_table["sample-id"] = metadata_table["sample-id"] + suffix
            metadata_table["basecalling_model"] = model
            metadata_table["basecalling_model_name"] = model_short_name
            metadata_table["database"] = database
            metadata_table["title"] = (
                "16S rRNA V1-V9 full gene amplicon sequencing of "
                + metadata_table["biosample_name"]
            )
            metadata_table["library_strategy"] = "AMPLICON"
            metadata_table["library_source"] = "METAGENOMIC"
            metadata_table["library_selection"] = "PCR"
            metadata_table["library_layout"] = "single"
            metadata_table["platform"] = "Oxford_Nanopore"
            metadata_table["instrument_model"] = "MinION Mk1C"
            metadata_table["design_description"] = (
                '16S rRNA V1-V9 full gene amplicon sequencing {"seq_run":"'
                + metadata_table["sequencing-run"]
                + '"}'
            )
            metadata_table["filetype"] = "fastq"

            # Suffixes in feature table
            feature_table.columns = [
                i + suffix if i in sample_ids else i for i in feature_table.columns
            ]

            # -- Add identifiers to tax --
            # this is in case the different databases have different taxids
            tax_preffix = f"{database}_"
            feature_table["feature_id"] = tax_preffix + feature_table["feature_id"]
            taxonomy_table["feature_id"] = tax_preffix + taxonomy_table["feature_id"]

            # -- Accumulate --
            feat_tables.append(feature_table)
            tax_tables.append(taxonomy_table)
            md_tables.append(metadata_table)

        # -- Merge --
        feature_table = reduce(
            lambda left, right: pd.merge(left, right, on=["feature_id"], how="outer"),
            feat_tables,
        ).fillna(0)
        taxonomy_table = pd.concat(tax_tables).drop_duplicates()
        metadata_table = pd.concat(md_tables)

        tax_columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        taxonomy_table.columns = ["feature_id"] + tax_columns

        # -- Save --
        feature_table.to_csv(output.feature_table, index=None, sep="\t")
        taxonomy_table.to_csv(output.taxonomy_table, index=None, sep="\t")
        metadata_table.to_csv(output.metadata_table, index=None, sep="\t")


rule merge_ont_and_illumina:
    input:
        ont_feature_table = rules.merge_all_emu_tables.output.feature_table,
        ont_taxonomy_table = rules.merge_all_emu_tables.output.taxonomy_table,
        ont_metadata_table = rules.merge_all_emu_tables.output.metadata_table,
        illu_feature_table = ILLUMINA_DATA_INPUT["feature_table"],
        illu_taxonomy_table = ILLUMINA_DATA_INPUT["taxonomy_table"],
        illu_metadata_table = ILLUMINA_DATA_INPUT["metadata_table"],
    output:
        folder = directory(Path(MERGE_ONT_AND_ILLUMINA)),
        feature_table = Path(MERGE_ONT_AND_ILLUMINA, "feature_table.tsv"),
        taxonomy_table = Path(MERGE_ONT_AND_ILLUMINA, "taxonomy_table.tsv"),
        metadata_table = Path(MERGE_ONT_AND_ILLUMINA, "metadata_table.tsv"),
    run:
        levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

        # -- Read data --
        illu_feature_table = pd.read_csv(input.illu_feature_table, sep="\t", skiprows=1)
        illu_taxonomy_table = pd.read_csv(input.illu_taxonomy_table, sep="\t")
        illu_metadata_table = pd.read_csv(input.illu_metadata_table, sep="\t")

        ont_feature_table = pd.read_csv(input.ont_feature_table, sep="\t")
        ont_taxonomy_table = pd.read_csv(input.ont_taxonomy_table, sep="\t")
        ont_metadata_table = pd.read_csv(input.ont_metadata_table, sep="\t")

        # -- Format ONT --
        ont_taxonomy_table[levels] = ont_taxonomy_table[levels].apply(
            lambda col: col.str.replace("_", " "), axis=1
        )

        # -- Format Illumina --
        illu_feature_table = illu_feature_table.drop("taxonomy", axis=1).rename(
            columns={"#OTU ID": "feature_id"}
        )

        illu_taxonomy_table.columns = [
            i.lower().replace(" ", "_") for i in illu_taxonomy_table.columns
        ]
        illu_taxonomy_table[levels] = (
            illu_taxonomy_table["taxon"]
            .str.replace(".__", "", regex=True)
            .str.replace("_", " ")
            .str.split("; ", expand=True)
        )
        illu_taxonomy_table = illu_taxonomy_table.drop(["taxon", "confidence"], axis=1)
        illu_taxonomy_table[levels] = illu_taxonomy_table[levels].apply(
            propagate_levels, axis=1
        )

        illu_metadata_table["database"] = "silva"

        # -- Merge --
        merged_feature_table = pd.merge(
            ont_feature_table, illu_feature_table, on=["feature_id"], how="outer"
        ).fillna(0)
        merged_taxonomy_table = pd.concat([ont_taxonomy_table, illu_taxonomy_table])
        merged_metadata_table = pd.concat([ont_metadata_table, illu_metadata_table])

        # -- Save --
        merged_feature_table.to_csv(output.feature_table, index=None, sep="\t")
        merged_taxonomy_table.to_csv(output.taxonomy_table, index=None, sep="\t")
        merged_metadata_table.to_csv(output.metadata_table, index=None, sep="\t")

rule biom_conversion:
    input:
        feature_table = rules.merge_ont_and_illumina.output.feature_table,
        taxonomy_table = rules.merge_ont_and_illumina.output.taxonomy_table,
        metadata_table = rules.merge_ont_and_illumina.output.metadata_table,
    output:
        folder = directory(Path(BIOM_CONVERSION)),
        biom_feature_table = Path(BIOM_CONVERSION, "feature_table.biom"),
        json_feature_table = Path(BIOM_CONVERSION, "feature_table.biom.json"),
        taxonomy_table = Path(BIOM_CONVERSION, "taxonomy_table.tsv"),
        metadata_table = Path(BIOM_CONVERSION, "metadata_table.tsv"),
    conda:
        CONDA_ENVS["qiime2"]
    shell:
        """    
        # Transform feature table so phyloseq can read it
        biom convert -i {input.feature_table} -o {output.biom_feature_table} --to-hdf5
        biom convert -i {output.biom_feature_table} -o {output.json_feature_table} --table-type="OTU table" --to-json

        # Simply copy taxonomy and metadata into output folder in order to have everything accessible
        cp {input.taxonomy_table} {input.metadata_table} {output.folder}

        """