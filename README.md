# cuteSV-OL:a real-time structural variation detection framework for nanopore sequencing devices

## Installation

**cuteSV-OL requires miniconda to create the runtime environment**

```bash
git clone https://github.com/gwmHIT/cuteSV-OL.git && cd cuteSV-OL && conda env create -f environment.yml -n <your_env_name> && python setup.py build_ext --inplace && python setup.py install
```

**you can also use conda to install, and be sure python < 3.13**

```bash
conda install -c conda-forge -c bioconda cutesv-ol
```

## Introduction

cuteSV-OL is a novel framework designed for real-time SV discovery, which can be embedded within nanopore sequencing instruments to analyze data concurrently with its generation.

## Dependence

```
python  3.8
scipy   1.10.1
pysam   0.22.1
Cigar   0.1.3
numpy   1.24.4
Biopython   1.83
pyvcf3  1.0.3
scikit-learn    1.3.2
Cython  3.0.11
minimap2    2.28
samtools    1.21
watchdog	4.0.1
```

## Usage

```
cuteSV_ONLINE <monitored_dir> <reference.fa> <work_dir> <output_vcf_dir> 
```

| Optional Parameter | Description                                                  | Default |
| ------------------ | ------------------------------------------------------------ | ------- |
| threads            | Number of threads to use.                                    | 4       |
| mmi_path           | The path of index of reference used in minimap2 to accelerate alignment. | NULL    |
| monitor_fade       | Monitor will close if no new files are detected after monitor_fade second. | 600     |
| target_set         | The path of high frequence SV file or user-defined target recall set[vcf] as the ground truth set. | NULL    |
| sv_freq            | Specify a high frequency variation threshold for the population to detect.It doesn't need if target_set doesn't have the attribute of AF. | 1.0    |
| user_defined       | The target recall set[vcf] is user-defined.                  | False   |
| pctsize            | Min pct allele size similarity between high_freq_file SV set and call set. | 0,9     |
| ref_dist           | Max reference location distance between high_freq_file SV set and call set. | 1000    |
| target_rate        | Stop sequency if the detected rate is higher than target_rate. | 100     |
| batch_interval     | Real-time results are generated every batch_interval batches. | 4       |

### **notice**

If you are using a population SV dataset as ground truth and want to specify a threshold for high frequency variation, you can use the sv_freq parameter and use the following command to generate the AF field for the vcf file:

```bash
bcftools +fill-tags input.vcf -Ov -- -t AF -o output.vcf  # <vcf> format
bcftools +fill-tags input.vcf.gz -Ou -- -t AF | bcftools view -Oz -o output.vcf.gz  # <vcf.gz> format
```

If you are using a custom SV dataset as the ground truth, use the user_defined parameter, which will ensure that all SV will be the target to be detected.

### Dataset

| **Dataset**                         | **Link**                                                     |
| ----------------------------------- | ------------------------------------------------------------ |
| HG002 ONT fastq.gz                  | https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.4.5/HG002_ONT-UL_GIAB_20200204.fastq.gz |
| GRCH38  HG002 T2T data              | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz |
| GRCH38  HGSVC population SV INS&DEL | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_sv_insdel_sym_HGSVC2024v1.0.vcf.gz |
| GRCH38 reference                    | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26 |

We use the following script to bench the output vcf file. it can be a little slow to refine the result. Truvari version is v5.1.1 .

```bash
#!/bin/bash
truvari_dir=$1
vcf_path=$2
source /home/user/guoweimin/miniconda3/etc/profile.d/conda.sh
rm -r /home/user/guoweimin/data/vcf_refine/$truvari_dir
conda activate truvari_git
export PATH="$PATH:/home/user/guoweimin/tools/mafft/bin" # a external tool mafft should be install.
truvari bench \
    --reference /home/user/guoweimin/data/standard_vcf/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --includebed /home/user/guoweimin/data/standard_t2t_vcf/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed \
    --base /home/user/guoweimin/data/standard_t2t_vcf/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz \
    --comp $vcf_path \
    --output /home/user/guoweimin/data/vcf_refine/$truvari_dir \
    --passonly \
    --pick ac \
    --refdist 2000 \
    -C 6000 \
    --dup-to-ins

truvari refine \
    --reference /home/user/guoweimin/data/standard_vcf/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --regions /home/user/guoweimin/data/vcf_refine/$truvari_dir/candidate.refine.bed \
    --use-original-vcfs \
    --align mafft \
    --mafft-params '--auto --thread 32' \
    -t 32 \
    /home/user/guoweimin/data/vcf_refine/$truvari_dir
```

### **An example**

```bash
export MONITORED_DIR=~/data/experiment/monitor_dir/
export REFPATH=~/data/hg38/hg38.fa
export WORK_DIR=~/data/experiment/work_dir/
export OUTPUTVCF=~/data/experiment/output_vcf/
export CONDAENV=online

# basic usage, and you can get real-time vcf files in OUTPUTVCF directory.
conda activate CONDAENV
cuteSV_ONLINE $MONITORED_DIR $REFPATH $WORK_DIR $OUTPUTVCF

# full usage. Use a human common SVs from the HGSVC dataset (Ebert et al. 2021) as the ground truth.
export THREADS=16
export MMI-PATH=~/data/hg38/hg38_ref.mmi 
export MONITOR_FADE=300 
export POP_FILE=~/data/HGSVC/GRCH38_HGSVC2024v1.0_insdel.vcf # a built-in population SV file in src/data, you can also defined target recall set in vcf format.
export SV_FREQ=0.1 # Specify a high frequency variation threshold for the population to detect. Don't use it if use a self-defined target set as the ground truth.
export PCTSIZE=0.9
export REF_DIST=1000
export TARGET_RATE=25
export BATCH_INTERVAL=4

conda activate CONDAENV
cuteSV_ONLINE $MONITORED_DIR $REFPATH $WORK_DIR $OUTPUTVCF --mmi_path $MMI-PATH --threads $THREADS --monitor_fade $MONITOR_FADE --target_set $POP_FILE --sv_freq $SV_FREQ --pctsize $PCTSIZE --ref_dist $REF_DIST --target_rate $TARGET_RATE --batch_interval $BATCH_INTERVAL

# full usage. Use a user-defined SV file as target recall set.
export DEFINED_FILE=~/data/experiment/self_defined.vcf

conda activate CONDAENV
cuteSV_ONLINE $MONITORED_DIR $REFPATH $WORK_DIR $OUTPUTVCF --mmi_path $MMI-PATH --threads $THREADS --monitor_fade $MONITOR_FADE --target_set $DEFINED_FILE --user_defined --pctsize $PCTSIZE --ref_dist $REF_DIST --target_rate $TARGET_RATE --batch_interval $BATCH_INTERVAL
```

**output result:**

```
1.vcf_file:In <output_vcf_dir>, you can get real-time result in vcf format, and it also retain old result. File name will indicate its sequence depth.
2.Recall file : The recall result between target recall set and cuteSV-OL call set. Its path is <work_dir>/recall_file.txt.
```
