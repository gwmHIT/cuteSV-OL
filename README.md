# cuteSV-OL

## Installation

```
git clone https://github.com/120L022331/cuteSV-OL.git && cd cuteSV-OL
```

## Introduction

cuteSV-OL is a novel framework designed for real-time SV discovery, which can be embedded within nanopore sequencing instruments to analyze data concurrently with its generation.

## Dependence

```
scipy	1.10.1
pysam	0.22.1
Cigar	0.1.3
numpy	1.24.4
Biopython	1.83
pyvcf3	1.0.3
scikit-learn	1.3.2
Cython	3.0.11
minimap2	2.28
samtools	1.21
gzip	1.5
```

## Usage

```
cuteSV_ONLINE <monitored_dir> <reference.fa> <work_dir> <output_vcf_dir> 
```

| Optional Parameter | Description                                                  | Default |
| ------------------ | ------------------------------------------------------------ | ------- |
| threads            | Number of threads to use.                                    | 8       |
| mmi_path           | The path of index of reference used in minimap2 to accelerate alignment. | NULL    |
| monitor_fade       | Monitor will close if no new files are detected after monitor_fade second. | 600     |
| high_freq_file     | The path of high frequence SV file or user-defined recall set[vcf]. | NULL    |
| sv_freq            | Specify a high frequency variation threshold for the population to detect. | NULL    |
| user_defined       | The recall set[vcf] is user-defined.                         | False   |
| pctsize            | Min pct allele size similarity between high_freq_file SV set and call set. | 0,9     |
| ref_dist           | Max reference location distance between high_freq_file SV set and call set. | 1000    |
| target_rate        | Stop sequency if the detected rate is higher than target_rate. | 100     |
| batch_interval     | Real-time results are generated every batch_interval batches. | 4       |

**An example**

```bash
export MONITORED_DIR=~/data/experiment/monitor_dir/
export REFPATH=~/data/hg38/hg38.fa
export WORK_DIR=~/data/experiment/work_dir/
export OUTPUTVCF=~/data/experiment/output_vcf/
export CONDAENV=online

# basic usage, and you can get real-time vcf file in OUTPUTVCF directory.
conda activate CONDAENV
cuteSV_ONLINE $MONITORED_DIR $REFPATH $WORK_DIR $OUTPUTVCF

# full usage. Use a population SV file as target recall set.
export THREADS=16
export MMI-PATH=~/data/hg38/hg38_ref.mmi 
export MONITOR_FADE=300 
export POP_FILE=~/data/HGSVC/GRCH38_HGSVC2024v1.0_insdel.vcf # a population SV file, you can defined target recall set in vcf format.
export SV_FREQ=0.1 # Specify a high frequency variation threshold for the population to detect. Don't use if use a self-defined target fille.
export PCTSIZE=0.9
export REF_DIST=1000
export TARGET_rATE=25
export BATCH_INTERVAL=4

conda activate CONDAENV
cuteSV_ONLINE $MONITORED_DIR $REFPATH $WORK_DIR $OUTPUTVCF --mmi_path $MMI-PATH --threads $THREADS --monitor_fade $MONITOR_FADE --high_freq_file $POP_FILE --sv_freq $SV_FREQ --pctsize $PCTSIZE --ref_dist $REF_DIST --target_rate $TARGET_rATE --batch_interval $BATCH_INTERVAL

# full usage. Use a user-defined SV file as target recall set.
export DEFINED_FILE=~/data/experiment/self_defined.vcf

conda activate CONDAENV
cuteSV_ONLINE $MONITORED_DIR $REFPATH $WORK_DIR $OUTPUTVCF --mmi_path $MMI-PATH --threads $THREADS --monitor_fade $MONITOR_FADE --high_freq_file $DEFINED_FILE --user_defined --pctsize $PCTSIZE --ref_dist $REF_DIST --target_rate $TARGET_rATE --batch_interval $BATCH_INTERVAL
```

**output result:**

```
1.vcf_file:In <output_vcf_dir>, you can get real-time result in vcf format, and it also retain old result. File name will indicate its sequence depth.
2.Recall file : The recall result between high_freq_file SV set and call set. Its path is <work_dir>/recall_file.txt.
```
