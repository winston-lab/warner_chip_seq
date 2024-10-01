
#ChIP-seq analysis pipeline

## description

An analysis pipeline for paired-end ChIP-seq data with the following major steps:

- alignment with bowtie2 (.bam)
- indexing and sorting alignment files (.bam and .bai)
- summary of fragment sizes with deeptools PEFragmentSize
- counting reads aligning to experimental (S. cerevisiae) and spike-in (S.pombe) genomes with samtools
- calculation of per-library spike-in normalization scaling factors using a custom python script
- generation of coverage tracks (.bw) scaled by spike-in normalization using deeptools bamCoverage
- log2 fold enrichment of IP over input coverage (.bw) using deeptools bigwigCompare
- generation of matrices for plotting scaled coverage data using deeptools (.gz) or directly in python (.tab) using deeptools computeMatrix
- data visualization (heatmaps and metagenes) using deeptools plotting functions

## requirements

- designed to be run on O2, the high-performance computing cluster at HMS
- uses slurm job scheduler to batch submit jobs
- Paired-end FASTQ files from ChIP-seq libaries
	- FASTQ files should be demultiplexed and can (should) be compressed (.gz)
	- FASTQ filenames are used by the scripts throughout this pipeline, and should easily identify the sample. For example, my filenames usually take the format: strain_treatment_IP_replicate (e.g. 93_D_8WG16_rep2). After demultiplexing, there is usually trailing info added to the filename: e.g. _S8_R1_001.fastq.gz. The S indicates the index number on your sample sheet, R1 and R2 indicate the paired reads from the sequencer, and 001 is a trailing number added for reasons beyond my comprehension.
- FASTA files for genome alignment, included for S. cerevisiae and S. pombe in the 'genomes/bowtie2_index/' directory in this repository.
- BED files to tell deeptools which portions of the genome you would like to plot. I have included the standard non-overlapping ORF BED file from James Chuang in the 'genomes/annotations/' directory of this repository.
- Some patience. It will likely take some troubleshooting and path editing in the slurm scripts to analyze your data. My goal is to make this process as pain-free as possible, so I will try to explain how each step functions so that you can debug with confidence.

## instructions

**1.** Clone this repository into a clean directory on O2.

```bash
# clone the repository
git clone <INSERT URL HERE>

# navigate to root directory
cd chip_seq
```


**2.** Download or transfer your .fastq.gz files into the 'fastq/' directory.

**3.** Align your libraries to the experimental and spike-in genomes.

Run all commands from the root 'chip_seq' directory.

```bash
# use a for loop to submit alignments for each set of paired reads separately
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/batch_aligner.sh $name; done
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/spike_batch_aligner.sh $name; done
```

This will generate two sets (experimental and spike-in) of three files for each library:
- in bam/
	- filename_unsorted.bam
	- filename_sorted.bam
	- filename_sorted.bam.bai
- in bam/spike-in/
	- filename_spikein_unsorted.bam
	- filename_spikein_sorted.bam
	- filename_spikein_sorted.bam.bai

```bash
# check to see if the alignment files are there
ls bam/
ls bam/spike-in
```

Also generated is a summary of each alignment in the logs/ directory:
	- filename_bowtie2.txt

```bash
# take a peek
vim logs/<FILE_NAME>_bowtie2.txt
```

We will use the 'sorted.bam' and 'sorted.bam.bai' files in subsequent steps. I don't think tha the 'unsorted.bam' files need to be saved, but I have not made a habit of deleting them.


