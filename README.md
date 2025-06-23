
# ChIP-seq analysis pipeline

## description

An analysis pipeline for paired-end ChIP-seq data with the following major steps:

- alignment with [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) (.bam)
- indexing and sorting alignment files with [`samtools`](http://www.htslib.org/) (.bam and .bai)
- summary of fragment sizes with [`deeptools PEFragmentSize`](https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html)
- counting reads aligning to experimental (*S. cerevisiae*) and spike-in (*S.pombe*) genomes with [`samtools view`](http://www.htslib.org/doc/samtools-view.html)
- calculation of per-library spike-in normalization scaling factors using a custom python script
- generation of coverage tracks (.bw) scaled by spike-in normalization using [`deeptools bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
- generation of matrices for calculating correlation between datasets using [`multiBigwigSummary`](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html)
- calculation and plotting of correlation between datasets using [`plotCorrelation`](https://deeptools.readthedocs.io/en/stable/content/tools/plotCorrelation.html)
- enrichment of IP over input coverage (.bw) using [`deeptools bigwigCompare`](https://deeptools.readthedocs.io/en/develop/content/tools/bigwigCompare.html)
- generation of matrices for plotting scaled coverage data using deeptools (.gz) or directly in python (.tab) using [`deeptools computeMatrix`](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#reference-point)
- data visualization ([heatmaps](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html) and [metagenes](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html)) using `deeptools` plotting functions and custom python scripts

## requirements

- Designed to be run on O2, the high-performance computing cluster at HMS
- Uses slurm job scheduler to batch submit jobs
- Paired-end FASTQ files from ChIP-seq libaries
	- FASTQ files should be demultiplexed and can (should) be compressed (.gz)
	- FASTQ filenames are used by the scripts throughout this pipeline, and should easily identify the sample. For example, my filenames usually take the format: strain_treatment_IP_replicate (e.g. 93_D_8WG16_rep2). After demultiplexing, there is usually trailing info added to the filename (e.g. _S8_R1_001.fastq.gz). The S indicates the index number on your sample sheet, R1 and R2 indicate the paired reads from the sequencer, and 001 is a trailing number added for reasons beyond my comprehension.
- Bowtie2 index files (.bt2) for genome alignment, included for *S. cerevisiae* and *S. pombe* in the `genomes/bowtie2_index/` directory in this repository.
- BED files to tell deeptools which portions of the genome you would like to plot. I have included the standard non-overlapping ORF BED file from James Chuang in the `genomes/annotations/` directory of this repository.
- Patience. It will likely take some troubleshooting and path editing in the slurm scripts to analyze your data. My goal is to make this process as pain-free as possible, so I will try to explain how each step functions so that you can debug with confidence.

## instructions

**1. Clone this repository into a clean directory on O2.**

```bash
# clone the repository
git clone https://github.com/winston-lab/warner_chip_seq.git

# navigate to the newly created directory
cd warner_chip_seq
```

Run all commands from the base `warner_chip_seq/` directory.

All Slurm submission scripts have a header that tells Slurm the parameters of the job.  It looks like this:

```
#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-03:00                    # Runtime in D-HH:MM format
#SBATCH --mem=2GB                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent
```

This job uses the 'short' partition, requests 4 cores, has a maximum runtime of 3 hours, requests 2 GB in memory, writes two additonal files (<JOB_ID>.out and <JOB_ID>.err) in the base directory that are useful for troubleshooting if jobs fail, and sends update emails on the jobs to the email address provided. I highly suggest you put your email in here so that you get pinged when jobs complete.

To edit:

```
vim scripts/<script_name_here.sh>

# press 'a' to enter -INSERT- mode
# navigate to the area to edit and make your edits
# press 'esc' to exit -INSERT- mode
# type ':wq' (you should see this in the bottom left) to save (w) and quit (q)

# if you mistype something and want to quit without saving, type ':q!'
```

For more information on Slurm, [see here](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1586793632/Using+Slurm+Basic).

**2. Prepare virtual environments.**

The code written for these analyses relies on two virtual environments `deeptools` and `spike_in` that you will have to set up using `virtualenv` and `pip`. The following commands should be run line-by-line.
(More on virtual environments on O2 [here](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1588662166/Personal+Python+Packages).)

```bash
# start an interactive session
srun --pty -p interactive -t 0-1:00 --mem=1G bash
```

To create the `deeptools` environment:
```bash
# load modules
module load gcc/9.2.0

module load python/3.10.11

# create a local virtual environment
virtualenv env/deeptools --system-site-packages

# activate the virtual environment
source env/deeptools/bin/activate

# install packages to the virtual environment
pip3 install deepTools

#deactivate the environment
source deactivate
```

To create the `spike_in` environment:
```bash
# create a local virtual environment
virtualenv env/spike_in --system-site-packages

# activate the virtual environment
source env/spike_in/bin/activate

# install packages to the virtual environment
pip3 install numpy

pip3 install pandas

pip3 install matplotlib

#deactivate the environment
source deactivate
```

There should now be two folders, one for each environment, in the `env/` directory.

**3. Download or transfer your .fastq.gz files into the `fastq/` directory.**


**4. Align your libraries to the experimental and spike-in genomes.**

We will submit two alignment jobs for each library: one to the experimental genome and the othe to the spike-in genome. Submitting the jobs separately allows all of the alignments to run in parallel.

You will have to edit the `dir` variable in `scripts/batch_aligner_new.sh` and `scripts/spike_batch_aligner_new.sh` to point towards the directory containing your fastq files.

```bash
# use a for loop to submit alignments for each set of paired reads separately
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/batch_aligner_new.sh $name; done
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/spike_batch_aligner_new.sh $name; done
```

Alignment is done using [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) using the following options:
- `--sensitive`
- `--no-unal`
- `--no-mixed`
- `--no-discordant`

This will generate two sets (experimental and spike-in) of three files for each library:
- in `bam/`
	- filename_unsorted.bam
	- filename_sorted.bam
	- filename_sorted.bam.bai
- in `bam/spike-in/`
	- filename_spikein_unsorted.bam
	- filename_spikein_sorted.bam
	- filename_spikein_sorted.bam.bai

```bash
# check to see if the alignment files are there
ls bam/
ls bam/spike-in
```

Also generated is a summary of each alignment in the `logs/` directory:
- filename_bowtie2.txt

```bash
# take a peek
vim logs/<FILE_NAME>_bowtie2.txt
```

We will use the 'sorted.bam' and 'sorted.bam.bai' files in subsequent steps. The 'unsorted.bam' files can be safely deleted.
```bash
rm bam/*unsorted*
rm bam/spike_in/*unsorted*
```

**5. Determine the distribution of insert sizes in your ChIP samples.**

Since the reads are paired, we can determine the size of each fragment that was sequenced from its two ends.

```bash
# use deeptools to generate summary statistics for each sorted.bam file
sbatch scripts/PEFragmentSize.sh
```

This script will look at all of the 'sorted.bam' files and will generate two new files in the `fragment_sizes/` directory:
- a histogram showing the distribution of the insert sizes for each library
- a CSV file with this information in tabular format. This is usually easier to interpret.

These files will be generated for the experimental and spike_in alignments separately.

**6. Count aligned reads for experimental and spike-in genomes.**

This step uses [samtools view](http://www.htslib.org/doc/samtools-view.html) to count the number of reads in each 'sorted.bam' file. 
- `-c` makes `samtools` output only the count and not the reads themselves
- `-F` filters the BAM files to EXCLUDE reads that fit the [flag condition](https://broadinstitute.github.io/picard/explain-flags.html)
	- the flag here is `388` which = unmapped OR not primary alignment OR second in pair (which confirms read pairs are only counted once and not twice)

What the script looks like under the hood:

> ```bash
> for name in bam/*_sorted.bam; do
> 	(basename ${name} _sorted.bam) >> logs/experimental_counts.log
> 	samtools view -c -F 388 ${name} >> logs/experimental_counts.log
> done
> ```

To run:

```bash
sbatch scripts/read_counter.sh
```

This scripts creates two files in the `logs/` directory: 'experimental_counts.log' and 'spikein_counts.log'. These files contain pairs of lines:
- the first line is the BAM filename with '_sorted.bam' stripped
- the second line is the number of paired reads that are mapped to the respective genome in each 'sorted.bam' file
The lines continue to alternate for each BAM file processed.


**7. Do spike-in normalization math.**

Spike-in normalization math is not intuitive (at least it isn't to me). I have included a PDF ('chip_spikeins.pdf') in this repository that was written by James Chuang and explains all of the logic and algebra behind this step. In short, the 'input' libraries allow us to empirically determine the proportion of experimental to spike-in material that went into each IP. 

Essentially, if we wanted to normalize the experimental input signal between libraries, we could just normalize by library size (number of experimental reads). Scaling by library size does not work for the IPs, however, as getting half as many experimental reads is meaningful (as long as the number of spike-in reads is the same). Therefore, if we use the spike-in reads to "link" the IP and input samples, we can scale each IP to the same scale as the input signal. Then you can calculate the IP enrichment over input, which is exactly what we want to do! (And is what we do in step 8.)

For the python script in this step to work for any data sets other than the ones used in my ChIP experiment, you will probably need to change indexing the for loop in lines 97-98. I haven't been able to figure out how to make the indexing step run autonomously, and so you need to figure out the indexing math to line up your IP samples with the appropriate inputs.

If you look carefully at the math that is being done, I multiply the calculated scaling factor by 10000000. This is purely to make the resulting coverage numbers human readable: since each scaling factor is multiplied by the same constant, it does not affect the ratios between libraries.

Run the following commands line-by-line:

```bash
# start an interactive session (if not already interactive)
srun --pty -p interactive -t 0-1:00 --mem=1G bash

# load modules
module load gcc/9.2.0

module load python/3.10.11

# activate the virtual environment
source env/spike_in/bin/activate

# run the python script
python scripts/mnase_spikein_norm.py

# deactivate the environment
deactivate
```

(More on virtual environments on O2 [here](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1588662166/Personal+Python+Packages).)

The output of all of this is a file called 'normalization_table.csv' in `logs/` that consists of two columns:
- column 1 is the library name
- column 2 is the scaling factor that will be used for normalization

```bash
# take a peek at the file
vim logs/normalization_table.csv
```

This step also generates a plot of the proportion of each library that aligned to either the *S. pombe* or *S. cerevisiae* genome. It's called 'proportion_reads_mapped.png' and can also be found in `logs/`.


**8. Generate coverage tracks for each library scaled by spike-in normalization.**

Before we can calculate IP enrichment over input, we have to generate coverage tracks. We can use the scaling factors that we calculated in the previous step to normalize the coverage tracks at this stage.

Similarly to how we submitted parallel jobs for each alignment, we'll do the same here. But how do we know which scaling factor to use? We need to match the BAM filename to the scaling factor in the 'normalization_table.csv' file. The first few lines of code handle this:

> ```bash
> base=$(basename ${1%} _sorted.bam)
> 
> alpha=$(grep ${base} logs/normalization_table.csv | cut -f2 -d,)
> ```

The first line strips the '_sorted.bam' from the filename and stores is as the variable 'base'. We did the same when we counted the reads. This means the filenames should match.

The second line defines 'alpha' which is the scaling factor. It uses [`grep`](https://man7.org/linux/man-pages/man1/grep.1.html) to look in 'normalization_table.csv' for lines that match 'base'. There should only be one line that matches. `|` pipes that line into [`cut`](https://man7.org/linux/man-pages/man1/cut.1.html), which prints everything past the delimiter, here defined as ','. Since the second column (meaning everything past the ',') on that line is the matched scaling factor, this gets stored as 'alpha'.

The rest of the script just plugs 'alpha' in as the scaling factor to `bamCoverage`. There are some parameters you can change here, read the [documentation](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) if you feel like tweaking any of them. The parameters here are what I used for my Spt6 and Rpb1 ChIP-seq.

> ```bash
> bamCoverage -b ${1%} -o deeptools/si/${base}_si.bw \
>        -bs 20 \
>        --scaleFactor ${alpha} \
>        --extendReads \
>        --centerReads \
>        -p max \
>        --smoothLength 60 \
>        --ignoreForNormalization chrM \
>```

The `--ignoreForNormalization` is probably superfluous since you are explicitly passing a scaling factor, but I have left it in for now (since I was previously using library size normalization for MNase-seq). I may remove it and this text at a later date.

```bash
# use a for loop to submit each coverage job in parallel
for name in bam/*_sorted.bam; do sbatch scripts/bamCoverage_si.sh $name; done
```

This step outputs BIGWIG (.bw) files.

```bash
# check to see the files were made
ls deeptools/si/
```


**9. Calculate and plot correlation scores.**

```bash
sbatch scripts/multiBigwigSummary_Correlation.sh
```

This script calculates coverage scores in non-overlapping 200 bp bins genome-wide (skipping the mitochondrial genome).

The binned coverage scores are used to calculate Spearman correlation coeeficients between all libraries. The results are plotted as a heatmap in the `correlation/plots` directory. It also outputs the results in tabular (.tab) format in the `correlation/tab` directory if you would like to plot it yourself. Sample code to plot for yourself in python is provided as `scripts/plot_correlation.py`.

**10. Average biological replicates to generate coverage tracks.**

If the correlation between biological replicates looks good, you can average your coverage tracks before plotting your results.

```bash
# calculate average coverage across all replicates in 10 bp windows
sbatch scripts/average_bigwigAverage.sh
```

The averaged files (.bw) can be found in `deeptools/averaged/`. These files can be loaded by IGV to visualize nucleosome dyad coverage at individual loci.

**11. Calculate enrichment of IP/input coverage.**

This step calculates the ratio of IP to input for each bin.

```bash
sbatch scripts/average_bigwigCompare_chip.sh
```

There should now be new .bw files in the `deeptools/averaged/ratio` directory with `<IP>vinput_averaged_midpoints_si_ratio.bw` style names. These show the enrichment of IP coverage over input coverage.
Also generated are .bw files in the `deeptools/averaged/perpol` directory with `<IP>v8WG16_averaged_midpoints_si_ratio.bw` style names. These show the Rpb1-normalized enrichment of IP coverage.

**12. Generate matrices to plot data.**

Using these new enrichment files, generate matrices with which to plot the data.

```bash
# Sorted by Rpb1 ChIP-seq occupancy in wild-type, with all genes are scaled to the same length. Plotting data is generated for 250 bp upstream of the TSS and 250 bp downstream of the CPS.
sbatch average_computeMatrix_midpoints_scale.sh

# Sorted by gene length. Plotting data is generated for 250 bp upstream and up to 4500 bp downstream. Positions after the end of the gene are filled with nan.
sbatch average_computeMatrix_midpoitns_reference.sh
```

There should now be new .tab files in the `deeptools/averaged/<ratio OR perpol>/tab` directory that will be used in the next step to plot the results. The compressed matrices are also saved in `deeptools/averaged/<ratio OR perpol>/gz` if you would like to use deeptools to generate plots.


**13. Import and plot data.**

The .tab files can be used to plot the results in python. I do this by moving the files to my local machine, though you could certainly also try to do the plotting on O2. Sample codes to plot for yourself in python are provided as `scripts/plot_chip_<ETC>.py`.

To transfer .tab files to your machine:
```bash
# from a local session
scp -r <YOUR_USER_ID>@transfer.rc.hms.harvard.edu <PATH/TO/DIR/>deeptools/averaged/ratio/tab <./PATH/TO/LOCAL/DIR>
scp -r <YOUR_USER_ID>@transfer.rc.hms.harvard.edu <PATH/TO/DIR/>deeptools/averaged/perpol/tab <./PATH/TO/LOCAL/DIR>
```
To fix the first few lines of the files and make them python-readable, run the following code (provided as `scripts/tab_converter.sh` on the directory containing the .tab files:
```bash
for name in <FILE/PTAH/>*.tab; do sed '3s/genes:3087\t/#genes:\n/g' $name > lala.tab && mv lala.tab $name; done
```


