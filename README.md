IsoPrimer is a pipeline built to automate the designing process of PCR primer pairs, targeting specific sets of expressed splicing variants.
IsoPrimer works with any annotation in the gtf format and allows the prioritization of the designed primers pairs, according to the level of expression of the splicing variants of a gene in an RNA-seq dataset.
The pipeline leverages Kallisto, Primer3 and EMBOSS PrimerSearch tools, respectively to:
1. identify the most expressed isoforms per gene in RNA-seq data;
1. design and consider all primer pairs overlapping exon-exon junctions common to the expressed variants;
1. verify the specificity of the designed primer pairs.

## Requirements

IsoPrimer was designed to run on GNU+Linux systems and it requires the following programs:

-   [R](https://www.r-project.org/) (at least 3.6.x)
-   [Bash](https://www.gnu.org/software/bash/) (at least 4.x.x)
-   [Kallisto](https://www.nature.com/articles/nbt.3519)
-   [Primer3](https://github.com/primer3-org/primer3)
-   [EMBOSS PrimerSearch](https://www.bioinformatics.nl/cgi-bin/emboss/help/primersearch)

Furthermore, the following R packages (in the versions indicated or newer) are
necessary to run the pipeline. 

-   dplyr 0.8.3
-   BiocManager 1.30.10
-   Biostrings 2.54.0
-   doParallel 1.0.16
-   openxlsx 4.1.5
-   stringr 1.4.0
-   msa 1.18.0

R, Kallisto, Primer3, PrimerSearch and the necessary R packages are
available as part of the respective **conda** packages. 
If necessary, to install conda, it is important to have a compatible Python version
and we suggest to refer to conda's [installation guide](https://www.anaconda.com/docs/getting-started/miniconda/install#macos-linux-installation).
If `conda` is available on the system running IsoPrimer, the dependency
installation process via conda may be automated by using the helper script `install.sh`
bundled with IsoPrimer:

```bash
bash install.sh
conda activate isoprimer
```

Or you may install the IsoPrimer dependencies through conda manually by creating
an ad-hoc environment (changing `<my-env>` with the desired name of the environment) as:

```bash
conda create --name <my-env>
```

After the activation of the environment with

```bash
conda activate <my-env>
```

the necessary dependencies may be installed via the following commands:

```bash
conda install conda-forge::r-base=4.3.3
conda install bioconda::emboss
conda install bioconda::primer3
conda install bioconda::kallisto=0.46.2
conda install conda-forge::r-dplyr
conda install conda-forge::r-biocmanager
conda install bioconda::bioconductor-biostrings
conda install conda-forge::r-doparallel
conda install conda-forge::r-openxlsx
conda install conda-forge::r-stringr
conda install bioconda::bioconductor-msa
```

Upon successful loading of the dependencies, a dedicated message will be printed both
to stdout (see the `nohup.out` file) and the `IP_Log.out` file.

### Alternatively

Pre-compiled versions of Kallisto are available
[here](https://pachterlab.github.io/kallisto/download) or the program
may be compiled from source following the instructions available at the
author's [github page](https://github.com/pachterlab/Kallisto/blob/master/INSTALL.md).

Primer3 may be compiled from source by following the instructions
available at the dedicated [github page](https://github.com/primer3-org/primer3).

The latest version of EMBOSS PrimerSearch may be obtained by first
downloading source code:

```bash
wget ftp://emboss.open-bio.org:21/pub/EMBOSS/EMBOSS-6.6.0.tar.gz
```

and then by compiling it as described in the dedicated
[page](https://emboss.sourceforge.net/download/).

IsoPrimer will automatically download and install the required R packages if necessary.

## Setting up the pipeline

### Pipeline setup outline

The outline of the pipeline setup is as follows:

Clone the IsoPrimer repository:

```bash
git clone https://github.com/ermesfil/IsoPrimer.git
```

Then move to the IsoPrimer folder and launch the pipeline without modifying any files to run a test:

```bash
cd IsoPrimer
bash launcher.sh
```

The progress of IsoPrimer may be monitored via the `nohup.out` (stderr and
stdout) and `IP_Log.out` (IsoPrimer messages) files.
If the test run completed successfully all necessary dependencies are correctly installed
and the message 'Primer design complete' will be printed in the `IP_Log.out` file.

You may now remove the test specific files by running the following in the IsoPrimer folder:

```bash
bash purger.sh
```

Finally, modify the `targets.txt`, `quantification/sample_list.txt` and `launcher.sh` text files appropriately
as described in the following paragraphs then start the primer design with:

```bash
bash launcher.sh
```

### Preparing the transcriptome headers

To ensure consistent processing of every transcript of a given target gene,
use an annotation file and a transcriptome FASTA whose sequence headers match
the transcript identifiers in that annotation.

The headers of the transcriptome FASTA file must conform to the
following format:

```bash
>transcript-id
```

As an example:

```bash
>ENST00000456328.2
```

IsoPrimer will prepare a complying copy of the transcriptome FASTA file:
`mod_transcriptome.fa` which will be used for the subsequent steps of the
IsoPrimer pipeline.

### Setting up Kallisto

Kallisto requires an index to be created to pseudoquantify transcripts.
You may adjust the parameters in the indexing section of the
quantification/kallister.sh script if necessary. See the [Kallisto manual](https://pachterlab.github.io/kallisto/manual) for more details.

Kallisto requires information about the RNA-seq library to quantify its
transcripts. Specifically, if the sequencing is stranded and paired-ended, its strandedness should be
specified via the dedicated parameters `--rf-stranded` or `--fr-stranded`.
A single-ended library requires the `--single`, `-s` and `-l` parameters instead. Please
refer to the [Kallisto manual](https://pachterlab.github.io/kallisto/manual)
for a complete description of the options available.

To set up Kallisto, open the `quantification/kallister.sh` script and
change the parameters as necessary (lines 15-23). Additional
quantification options may be specified here. More information is
provided in the script comments.

By default, IsoPrimer launches the Kallisto quantification process in
paired-end mode and with the `--rf-stranded` option.

### Specifying input files

Use targets.txt to specify what target genes to design primer pairs for.
The file should contain a list of the gene IDs. E.g.

```bash
ENSG00000224609
ENSG00000227811
ENSG00000225937
```

if using the ENSEMBL annotation or 

```bash
FGGY-DT
INKA2-AS1
PCA3
```

for the ncbi annotation.

Specify the paths to the RNA-seq sample **folders** in
`quantification/sample_list.txt`. E.g.

```bash
/home/username/sample1
/home/username/sample2
/home/username/samplen
```

For paired end libraries, each folder should contain two files: a `*1.fastq.gz` for the forward reads and a
`*2.fastq.gz` for the reverse. For single end designs, each sample folder
should contain the relative fastq file.

Please modify the `launcher.sh` script to pass the options necessary to
run the pipeline as follows:

```bash
nohup Rscript IP.R \
    number of threads \
    path to the Kallisto executable \
    path to the primer3_core executable \
    path to the EMBOSS PrimerSearch executable \
    Isoform significance threshold percentile [0,100 integer] \
    PrimerSearch mismatch percentage allowed [0,100 integer] \
    path to the transcriptome \
    path to the annotation \
    path to the genome \
    debug flag (T/F) &
```

If the Kallisto, Primer3 and Primersearch executables are available via
the `PATH`, environment variable (which is the case if a **dedicated conda
environment** was created and activated as mentioned above), only the transcriptome,
annotation and genome fields may be modified.

Please remove any metacharacters and/or trailing white spaces from the
script to avoid runtime errors. Remove `nohup` and `&` to launch the command
in foreground.

## Launching the pipeline

Start the primer design:

```bash
bash launcher.sh
```

## Results

The main final product of the primer design process consists in a primers.xlsx file with
3 spreadsheets: I) validation_candidates which contains a list of the primer pairs that satisfy the predefined quality criteria;
II) primers_order which reports primer pairs in a ready-for-order format;
III) primer_omnibus which reports all primer pairs designed for all target genes.

The output spreadsheet I) and III) report for each primer pair the
following information:

-   gene ID
-   IsoPrimer score: the higher, the better
-   transcript ID
-   Forward primer sequence
-   Reverse primer sequence
-   *Expressed amplified*, the number of isoforms whose expression is
    above the percentile threshold that are amplified by the primer pair
-   *Expressed percentage*, the percentage of the gene's total
    expression, quantified in TPM, that is amplified by the primer pair
-   The amplicons predicted by the PrimerSearch in-silico PCR expressed
    as a space separated list of
    `<transcript_id-amplicon_length>` entries

IsoPrimer was designed to address all the expressed splicing variants,
but this may not be possible with a single primer pair if the variants
of interest do not share a common splicing junction. For this reason,
IsoPrimer may return more than one primer couple for each target gene.
However all primer pairs must respect the quality criteria defined.

Finally, for each target gene, a non tabular detailed report is
generated in the *outputs* folder with the following information:

-   quantitation of the expression of the isoforms of a gene obtained
    with Kallisto and most expressed isoforms (isoforms more expressed
    than the average reported by Kallisto);
-   primer sequences;
-   isoform used to design the primer pair and splicing junction(s)
    overlapped by primers;
-   representation of the hybridization of the primers with the matching
    multialigned isoforms of a gene;
-   FASTA sequences of the isoforms matched by the primers designed;
-   relevant information returned by Primer3 (e.g. putative secondary
    structures and hairpins generated within the primers and relative
    estimates about their stability)

## Primer design for manually picked isoforms

IsoPrimer allows for designing primer pairs targeting known variants even without RNA-seq data.
In fact, given the appropriate annotation files and a customized
quantification output specified in the `IsoPrimer/quantification/KA_CountingOutput/kalcounts.tsv`,
prepared by manually entering high expression values for one or more isoforms of interest,
IsoPrimer will provide primer pairs designed to amplify the manually prioritized isoforms.
For example, to instruct IsoPrimer to specifically design a primer pair for the *Odc1-207* transcript variant
(ENSEMBL ID ENSMUST00000222250) of the murine *Odc* gene, without performing the transcript quantification,
the user may provide a `kalcounts.tsv` file structured as exemplified in the table below.

| t_name            | tseq  |
|-------------------|-------|
| ENSMUST00000171737 | 0     |
| ENSMUST00000222617 | 0     |
| ENSMUST00000221701 | 0     |
| ENSMUST00000220947 | 0     |
| ENSMUST00000221613 | 0     |
| ENSMUST00000220849 | 0     |
| ENSMUST00000221354 | 0     |
| ENSMUST00000222250 | 1e+06 |

## Troubleshooting

IsoPrimer progressively saves the results it produces to facilitate the
debug process and to restart from the last step it performed upon
interruption.

When changing annotation, launch `purger.sh` in the IsoPrimer folder
to delete the previous results and reset the pipeline (except for the
xlsx file which gets overwritten at each primer pair design run).
