# rRNA_pipeline
#### Pipeline for read filtering, swarm OTU clustering, and taxonomic classification of ssu-rRNA

![rrna_pipeline_procedure](https://cloud.githubusercontent.com/assets/14023091/15032729/e97deba4-1218-11e6-870f-1715e8b116a0.jpg)

Usage
-----

| File | Description |
|------|-------------|
| rRNA_pipeline.py | Full rRNA pipeline |
|  |  |
| init.txt | optional paths to alternate databases |
| [db/](./db/) | ssu-rRNA databases |
| fastq_filter.py | FASTQ filtering |
| fasta_dereplicate.py | FASTA dereplication |
| swarm_map.py | run swarm |
| swarm_classify_taxonomy.py | classify swarm OTUs |
| group_taxa.py | aggregate taxonomic group counts
| purity_plot.py | classifify OTU content |
| plot_OTU_purity.r | produce OTU purity plot |
| plot_sample_correlations.r | produce sample tree plot |
| plot_taxa_groups.r | produce taxonomic groups plot |
| plot_diversity.r | produce sample diversity plot |
| plot_heatmap.r | produce sample/OTU heatmap |

```
rRNA_pipeline v0.3 (May 19, 2016)
Full ssu-rRNA, swarm OTU classification pipeline

Usage: rRNA_pipeline.py (options)
   -d name          : database name (16S, V4, V9)
   -q dir           : FASTQ folder
   -o file          : base filename for results (default: rrna)
   -n file          : sample names file (optional)
   -p               : calculate/plot OTU purity
   -m int           : minimum quality score for FASTQ (default: 35)
   -s, --steps list : run only the steps in list (default: All)
   -t, --cpus int   : number of processes (default: 1)
   -w               : no overwrite of files, skip completed steps (default)
   -W, --overwrite  : overwrite files (default if -s)
   -h, --help       : help
   -v, --verbose    : more information to stderr

   If -s is used, then -W overwrite is the default
   Steps list can include any of the following: 
      All, merge_fastq, chimera, filter_fasta, derep
      swarm, classify, plots, purity, split_plastid
      plastid_classify, plastid_plots, plastid_purity

Example: rRNA_pipeline.py -d 16S -q ./fastq
```

If using alternate databases for 16S, Plastid, 18S_V4, or 18S_V9, specify paths in init.txt.

**Basic operation for 16S:**
```bash
rRNA_pipeline.py -d 16S -q ./fastq
```

To replace FASTQ filenames with sample names in all output, use -n to specify tab-delimited file (sample_name, FASTQ base name).  FASTQ base names may be followed by any of [_R1, _R2, .filtered, .fastq, .fq] in the full FASTQ file name.  

The basic pipeline runs relatively quickly, however the extra calculation of OTU purity takes much longer.  Use -p to calculate and plot purity.

**Use the following for 18S V4, with sample names, run on 4 CPUs, with purity plot:**
```bash
rRNA_pipeline.py -d V4 -o rrna -n sample_names.txt -t 4 -p
```

The rRNA pipeline will skip previous steps if stopped and rerun, unless -W is specified.  To rerun starting at a particular step, delete the output files created during that step, and run again.  To run only specific steps, use -s and specify steps by name.  Spaces are not allowed unless the list is bound by quotes.  If -s is used then -W is the default and files will be overwritten unless -w is given.

**Use the following for 18S V9, running only derep, swarm, and classify steps:**
```bash
rRNA_pipeline.py -d V9 -s "derep,swarm,classify"
```

*Output files, created in the following order:*

| File | Description |
|------|-------------|
| fqbase1.assembled.fastq | Pear merged paired reads
| fqbase1.discarded.fastq | Pear unmerged reads
| fqbase1.unassembled.forward.fastq | Pear unmerged reads R1 
| fqbase1.unassembled.reverse.fastq | Pear unmerged reads R2
| fqbase1.uchime | Usearch -uchime_ref list of chimeric reads
| fqbase1.filtered.fa | final set of filtered reads
| ... | |
| | |
| rrna.derep.fa | dereplicated reads |
| rrna.derep.counts | read counts for dereplicated reads |
| rrna.swarm | swarm dereplicated reads in each swarm cluster |
| rrna.swarm.fa | representative swarm reads |
| rrna.swarm.counts | swarm OTU sample counts table |
| rrna.swarm.ggsearch | Fasta36 m8 output |
| rrna.swarm.tax | swarm OTU taxonomy and counts table |
| rrna.swarm.sample_corr.pdf | plot of sample correlation tree |
| rrna.swarm.diversity.pdf | plots of sample richness/diversity metrics |
| rrna.swarm.heatmap.pdf | plot of normalized sample OTU abundance |
| rrna.taxa_groups.txt | aggregated counts for taxonomic groups |
| rrna.taxa_groups.pdf | plot of taxonomic group counts |
| rrna.swarm.content.fa | dereplicated reads in largest swarm clusters |
| rrna.swarm.content.ggsearch | Fasta36 m8 output |
| rrna.swarm.content.tax | swarm content taxonomy |
| rrna.swarm.purity | swarm purity table |
| rrna.swarm.purity.pdf | plot of swarm purity |

(multiple paired FASTQ base filename *'fqbase1'*, and single output base filename *'rrna'*)

In the case of 16S, files are split into two groups: those that are classified as plastid, and those that are not (determined by *'rrna.swarm.ggsearch'*).  The non-plastid rows/sequences are written over the base *'rrna'* files, and the plastid 16S rows/sequences are written to a separate set of files, *'rrna.plastid'*, and swarm OTUs are re-classified by the phytoRef database.

Installation
------------

1. Make sure all dependencies are installed (see below), and make them accessible to your path.
2. Databases derived from PR2 with updated taxonomy are included in db/.  Use gunzip to uncompress these files.
3. Download SILVA NR database (http://www.arb-silva.de/), and add to db/ if you want to use 16S
4. Test the rRNA pipeline:

```bash
rRNA_pipeline.py --test
```

Dependencies
------------

* Python 2.7 (https://www.python.org/downloads/)
* R (https://cran.r-project.org/)
* PEAR (https://github.com/xflouris/PEAR.git)
* USEARCH (http://www.drive5.com/usearch/download.html)
* SWARM (https://github.com/torognes/swarm)
* FASTA36 (https://github.com/wrpearson/fasta36)
