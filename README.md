Metagenomics
============

The software in this repository is related to the paper titled “**Detecting Bacterial Genomes in a Metagenomic Sample Using NGS Reads**”.  The software here covers the statistical and computational methods mentioned in the above paper, along with miscellaneous programs, scripts, and notes.

## How To Get Started

- [Download the Code](https://github.com/camilo-v/Metagenomics/archive/master.zip) and replace the generic Path variables to match your environment/computer.

## Communication

- If you **found a bug**, open an issue and _please provide detailed steps to reliably reproduce it_.
- If you have **feature request**, open an issue.
- If you **would like to contribute**, please submit a pull request.

## Requirements

The software was developed on OS X (10.9) and covers four (4) programming languages:

- Python (2.7 variant)
- Perl (perl 5, v.16 or above)
- R (3.1 or above)
- Bash Shell

You will need a reasonably modern computer and operating system that can support all four.

### Bowtie2

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is required for the alignment step.  See the Bowtie2 [Getting Started](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example) guide for instructions on installing and running Bowtie2.

### Samtools

[Samtools](http://www.htslib.org/) is required for parsing the alignments and generating miscellaneous alignment statistics.  See the [Samtools Documentation](http://www.htslib.org/doc/) for instructions on installing and running Samtools.

### R

The latest version of [R](http://www.r-project.org/) is required for calculating the bacterial p-values, and visualizing the read-hit distribution.  You will also need the following R libraries & packages from [CRAN](http://cran.r-project.org/):

- [ggplot2](http://ggplot2.org/)
- [reshape2](http://cran.r-project.org/web/packages/reshape2/index.html)
- [plyr](http://cran.r-project.org/web/packages/plyr/index.html)
- [untb](http://cran.r-project.org/web/packages/untb/index.html)
- [vegan](http://cran.r-project.org/web/packages/vegan/index.html)
- [multtest (Bioconductor)](http://www.bioconductor.org/packages/release/bioc/html/multtest.html)



## Contact

Contact [Camilo Valdes](mailto:cvaldes3@miami.edu) for pull requests, bug reports, good jokes and coffee recipes. :smiley:

### Maintainers

- [Camilo Valdes](mailto:cvaldes3@miami.edu)


### Collaborators

- [Jennifer Clarke](mailto:jclarke3@unl.edu)
- [Bertrand Clarke](mailto:bclarke3@unl.edu)


## License

The software in this repository is available under the GNU GENERAL PUBLIC LICENSE, Version 3.  See the LICENSE file for more information.
