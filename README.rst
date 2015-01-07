Software to identify regions of interest with a higher than expected number of mutations than
the near-by regions. 

Data Format
===========

Input files must be in BED format(first columns are chrome, start, stop)::

    chr1    11873   14409

Simple Example:

    https://raw.githubusercontent.com/kylessmith/SASE-hunter/master/example/tier3.bed

More Complicated Example:

    https://github.com/kylessmith/SASE-hunter/blob/master/example/promoters.bed

NOTE: Most input files are assumed to be *sorted* BED files

Invocation
==========

Running the following command will result in a more detailed help message::

    $ python -m SASE_hunter -h

Gives::

      --upstream UPSTREAM   distance upstream of seed to look for flanking regions
                            to compare with `seed`
      --downstream DOWN     distance downstream to look for flanking regions
      --seed BED            regions of interest, e.g. promoters
      --exclude EXCLUDE     regions to be excluded when looking for flanks
      --include INCLUDE     regions to be included when looking for flanks
      --test {fisher,permutation,both}
      --shuffles SHUFFLES   number of shuffles to do for permutation analysis
      --genome GENOME       the name of the genome file for BEDTools
      --full                output full, dataset with per-sample p-values
      variants BED/VCF      variants to shuffle. can be multiple VCF files.

QuickStart
==========

If your files are in sorted BED format, you want to analyze 20000 base pairs upstream and downstream of seed regions,
 give an include file for the flanking analysis, and analyze with the fisher’s exact test.


mutations in promoters for melanoma
-----------------------------------
::

    $ python -m SASE_hunter \
        --upstream 20000 \
        --downstream 20000 \
        --seed example/promoters.bed \
        --include example/tier3.bed \
        --test fisher \
        --genome example/hg19.genome \
        example/mutations.bed \
        > output.txt

The output will be shown in the following columns::

    chrom  start   end flank_bases n_samples   info    fisher

Where the last column is the p-value from the Fisher's exact test with
a contingency table created from the number of variants in the seed region
compared to the number of variants in the flanking regions (derived from
nearby regions in the `-include` argument relative to their size in bases.

The above command will find the accelerated regions (Promoters) that are
mutated more often than the surrounding 20000 base pairs upstream and
downstream.

Installation
============

pip can be used to install by::

    pip install SASE_hunter

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  

    https://continuum.io/downloads

Dependencies can be installed by::

    pip install -r requirements.txt

SASE_hunter also depends on BEDTools which is available from https://github.com/arq5x/bedtools2/
