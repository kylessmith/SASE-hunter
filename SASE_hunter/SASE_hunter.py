from random import random
from collections import defaultdict
import numpy as np
from numpy import argsort
from scipy.stats import chisqprob
from pybedtools import BedTool, Interval, IntervalFile
from operator import attrgetter
import argparse, os, tempfile, fisher
import sys
#from bx.bbi.bigwig_file import BigWigFile



def fisher_combine(pvals):
    """ combined fisher probability with correction
    Use fdr correction for 25 comparisons using rpy2"""
    if all(p == "NA" for p in pvals): return np.nan
    pvals = [p for p in pvals if p != "NA"]
    if len(pvals) == 1: return pvals[0]
    s = -2 * np.sum(np.log(pvals))
    return chisqprob(s, 2 * len(pvals))

def shuff_score(n_variants_in_seed, n_variants_in_flank,
        n_shuffles, seed_scores, flank_scores):
    sims = []
    for i in xrange(n_shuffles):

        sims.append((
            np.random.choice(seed_scores, n_variants_in_seed,
                    replace=False).sum(),
            np.random.choice(flank_scores, n_variants_in_flank,
                    replace=False).sum()))

    return sims

def sim_score(n_variants_in_seed, n_variants_in_flank,
        n_shuffles, seed_scores, flank_scores, obs_sum_in_region,
        obs_sum_in_flanks):

    r = shuff_score(n_variants_in_seed, n_variants_in_flank,
        n_shuffles, seed_scores, flank_scores)

    assert obs_sum_in_flanks > 0
    obs_ratio = obs_sum_in_region / obs_sum_in_flanks
    pval = sum((sum_in_region / sum_out_region) >= obs_ratio for sum_in_region,
                      sum_out_region in r) / float(n_shuffles)
    return {'score': pval }


def shuff(n_total_variants, n_shuffles, total_size, select_size):
    """
    n_total_variants: the number of variants that were observed
                      in the region
    n_shuffles: number of shuffles to perform
    total_size: the total region size to consider (46kb
    select_size: the region-size of interest
    returns: a list of tuples of lenght n_shuffles where each tuple
             contains (n_in_region, n_out_region) for that shuffle.
             likely, this will be used for the null distribution.
    """

    sims = []
    for i in xrange(n_shuffles):

        # create empty array
        sim_variants = [0] * total_size
        # simulate the variants
        for v in xrange(n_total_variants):
            #pos = randint(0, total_size - 1) # slow
            pos = int(random() * total_size) # faster
            sim_variants[pos] += 1.0

        # count number of times it landed in first [select_size] posns
        sum_in_region = sum(sim_variants[npos] for npos in xrange(select_size))
        # subtract out the ones in the region
        sum_out_region = sum(sim_variants) - sum_in_region

        sims.append((sum_in_region, sum_out_region))

    return sims


def fisher_compare(var_in_region, var_in_flanks, region_size, total_size,
                   n_shuffles, test):
   '''
    impliment fisher to get pvalue and store in a dictionary in ref_count
    calls function shuff which conducts random shuffle
   '''

   pvalue = (fisher.pvalue(var_in_region,
                     var_in_flanks,
                     region_size,
                     total_size - region_size))
   record = {}
   if test in ('fisher', "both"):
       record['fisher'] = pvalue.right_tail

   if test in ('permutation', "both"):
       r = shuff(var_total, n_shuffles=n_shuffles,total_size=total_size,
                 select_size=region_size)
       r_pvalue = sum(sim_in_region >= var_in_region for sim_in_region,
                      sim_out_region in r) / float(n_shuffles)
       record['permutation'] = r_pvalue
   return record


def analyze_intervals(include_file, seed_variants, variant_intervals,
        left_distance, right_distance, test, n_shuffles, out_file, full_format,
        score=None):
    '''
    calculate p-value of obs in seed vs expected from flanks.
    '''

    out = open(out_file, 'w') if isinstance(out_file, basestring) else out_file

    # flanks is a generator
    flankiter = get_flanks(include_file, seed_variants, left_distance,
                        right_distance)
    fmt = "{chrom}\t{start}\t{end}\t{flank_bases}\t{seed_bases}\t{n_samples}\t{info}"
    if full_format:
        fmt += "\t{sample}\t{seed_mutations}\t{flank_mutations}"
    if test in ("fisher", "both"):
        fmt += "\t{fisher}\t{sig}:{total_analyzed}"
    if test in ("permutation", "both"):
        fmt += "\t{permutation}"
    fmt += "\n"

    out.write(fmt.replace("{", "").replace("}", ""))

    NA = dict(fisher="NA", permutation="NA", sample="NA", score="NA")
    for flank_hits, seed_hits, total_bases, orig_seed in flankiter:
        sig = 0
        seen = set()
        variants_in_seed_by_sample = defaultdict(list)
        seed_length = sum(seed.length for seed in seed_hits)
        for s in seed_hits:
            s.file_type = "bed"
            variants_in_seed = [x for x in
                    variant_intervals.all_hits(s) if not x in seen]
            seen.update(variants_in_seed)
            # for each sample, track the variants in flanks
            for variant in variants_in_seed:
                variants_in_seed_by_sample[variant[3]].append(variant)

        variants_in_flanks_by_sample = defaultdict(list)

        # since the same (multi-nucleotide) variant could overlap multiple
        # we use set() so we dont double-count
        seen = set()
        for interval in flank_hits:
            interval.file_type = "bed"
            variants_in_flanks = [x for x in
                    variant_intervals.all_hits(interval) if not x in seen]
            seen.update(variants_in_flanks)
            # for each sample, track the variants in flanks
            for variant in variants_in_flanks:
                variants_in_flanks_by_sample[variant[3]].append(variant)

        samples = set((variants_in_flanks_by_sample.keys() \
                       + variants_in_seed_by_sample.keys()))
        n_samples = len(variants_in_seed_by_sample.keys())
        if int(total_bases) > seed_length and len(samples) != 0:
            pvals = []
            # for all samples in the seed, see how many variants are in flank
            # from that sample, throw to analyze function and get pvalue
            if score is not None:
                seed_scores = np.hstack([
                    score.get_as_array(s.chrom, s.start, s.end) for s in
                    seed_hits])
                flank_scores = np.hstack([
                    score.get_as_array(f.chrom, f.start, f.end) for f in
                        flank_hits ])

            for sample in samples:
                total_variants = len(variants_in_seed_by_sample[sample]) + \
                                 len(variants_in_flanks_by_sample[sample])
                if score is not None:
                    # will need the actual variants, not just the counts.

                    vs = variants_in_seed_by_sample[sample]
                    # get returns [(start, end, value)]
                    obs_sum_in_region = sum(score.get(v.chrom,
                                                      v.start,
                                                      v.end)[0][2] for v in vs)

                    vsf = variants_in_flanks_by_sample[sample]
                    if len(vs) + len(vsf) < 3 or len(vsf) < 1:
                        pvals.append(NA.copy())
                        continue
                    obs_sum_in_flanks = sum(score.get(v.chrom,
                                                      v.start,
                                                      v.end)[0][2] for v in vsf)

                    pvals.append(sim_score(
                        len(variants_in_seed_by_sample[sample]),
                        len(variants_in_flanks_by_sample[sample]),
                        n_shuffles,
                        seed_scores, flank_scores,
                        obs_sum_in_region, obs_sum_in_flanks))
                elif total_variants < 2:
                    pvals.append(NA.copy())
                else:
                    p = fisher_compare(
                          len(variants_in_seed_by_sample[sample]),
                          len(variants_in_flanks_by_sample[sample]),
                          seed_length, total_bases, n_shuffles, test)
                    if p['fisher'] <= 0.05: sig += 1
                    pvals.append(p)

                pvals[-1]['sample'] = sample
            assert pvals, (pvals, samples)
        else:
            pvals = [NA.copy()]
            
        idict = dict(chrom=orig_seed.chrom, start=orig_seed.start,
                     end=orig_seed.end, n_samples=n_samples,
                     info=orig_seed[3], flank_bases=total_bases - seed_length,
                     seed_bases=seed_length, sig=sig, total_analyzed=len(samples))
        if full_format:
            for pdict in pvals:
                idict.update(dict(
                  seed_mutations=len(variants_in_seed_by_sample[pdict['sample']]),
                  flank_mutations=len(variants_in_flanks_by_sample[pdict['sample']])
                ))
                idict.update(pdict)
                out.write(fmt.format(**idict))
        else:
            if test in ("fisher", "both"):
                idict['fisher'] = fisher_combine([p['fisher'] for p in pvals])
            if test in ("permutations", "both"):
                idict['permutation'] = fisher_combine([p['permutation'] for p
                    in pvals])
            if score:
                idict['score'] = fisher_combine([p['score'] for p in pvals])
            out.write(fmt.format(**idict))

    return out.name


def get_flanks(include_file, seeds, up_distance, down_distance):

    # .intervals returns the tree
    include = include_file.intervals

    def bad_cov(intervals):
        coverage = set()
        for i in intervals:
            coverage.update(range(i.start, i.end))
        return len(coverage)

    fmt = "{s.chrom}\t{s.start}\t{s.end}\t{t}\t{left}\t{right}\t{n_intervals}\n"
    for seed in seeds:
        seed_hits = include.all_hits(seed)

        # and either continue or yield None?
        if seed_hits == []: continue
        seed_hits = sorted(seed_hits, key=attrgetter('start'))

        if seed.strand == "-":
            region_left = Interval(seed.chrom, max(0,seed.start - down_distance), seed.start)
            region_right = Interval(seed.chrom, seed.end, seed.end + up_distance)
        else:
            region_left = Interval(seed.chrom, max(0,seed.start - up_distance), seed.start)
            region_right = Interval(seed.chrom, seed.end, seed.end + down_distance)
        # this is needed, unfortunately
        region_right.file_type = region_left.file_type = "bed"

        # this won't handle overlapping include intervals...
        include_left = sorted(include.all_hits(region_left), key=attrgetter('start'))
        include_right = sorted(include.all_hits(region_right), key=attrgetter('start'))

        #assert bad_cov(include_right) == sum(i.length for i in include_right)
        #assert bad_cov(include_left) == sum(i.length for i in include_left)

        if include_left:
            # truncate to not include the seed point so we get unique for left, right
            # and add in the seed at the end
            include_left[-1].end = min(include_left[-1].end, seed.start)
            # adjust end-point so we get exactly pad
            include_left[0].start = max(include_left[0].start, region_left.start)

        if include_right:
            # truncate to region
            include_right[0].start = max(include_right[0].start, seed.end)
            include_right[-1].end = min(include_right[-1].end, region_right.end)

        #assert bad_cov(include_right) == sum(i.length for i in include_right)
        #assert bad_cov(include_left) == sum(i.length for i in include_left)

        assert include_left == [] or include_left[-1].end <= seed.start
        assert include_right == [] or seed.end <= include_right[0].start

        #null_bases = total_bases - l

        #assert bad_cov(include_left) == sum(i.length for i in include_left)
        #assert bad_cov(include_right) == sum(i.length for i in include_right)
        #assert bad_cov(seed_hits) == sum(i.length for i in seed_hits)

        # truncate seed_hits to actual seed region
        seed_hits[0].start = max(seed_hits[0].start, seed.start)
        seed_hits[-1].end = min(seed_hits[-1].end, seed.end)

        #assert bad_cov(seed_hits) == sum(i.length for i in seed_hits)

        flanks = include_left + include_right
        total_bases = sum(i.length for i in flanks) + sum(i.length for i in
                seed_hits)
        #coverage = bad_cov(flanks)
        #assert coverage == total_bases
        yield flanks, seed_hits, total_bases, seed

def bed_generator(in_vcfs, prefix="chr", filter=True, qual_min=1):

    for in_vcf in in_vcfs:
        for line in open(in_vcf, mode = 'r'):
            if line[0] == "#": continue
            line = line.strip().split('\t')

            if line[0] == '#CHROM':
                samples = line[9:]
            if line[6] not in (".", "PASS"): continue
            if line[3] == line[4]: continue

            if float(line[5]) < qual_min: continue
            for sample in samples:
                yield "{prefix}{chrom}\t{start}\t{end}\t{sample}".format(**locals())

def vcf_to_long_bed(*in_vcfs):

    intervals = bed_generator(in_vcfs)
    region_file = BedTool(intervals)

    return region_file.sort()
    
def create_genome_bed(genome, start = "0"):
    
    for line in open(genome,mode = 'r'):
        line = line.strip().split('\t')
        chrom = line[0]
        end = line[1]
        yield "{chrom}\t{start}\t{end}\n".format(**locals())

def main():
    '''
    main function
    takes arguments from the cammandline
    '''
    parser=argparse.ArgumentParser()

    parser.add_argument('--upstream', type=int, help='distance upstream of seed'
            ' to look for flanking regions to compare with `seed`, default = 20000', default=20000)
    parser.add_argument('--downstream', type=int, help='distance downstream to '
            'look for flanking regions, default = 20000', default=20000)
    parser.add_argument('--seed', help='regions of interest, e.g. promoters',
            required=True, metavar="BED")

    parser.add_argument('--exclude', help='regions to be excluded when looking for flanks')
    parser.add_argument('--include', help='regions to be included when looking for flanks')
    parser.add_argument('--test', default = 'fisher', choices = ['fisher','permutation','both'])
    parser.add_argument('--shuffles',type = int, help='number of shuffles to do for permutation analysis',  default = 1000)
    parser.add_argument('--genome', help='the name of the genome file for BEDTools')
    parser.add_argument('--full', help='output full, dataset with per-sample p-values',
            default=False, action='store_true')
    parser.add_argument('--score', metavar="BIGWIG/INT", help='')
    parser.add_argument('variants', help='regions to assign significance e.g.'
            'a list of variants', metavar="BED/VCF", nargs='+')

    args=parser.parse_args()

    seed_region = BedTool(args.seed)

    if args.variants[0][-4:] == '.vcf':
        region_file = vcf_to_long_bed(*args.variants)
    else:
        assert len(args.variants) == 1
        region_file = IntervalFile(args.variants[0])

    out_file = sys.stdout
    genome = args.genome

    if args.exclude:
        include_file = BedTool(args.exclude).complement(g=genome)
    elif args.include:
        include_file = BedTool(args.include)
    else:
        include_file = BedTool(create_genome_bed(genome))

    bw = None
    if args.score:
        bw = BigWigFile(open(args.score))

    analyze_intervals(include_file, seed_region, region_file, args.upstream,
            args.downstream, args.test, args.shuffles, out_file, args.full,
            score=bw)

if __name__=='__main__':
    main()
