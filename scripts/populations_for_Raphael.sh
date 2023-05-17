#!/bin/bash
#SBATCH -A snic2017-7-413
#SBATCH -p core -n 20
#SBATCH -t 240:00:00
#SBATCH -J Raphael_populations_vs_std
#SBATCH -o std_out_%j.txt
#SBATCH -e std_err_%j.txt
#SBATCH -C mem1TB


# populations 2.41
# Usage:
# populations -P dir [-O dir] [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)
# populations -V vcf -O dir [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)

#   -P,--in-path: path to a directory containing Stacks ouput files.
#   -V,--in-vcf: path to a standalone input VCF file.
#   -O,--out-path: path to a directory where to write the output files. (Required by -V; otherwise defaults to value of -P.)
#   -M,--popmap: path to a population map. (Format is 'SAMPLE1 \t POP1 \n SAMPLE2 ...'.)
#   -t,--threads: number of threads to run in parallel sections of code.
#   --batch-size [int]: the number of loci to process in a batch (default: 10,000 in de novo mode; in reference mode, one chromosome
#                       per batch). Increase to speed analysis, uses more memory, decrease to save memory).


# These steps, before loading the bioinfo-tools module, are because
# the existing .tags. .snps. .alleles. and .matches. files for each sample
# were created more or less individually, so their sample IDs would collide
# later when running them through the pipeline towards populations.
# As of 09-06-2020 these steps seem to have been successfully performed and are
# not needed anymore. Keeping them for science. WRONGGGGGGGG

#######################################
# trap "chmod 555 $out" EXIT
# out="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/inputs"
# # out="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/inputs/Corrected_IDs"
# chmod 755 "$out" ; rm -f "$out"/*
# cd "$out"

# counter=0 # Should this start on 0?
# while read -d '' -r sample ; do

# base=$(basename "$sample")

# if [[ "$base" == 'BO_rep4' ]] ; then continue ; fi # We don't include this one

# # # alleles
# sed "s/^[0-9][0-9]*/$counter/" <(gzip -cd "$sample".alleles.tsv.gz) > "$out"/"$base".alleles.tsv
# echo "Done with alleles for $base"

# # # snps
# sed "s/^[0-9][0-9]*/$counter/" <(gzip -cd "$sample".snps.tsv.gz) > "$out"/"$base".snps.tsv
# echo "Done with snps for $base"

# # # tags
# sed "s/^[0-9][0-9]*/$counter/" <(gzip -cd "$sample".tags.tsv.gz) > "$out"/"$base".tags.tsv
# echo "Done with tags for $base"

# # # matches
# path2file="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/sstacks/outputs/std_both"/"$base".matches.tsv.gz
# sed "s/^\([0-9][0-9]*\)	[0-9][0-9]*/\1	$counter/" <(gzip -cd "$path2file") > "$out"/"$base".matches.tsv
# echo "Done with matches for $base"

# counter=$((counter+1))
# done < /crex/proj/uppstore2017173/b2014326/Raphael/Analysis/.abspaths_to_samples_without_extensions.nuldelimited
# gzip *tsv
#######################################



module load bioinfo-tools
module load Stacks/2.41

reads="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/reads"

projdir="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/inputs"
# projdir="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/inputs/Corrected_IDs"

# outdir="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/outputs"

popmap="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/bestpopmap0.tsv"
# popmap="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/bestpopmap1.tsv"
# popmap="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/alt_popmap.tsv"
# popmap="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/popmap.tsv"

outdirtec="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/outputs_tec"
popmaptec="/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/populations/tec_popmap.tsv"


# trap "chmod 555 $projdir $projdir $outdir" EXIT
# trap "chmod 555 $projdir $outdir $outdirtec" EXIT
trap "chmod 555 $projdir" EXIT

# populations -P "$projdir" -O "$outdir" -M "$popmap" (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)
# populations -P "$projdir" -O "$outdir" -M "$popmap" -p 73 --vcf -f p_value -t 10 -p 3 -r "$r" --fstats --structure --write_single_snp
# populations -P "$projdir" -O "$outdir" -M "$popmap" -p 73 --vcf

chmod 755 "$projdir" 
# chmod 755 "$projdir" "$projdir" "$outdir"
# chmod 755 "$projdir" "$outdir" "$outdirtec"

# tsv2bam -P "$projdir" -M "$popmap" -t 10 -R "$reads"
# gstacks -P "$projdir" -M "$popmap" -t 10 --write-alignments
/crex/proj/uppstore2017173/b2014326/Raphael/Analysis/STACKS_2.53/compiled/bin/gstacks -P "$projdir" -M "$popmap" -t 1 --write-alignments --details
# gstacks -P "$projdir" -M "$popmap" -t 1 --write-alignments
# populations -P "$projdir" -O "$outdir" -M "$popmap" -t 10 --vcf --plink --hzar --genepop
# populations -P "$projdir" -O "$outdirtec" -M "$popmaptec" -t 10 --vcf --plink --hzar --genepop


# populations -P "$projdir" -O "$outdir" -M "$popmap" -p 80 -t 10 --vcf