#!/usr/bin/env python

localrules:
    get_random_motifs,

#bedtools intersect peaks with fimo motifs
#0. with the summit of the peak as reference and strand information from annotation, extend annotation to upstream and 'downstream' distances
#1. merge overlapping (but not book-ended) features
#2. intersect with motif file
rule get_overlapping_motifs:
    input:
        peaks = "diff_binding/peaks/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-{category}-{direction}.narrowpeak",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        motifs = build_annotations("motifs/" + config["genome"]["name"] + "_all_dna_motifs.bed")
    output:
        "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{category}-{direction}-allFIMOresults.tsv.gz",
    params:
        distance = config["motifs"]["search-distance"],
    log:
        "logs/get_upstream_motifs/get_upstream_motifs-{condition}-v-{control}-{norm}-{category}-{direction}-{factor}.log"
    shell: """
        (sed 's/Inf/255/g' {input.peaks} | \
        awk 'BEGIN{{FS=OFS="\t"}}{{$2=$2+$10; $3=$2+1; print $0}}' | \
        bedtools slop -b {params.distance} -i stdin -g <(faidx {input.fasta} -i chromsizes) | \
        sort -k1,1 -k2,2n | \
        bedtools merge -d -1 -c 4,5,6,7,8,9 -o collapse,max,first,max,max,max -i stdin | \
        sort -k1,1 -k2,2n | \
        bedtools intersect -a stdin -b {input.motifs} -sorted -F 1 -wao | \
        cut -f18 --complement | \
        cat <(echo -e "chrom\tregion_start\tregion_end\tpeak_id\tpeak_score\tpeak_strand\tpeak_lfc\tpeak_logpval\tpeak_logqval\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | \
        pigz -f > {output}) &> {log}
        """

rule get_random_motifs:
    input:
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        motifs = build_annotations("motifs/" + config["genome"]["name"] + "_all_dna_motifs.bed"),
    output:
        "motifs/random_sequences-allFIMOresults.tsv.gz"
    params:
        window = 2*(config["motifs"]["search-distance"])+1,
        n = 6000
    log:
        "logs/get_random_motifs.log"
    shell: """
        (bedtools random -l {params.window} -n {params.n} -g <(faidx {input.fasta} -i chromsizes) | \
        sort -k1,1 -k2,2n | \
        bedtools merge -d -1 -c 4,5,6 -o collapse,first,first -i stdin | \
        sort -k1,1 -k2,2n | \
        bedtools intersect -a stdin -b {input.motifs} -sorted -F 1 -wao | \
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, 0, 0, 0, $7, $8, $9, $10, $11, $12, $13, $14}}' | \
        cat <(echo -e "chrom\tregion_start\tregion_end\tpeak_id\tpeak_score\tpeak_strand\tpeak_lfc\tpeak_logpval\tpeak_logqval\tmotif_chrom\tmotif_start\tmotif_end\tmotif_id\tmotif_logpval\tmotif_strand\tmotif_alt_id\tmatch_sequence") - | \
        pigz -f > {output}) &> {log}
        """

rule test_motif_enrichment:
    input:
        fimo_pos = "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{category}-{direction}-allFIMOresults.tsv.gz",
        fimo_neg = lambda wc: "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{category}-unchanged-allFIMOresults.tsv.gz".format(**wc) if wc.negative=="unchanged" else "motifs/random_sequences-allFIMOresults.tsv.gz",
    output:
        tsv = "motifs/{condition}-v-{control}/{norm}/{category}/{negative}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{category}-{direction}-v-{negative}-motif_enrichment.tsv",
        plot = "motifs/{condition}-v-{control}/{norm}/{category}/{negative}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{category}-{direction}-v-{negative}-motif_enrichment.svg",
    params:
        fdr = -log10(config["motifs"]["enrichment-fdr"]),
        direction = lambda wc: "upregulated" if wc.direction=="up" else "downregulated"
    conda:
        nexuspipe("envs/tidyverse.yaml")
    script:
        "../scripts/motif_enrichment.R"

##TODO: this all changes with MEME suite 5.0
##0. extend peak summit annotation to upstream and downstream distances
##1. if multiple annotations overlap on same strand, keep the one that is the most significant (avoid multiple-counting poorly called peaks erroneously split into multiple peaks)
#rule get_meme_sequences:
#    input:
#        peaks = "diff_exp/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}-summits.bed",
#        fasta = config["genome"]["fasta"]
#    output:
#        "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.fa"
#    params:
#        upstr = config["motifs"]["meme-chip"]["upstream"],
#        dnstr = config["motifs"]["meme-chip"]["downstream"]
#    log: "logs/get_meme_sequences/get_meme_sequences_{condition}-v-{control}-{norm}-{category}-{direction}.log"
#    shell: """
#        (bedtools slop -l {params.upstr} -r {params.dnstr} -s -i {input.peaks} -g <(faidx {input.fasta} -i chromsizes) | bedtools cluster -s -d 0 -i stdin | bedtools groupby -g 7 -c 5 -o max -full -i stdin | sort -k4,4V | bedtools getfasta -name+ -s -fi {input.fasta} -bed stdin > {output}) &> {log}
#        """

#rule meme_chip:
#    input:
#        seq = "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.fa"
#        genome_fasta = config["genome"]["fasta"],
#        dbs = "motifs/allmotifs.meme"
#    output:
#        "motifs/{condition}-v-{control}/{norm}/{category}/{background}/meme_chip/meme-chip.html"
#    params:
#        nmeme = lambda wc: int(1e5//(config["motifs"]["meme-chip"]["upstream"] + config["motifs"]["meme-chip"]["downstream"])) if wc.region=="upstream" else int(1e5//config["motifs"]["meme-chip"]["peak-ccut"]),
#        ccut = lambda wc: int(config["motifs"]["meme-chip"]["upstream"] + config["motifs"]["meme-chip"]["downstream"]) if wc.region=="upstream" else config["motifs"]["meme-chip"]["peak-ccut"],
#        meme_mode = config["motifs"]["meme-chip"]["meme-mode"],
#        meme_nmotifs = config["motifs"]["meme-chip"]["meme-nmotifs"],
#    conda: "envs/meme_chip.yaml"
#    shell: """
#        meme-chip -oc motifs/meme/{wildcards.condition}-v-{wildcards.control}/{wildcards.norm}/{wildcards.region}/{wildcards.condition}-v-{wildcards.control}-results-{wildcards.norm}-{wildcards.direction}-{wildcards.category}-{wildcards.region}-meme_chip -bfile <(fasta-get-markov {input.genome_fasta} -m 1) -nmeme {params.nmeme} -norand -ccut {params.ccut} -meme-mod {params.meme_mode} -meme-nmotifs {params.meme_nmotifs} -centrimo-local {params.dbs} {input.seq}
#        """

##NOTE: below are rules for visualizing motif occurrences, but need to have a way to do this efficiently/in an interpretable way for thousands of motifs
## rule get_motif_coverage:
##     input:
##         bed = "motifs/.{motif}.bed", #this is sorted when created
##         chrsizes = config["genome"]["chrsizes"]
##     output:
##         bg = "motifs/coverage/{motif}.bedgraph",
##         bw = "motifs/coverage/{motif}.bw",
##     shell: """
##         cut -f1-6 {input.bed} | bedtools genomecov -bga -i stdin -g {input.chrsizes} | sort -k1,1 -k2,2n > {output.bg}
##         bedGraphToBigWig {output.bg} {input.chrsizes} {output.bw}
##         """

## rule motif_matrix:
##     input:
##         annotation = "diff_exp/{condition}-v-{control}/{category}/{condition}-v-{control}-results-{norm}-{direction}-{category}.bed",
##         bw = "motifs/coverage/{motif}.bw"
##     output:
##         dtfile = temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.mat"),
##         matrix = temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv"),
##         matrix_gz = "motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv.gz",
##     params:
##         refpoint = "TSS",
##         upstream = config["motifs"]["upstream"] + config["motifs"]["binsize"],
##         dnstream = config["motifs"]["freq-downstream"] + config["motifs"]["binsize"],
##         binsize = config["motifs"]["binsize"],
##         sort = "keep",
##         binstat = "sum"
##     threads : config["threads"]
##     log: "logs/deeptools/computeMatrix-{motif}-{condition}-{control}-{norm}-{direction}-{category}.log"
##     shell: """
##         (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --averageTypeBins {params.binstat} -p {threads}) &> {log}
##         pigz -fk {output.matrix}
##         """

## rule melt_motif_matrix:
##     input:
##         matrix = "motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}.tsv.gz",
##     output:
##         temp("motifs/datavis/{motif}_{condition}-v-{control}_{norm}-{direction}-peaks-{category}-melted.tsv.gz"),
##     params:
##         refpoint = "TSS",
##         binsize = config["motifs"]["binsize"],
##         upstream = config["motifs"]["upstream"],
##     script:
##         "scripts/melt_motif_matrix.R"

## rule cat_motif_matrices:
##     input:
##         expand("motifs/datavis/{motif}_{{condition}}-v-{{control}}_{{norm}}-{direction}-peaks-{category}-melted.tsv.gz", category=CATEGORIES, motif=MOTIFS, direction=["up","down","unchanged"]),
##     output:
##         "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.tsv.gz"
##     log: "logs/cat_matrices/cat_matrices-{condition}-{control}-{norm}.log"
##     shell: """
##         (cat {input} > {output}) &> {log}
##         """

#rule plot_motif_freq:
#    input:
#        "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.tsv.gz"
#    output:
#        "motifs/datavis/allmotifs-{condition}-v-{control}-{norm}.svg"
#    script: "scripts/motif_metagenes.R"

