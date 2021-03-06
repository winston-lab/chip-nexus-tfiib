#!/usr/bin/env python

localrules:
    classify_genic_peaks,
    classify_intragenic_peaks,
    classify_intergenic_peaks,
    classify_genic_diffbind_peaks,
    classify_intragenic_diffbind_peaks,
    classify_intergenic_diffbind_peaks

peak_fields = "peak_chrom\tpeak_start\tpeak_end\tpeak_name\tpeak_score\tpeak_strand\tpeak_enrichment\tpeak_logpval\tpeak_logqval\tpeak_summit\t"

rule classify_genic_peaks:
    input:
        annotation = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        peaks = nexuspipe("peakcalling/{group}/{group}_experimental-{factor}-chipnexus-idrpeaks.narrowPeak"),
    output:
        table = "peakcalling/{group}/genic/{group}_experimental-{factor}-chipnexus_peaks-genic.tsv",
        narrowpeak = "peakcalling/{group}/genic/{group}_experimental-{factor}-chipnexus_peaks-genic.narrowpeak",
        bed = "peakcalling/{group}/genic/{group}_experimental-{factor}-chipnexus_peaks-genic-summits.bed",
    log:
        "logs/classify_peaks/classify_genic_peaks-{group}-{factor}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo | \
         cut --complement -f17 | \
         cat <(echo -e "{peak_fields}\tgenic_chrom\tgenic_start\tgenic_end\tgenic_name\tgenic_score\tgenic_strand") - | \
         tee {output.table} | \
         scripts/peak_tsv_to_unique_narrowpeak.sh | \
         tee {output.narrowpeak} | \
         scripts/peak_narrowpeak_to_summit_bed.sh > {output.bed}) &> {log}
        """

rule classify_intragenic_peaks:
    input:
        genic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        orf_anno = os.path.abspath(build_annotations(config["genome"]["orf_annotation"])),
        peaks = nexuspipe("peakcalling/{group}/{group}_experimental-{factor}-chipnexus-idrpeaks.narrowPeak"),
    output:
        table = "peakcalling/{group}/intragenic/{group}_experimental-{factor}-chipnexus_peaks-intragenic.tsv",
        narrowpeak = "peakcalling/{group}/intragenic/{group}_experimental-{factor}-chipnexus_peaks-intragenic.narrowpeak",
        bed = "peakcalling/{group}/intragenic/{group}_experimental-{factor}-chipnexus_peaks-intragenic-summits.bed",
    log:
        "logs/classify_peaks/classify_intragenic_peaks-{group}-{factor}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v | \
         bedtools intersect -a stdin -b <(cut -f1-6 {input.orf_anno}) -f 1 -wo | \
         awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$10}} $16=="+"{{$17=summit-$12}} $16=="-"{{$17=$13-(summit+1)}} {{print $0}}' | \
         cat <(echo -e "{peak_fields}\torf_chrom\torf_start\torf_end\torf_name\torf_score\torf_strand\tatg_to_peak_dist") - | \
         tee {output.table} | \
         scripts/peak_tsv_to_unique_narrowpeak.sh | \
         tee {output.narrowpeak} | \
         scripts/peak_narrowpeak_to_summit_bed.sh > {output.bed}) &> {log}
        """

rule classify_intergenic_peaks:
    input:
        intergenic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_intergenic-regions.bed"),
        transcript_anno = os.path.abspath(build_annotations(config["genome"]["transcript_annotation"])),
        orf_anno = os.path.abspath(build_annotations(config["genome"]["orf_annotation"])),
        genic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        peaks = nexuspipe("peakcalling/{group}/{group}_experimental-{factor}-chipnexus-idrpeaks.narrowPeak"),
    output:
        table = "peakcalling/{group}/intergenic/{group}_experimental-{factor}-chipnexus_peaks-intergenic.tsv",
        narrowpeak = "peakcalling/{group}/intergenic/{group}_experimental-{factor}-chipnexus_peaks-intergenic.narrowpeak",
        bed = "peakcalling/{group}/intergenic/{group}_experimental-{factor}-chipnexus_peaks-intergenic-summits.bed",
    log:
        "logs/classify_peaks/classify_intergenic_peaks-{group}-{factor}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | \
         bedtools intersect -a stdin -b {input.intergenic_anno} -u | \
         cat <(echo -e "{peak_fields}") - | \
         tee {output.table} | \
         scripts/peak_tsv_to_unique_narrowpeak.sh | \
         tee {output.narrowpeak} | \
         scripts/peak_narrowpeak_to_summit_bed.sh > {output.bed}) &> {log}
        """

# rule peakstats:
#     input:
#         expand("peakcalling/macs/{category}/{group}-exp-peaks-{category}.tsv", group=GROUPS, category=CATEGORIES),
#     output:
#         table = "peakcalling/macs/{condition}-v-{control}-{factor}-chipnexus-peaknumbers.tsv",
#         size = "peakcalling/macs/{condition}-v-{control}-{factor}-chipnexus-peaksizes.svg",
#         dist = "peakcalling/macs/{condition}-v-{control}-{factor}-chipnexus-peakdistances.svg"
#     params:
#         groups = lambda wc: [g for sublist in zip(controlgroups, conditiongroups) for g in sublist] if wc.condition=="all" else [wc.control, wc.condition],
#         prefix = config["combinedgenome"]["experimental_prefix"]
#     conda: "../envs/tidyverse.yaml"
#     script:
#         "scripts/peakstats.R"

rule classify_genic_diffbind_peaks:
    input:
        annotation = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        narrowpeak = nexuspipe("diff_binding/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-{direction}.narrowpeak"),
        results = nexuspipe("diff_binding/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-{direction}.tsv"),
    output:
        results = "diff_binding/peaks/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-genic-{direction}.tsv",
        narrowpeak = "diff_binding/peaks/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-genic-{direction}.narrowpeak",
        bed = "diff_binding/peaks/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-genic-{direction}-summits.bed",
    log :
        "logs/classify_diffbind_peaks/classify_genic_diffbind_peaks-{condition}-v-{control}_{norm}-{direction}-{factor}.log"
    shell: """
        (tail -n +2 {input.results} | \
         paste - <(cut -f10 {input.narrowpeak}) | \
         bedtools intersect -a stdin -b {input.annotation} -wo | \
         cut --complement -f22 | \
         cat <(paste <(head -n 1 {input.results}) <(echo -e "peak_summit\tgenic_chrom\tgenic_start\tgenic_end\tgenic_name\tgenic_score\tgenic_strand")) - | \
         tee {output.results} | \
         scripts/diffexp_tsv_to_unique_narrowpeak.sh | \
         tee {output.narrowpeak} | \
         scripts/peak_narrowpeak_to_summit_bed.sh > {output.bed}) &> {log}
        """

rule classify_intragenic_diffbind_peaks:
    input:
        genic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        orf_anno = os.path.abspath(build_annotations(config["genome"]["orf_annotation"])),
        narrowpeak = nexuspipe("diff_binding/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-{direction}.narrowpeak"),
        results = nexuspipe("diff_binding/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-{direction}.tsv"),
    output:
        results = "diff_binding/peaks/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-intragenic-{direction}.tsv",
        narrowpeak = "diff_binding/peaks/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-intragenic-{direction}.narrowpeak",
        bed = "diff_binding/peaks/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-intragenic-{direction}-summits.bed",
    log:
        "logs/classify_diffbind_peaks/classify_intragenic_diffbind_peaks-{condition}-v-{control}_{norm}-{direction}-{factor}.log"
    shell: """
        (tail -n +2 {input.results} | \
         paste - <(cut -f10 {input.narrowpeak}) | \
         bedtools intersect -a stdin -b {input.genic_anno} -v | \
         bedtools intersect -a stdin -b <(cut -f1-6 {input.orf_anno}) -f 1 -wo | \
         awk 'BEGIN{{FS=OFS="\t"}} {{summit=$2+$15}} $21=="+"{{$22=summit-$17}} $21=="-"{{$22=$18-(summit+1)}} {{print $0}}' | \
         cat <(paste <(head -n 1 {input.results}) <(echo -e "peak_summit\torf_chrom\torf_start\torf_end\torf_name\torf_score\torf_strand\tatg_to_peak_dist")) - | \
         tee {output.results} | \
         scripts/diffexp_tsv_to_unique_narrowpeak.sh | \
         tee {output.narrowpeak} | \
         scripts/peak_narrowpeak_to_summit_bed.sh > {output.bed}) &> {log}
        """

rule classify_intergenic_diffbind_peaks:
    input:
        intergenic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_intergenic-regions.bed"),
        transcript_anno = os.path.abspath(build_annotations(config["genome"]["transcript_annotation"])),
        orf_anno = os.path.abspath(build_annotations(config["genome"]["orf_annotation"])),
        genic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        narrowpeak = nexuspipe("diff_binding/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-{direction}.narrowpeak"),
        results = nexuspipe("diff_binding/peaks/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-{direction}.tsv"),
    output:
        results = "diff_binding/peaks/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-intergenic-{direction}.tsv",
        narrowpeak = "diff_binding/peaks/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-intergenic-{direction}.narrowpeak",
        bed = "diff_binding/peaks/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-peaks-diffbind-results-intergenic-{direction}-summits.bed",
    log:
        "logs/classify_diffbind_peaks/classify_intragenic_diffbind_peaks-{condition}-v-{control}_{norm}-{direction}-{factor}.log"
    shell: """
        (tail -n +2 {input.results} | \
         paste - <(cut -f10 {input.narrowpeak}) | \
         bedtools intersect -a stdin -b {input.transcript_anno} {input.orf_anno} {input.genic_anno} -v | \
         bedtools intersect -a stdin -b {input.intergenic_anno} -u | \
         cat <(paste <(head -n 1 {input.results}) <(echo "peak_summit")) - | \
         tee {output.results} | \
         scripts/diffexp_tsv_to_unique_narrowpeak.sh | \
         tee {output.narrowpeak} | \
         scripts/peak_narrowpeak_to_summit_bed.sh > {output.bed}) &> {log}
        """

