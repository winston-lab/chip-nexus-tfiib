#!/usr/bin/env python
import itertools

configfile: "config.yaml"

subworkflow nexuspipe:
    workdir: config["nexus_pipeline"]

configfile: nexuspipe("config.yaml")

subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

configfile: build_annotations("config.yaml")

FACTOR = config["factor"]

CATEGORIES = ["genic", "intragenic", "intergenic"]

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}
GROUPS = set(v["group"] for (k,v) in SAMPLES.items())

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"]]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"]]))

comparisons_si = config["comparisons"]["spikenorm"]
if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in config["comparisons"]["spikenorm"]]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in config["comparisons"]["spikenorm"]]))

include: "rules/chip-nexus_classify_peaks.smk"
include: "rules/chip-nexus_gene_ontology.smk"
# include: "rules/chip-nexus_motifs.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules:
    target

rule target:
    input:
        #categorise peaks
        expand("peakcalling/macs/{group}/{group}_experimental-{factor}-chipnexus_peaks-{category}.narrowpeak", group=GROUPS, factor=FACTOR, category=CATEGORIES),
        #expand(expand("peakcalling/macs/{condition}-v-{control}-{{factor}}-chipnexus-peaknumbers.tsv", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"]),
        #categorize DB peaks
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-diffbind-results-{{category}}-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups), factor=FACTOR, direction=["all","up","unchanged","down"], category=CATEGORIES),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-diffbind-results-{{category}}-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR, direction=["all","up","unchanged","down"], category=CATEGORIES) if comparisons_si else [],
        #genic vs. intragenic
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/genic_v_intragenic/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-genic-v-intragenic.tsv", zip, condition=conditiongroups, control=controlgroups), factor=FACTOR),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/genic_v_intragenic/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-genic-v-intragenic.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR) if comparisons_si else [],
        #DB summary
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-diffbind-mosaic.svg", zip, condition=conditiongroups, control=controlgroups), factor=FACTOR),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-diffbind-mosaic.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR),
        #gene ontology
        expand(expand("gene_ontology/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down", "unchanged"], category=["genic", "intragenic"], factor=FACTOR) if config["run_gene_ontology"] else [],
        expand(expand("gene_ontology/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down", "unchanged"], category=["genic", "intragenic"], factor=FACTOR) if config["run_gene_ontology"] and comparisons_si else [],
        ##enrichment of known motifs
        #expand(expand("motifs/{condition}-v-{control}/libsizenorm/{{category}}/{{negative}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES, factor=FACTOR) if config["motifs"]["run_motif_analyses"] else [],
        #expand(expand("motifs/{condition}-v-{control}/spikenorm/{{category}}/{{negative}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES, factor=FACTOR) if config["motifs"]["run_motif_analyses"] and comparisons_si else [],

rule summarise_diffbind_results:
    input:
        total = nexuspipe("diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-all.tsv"),
        genic = "diff_binding/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-genic-all.tsv",
        intragenic = "diff_binding/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-intragenic-all.tsv",
        intergenic = "diff_binding/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-intergenic-all.tsv",
    output:
        summary_table = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-summary.tsv",
        mosaic = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-mosaic.svg",
        maplot = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-maplot.svg",
        volcano = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-volcano.svg",
        volcano_free = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-volcano-freescale.svg",
    params:
        lfc = config["deseq"]["fold-change-threshold"],
        alpha = config["deseq"]["fdr"]
    conda:
        nexuspipe("envs/tidyverse.yaml")
    script:
        "scripts/plot_diffbind_summary.R"

rule genic_v_intragenic:
    input:
        genic = "diff_binding/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-genic-all.tsv",
        intragenic = "diff_binding/{condition}-v-{control}/{norm}/intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-intragenic-all.tsv",
    output:
        tsv = "diff_binding/{condition}-v-{control}/{norm}/genic_v_intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-genic-v-intragenic.tsv",
        lfc_v_lfc = "diff_binding/{condition}-v-{control}/{norm}/genic_v_intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-genic-v-intragenic-lfc-v-lfc.svg",
        lfc_v_expr = "diff_binding/{condition}-v-{control}/{norm}/genic_v_intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-genic-v-intragenic-lfc-v-expr.svg",
        expr_v_expr = "diff_binding/{condition}-v-{control}/{norm}/genic_v_intragenic/{condition}-v-{control}_{factor}-chipnexus-{norm}-genic-v-intragenic-expr-v-expr.svg",
    conda:
        nexuspipe("envs/tidyverse.yaml")
    script:
        "scripts/tfiib_genic_v_intragenic.R"

