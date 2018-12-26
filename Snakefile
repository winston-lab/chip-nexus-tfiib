#!/usr/bin/env python

configfile: "config.yaml"

subworkflow nexuspipe:
    workdir: config["nexus_pipeline"]

CATEGORIES = ["genic", "intragenic", "intergenic"]

include: "rules/chip-nexus_classify_peaks.smk"
include: "rules/chip-nexus_gene_ontology.smk"
include: "rules/chip-nexus_motifs.smk"

onsuccess:
    nexuspipe(shell("(./mogrify.sh) > mogrify.log"))

localrules:
    target

rule target:
        #categorise peaks
        expand("peakcalling/macs/{group}/{group}_experimental-{factor}-chipnexus_peaks-{category}.narrowpeak", group=GROUPS, factor=FACTOR, category=CATEGORIES),
        #expand(expand("peakcalling/macs/{condition}-v-{control}-{{factor}}-chipnexus-peaknumbers.tsv", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"]),
        ##categorize DB peaks
        #expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-diffbind-results-{{category}}-{{direction}}.narrowpeak", condition=conditiongroups, control=controlgroups), factor=FACTOR, direction=["all","up","unchanged","down"], category=CATEGORIES),
        #expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-diffbind-results-{{category}}-{{direction}}.narrowpeak", condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR, direction=["all","up","unchanged","down"], category=CATEGORIES) if comparisons_si else [],
        ###DB summary
        ##expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-libsizenorm-diffbind-summary.svg", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"]),
        ##expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-spikenorm-diffbind-summary.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"]),
        ##gene ontology
        #expand(expand("gene_ontology/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups, control=controlgroups), direction=["up", "down", "unchanged"], category=["genic", "intragenic"], factor=FACTOR) if config["run_gene_ontology"] else [],
        #expand(expand("gene_ontology/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-{{category}}-{{direction}}-gene-ontology-enriched-all.svg", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up", "down", "unchanged"], category=["genic", "intragenic"], factor=FACTOR) if config["run_gene_ontology"] and comparisons_si else [],
        ##enrichment of known motifs
        #expand(expand("motifs/{condition}-v-{control}/libsizenorm/{{category}}/{{negative}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups, control=controlgroups), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES, factor=FACTOR) if config["motifs"]["run_motif_analyses"] else [],
        #expand(expand("motifs/{condition}-v-{control}/spikenorm/{{category}}/{{negative}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-{{category}}-{{direction}}-v-{{negative}}-motif_enrichment.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), direction=["up","down"], negative=["unchanged", "random"], category=CATEGORIES, factor=FACTOR) if config["motifs"]["run_motif_analyses"] and comparisons_si else [],

# rule summarise_db_results:
#     input:
#         total = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-all.tsv",
#         genic = "diff_binding/{condition}-v-{control}/genic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-genic.tsv",
#         intragenic = "diff_binding/{condition}-v-{control}/intragenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-intragenic.tsv",
#         intergenic = "diff_binding/{condition}-v-{control}/intergenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-intergenic.tsv",
#     output:
#         summary = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-summary.svg",
#         maplot = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-maplot.svg",
#         volcano = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-volcano.svg",
#         volcano_free = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-volcano-freescale.svg",
#     params:
#         lfc = config["deseq"]["fold-change-threshold"],
#         alpha = config["deseq"]["fdr"]
#     conda: "../envs/tidyverse.yaml"
#     script: "../scripts/plot_diffbind_summary.R"

