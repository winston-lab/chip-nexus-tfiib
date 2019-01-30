library(tidyverse)
library(broom)
library(viridis)

main = function(condition, control,
                genic_path, intra_path,
                tsv_out, lfc_v_lfc_out, lfc_v_expr_out, expr_v_expr_out,
                factor){
    df = read_tsv(genic_path) %>%
        select(chrom, start, end, name, score, strand,
               log2_FC=log2_foldchange, lfc_SE, log10_padj,
               condition_expr, control_expr, peak_summit,
               feature_name=genic_name) %>%
        inner_join(read_tsv(intra_path) %>%
                       select(chrom, start, end, name, score, strand,
                              log2_FC=log2_foldchange, lfc_SE, log10_padj,
                              condition_expr, control_expr, peak_summit,
                              feature_name=orf_name, atg_to_peak_dist),
                   by = "feature_name", suffix = c("_genic", "_intragenic")) %>%
        write_tsv(tsv_out)

    lfc_v_lfc = ggplot() +
        geom_hline(yintercept = 0, size=0.5, color="grey65") +
        geom_vline(xintercept = 0, size=0.5, color="grey65") +
        geom_point(data = df,
                   aes(x=log2_FC_intragenic, y=log2_FC_genic),
                   size=0.5, alpha=0.6) +
        # stat_binhex(data = df,
        #             aes(x=log2_FC_intragenic, y=log2_FC_genic,
        #                 color=..count..),
        #             geom="point", binwidth=c(0.1, 0.1),
        #             size=0.5, alpha=0.6, fill=NA) +
        geom_smooth(data = df,
                    aes(x=log2_FC_intragenic, y=log2_FC_genic),
                    method="lm", size=0.8, color="#114477", alpha=0.4) +
        geom_label(data = df %>%
                       do(tidy(lm(.$log2_FC_genic ~
                                      .$log2_FC_intragenic))) %>%
                       filter(term != "(Intercept)") %>%
                       mutate(label = paste("slope =", signif(estimate, 2),
                                            "\np =",
                                            scales::scientific(p.value))),
                   aes(label=label),
                   x = min(df[["log2_FC_intragenic"]]),
                   y = max(df[["log2_FC_genic"]]),
                   hjust=0, vjust=1, alpha=0.6, size=10/72*25.4,
                   label.size=0) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        xlab(bquote(intragenic ~ .(factor) ~ log[2] ~ textstyle(frac(.(condition), .(control))))) +
        ylab(bquote(atop(genic ~ .(factor), ~ log[2] ~ textstyle(frac(.(condition), .(control)))))) +
        theme_light() +
        theme(text = element_text(size=12, color="black"),
              plot.title = element_text(size=12),
              strip.background = element_blank(),
              strip.text = element_text(color="black"),
              axis.text = element_text(color="black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
    ggsave(lfc_v_lfc_out, plot=lfc_v_lfc, width=16, height=12, units="cm")

    lfc_v_expr = ggplot() +
        geom_hline(yintercept = 0, size=0.5, color="grey65") +
        geom_vline(xintercept = 0, size=0.5, color="grey65") +
        geom_point(data = df,
                   aes(x=control_expr_genic+1, y=log2_FC_intragenic),
                   size=0.5, alpha=0.6) +
        # stat_binhex(data = df,
        #             aes(x=control_expr_genic+1, y=log2_FC_intragenic,
        #                 color=..count..),
        #             geom="point", binwidth=c(0.05, 0.1),
        #             size=0.5, alpha=0.6, fill=NA) +
        geom_smooth(data = df,
                    aes(x=control_expr_genic+1, y=log2_FC_intragenic),
                    method="lm", size=0.8, color="#114477", alpha=0.4) +
        geom_label(data = df %>%
                       do(tidy(lm(.$log2_FC_intragenic ~
                                      log10(.$control_expr_genic+1)))) %>%
                       filter(term != "(Intercept)") %>%
                       mutate(label = paste("slope =", signif(estimate, 2),
                                            "\np =",
                                            scales::scientific(p.value))),
                   aes(label=label),
                   x = log10(min(df[["control_expr_genic"]])+1),
                   y = max(df[["log2_FC_intragenic"]]),
                   hjust=0, vjust=1, alpha=0.6, size=10/72*25.4,
                   label.size=0) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        scale_x_log10(name = paste("genic", factor, "occupancy in", control)) +
        ylab(bquote(atop(intragenic ~ .(factor), ~ log[2] ~ textstyle(frac(.(condition), .(control)))))) +
        theme_light() +
        theme(text = element_text(size=12, color="black"),
              plot.title = element_text(size=12),
              strip.background = element_blank(),
              strip.text = element_text(color="black"),
              axis.text = element_text(color="black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
    ggsave(lfc_v_expr_out, plot=lfc_v_expr, width=16, height=12, units="cm")

    expr_v_expr = ggplot() +
        geom_hline(yintercept = 0, size=0.5, color="grey65") +
        geom_vline(xintercept = 0, size=0.5, color="grey65") +
        geom_point(data = df,
                    aes(x=control_expr_genic+1, y=condition_expr_intragenic+1),
                   size = 0.5, alpha=0.6) +
        # stat_binhex(data = df,
        #             aes(x=control_expr_genic+1, y=condition_expr_intragenic+1,
        #                 color=..count..),
        #             geom="point", binwidth=c(0.01, 0.01),
        #             size=0.5, alpha=0.6, fill=NA) +
        geom_smooth(data = df,
                    aes(x=control_expr_genic+1, y=condition_expr_intragenic+1),
                    method="lm", size=0.8, color="#114477", alpha=0.4) +
        geom_label(data = df %>%
                       do(tidy(lm(log10(.$condition_expr_intragenic+1) ~
                                      log10(.$control_expr_genic+1)))) %>%
                       filter(term != "(Intercept)") %>%
                       mutate(label = paste("slope =", signif(estimate, 2),
                                            "\np =",
                                            scales::scientific(p.value))),
                   aes(label=label),
                   x = log10(min(df[["control_expr_genic"]])+1),
                   y = log10(max(df[["condition_expr_intragenic"]])+1),
                   hjust=0, vjust=1, alpha=0.6, size=10/72*25.4,
                   label.size=0) +
        scale_color_viridis(option="inferno", guide=FALSE) +
        scale_x_log10(name = paste("genic", factor, "occupancy in", control)) +
        scale_y_log10(name = bquote(atop(intragenic ~ .(factor) ~ occupancy, "in" ~ .(condition)))) +
        theme_light() +
        theme(text = element_text(size=12, color="black"),
              plot.title = element_text(size=12),
              strip.background = element_blank(),
              strip.text = element_text(color="black"),
              axis.text = element_text(color="black"),
              axis.title.y = element_text(angle=0, hjust=1, vjust=0.5))
    ggsave(expr_v_expr_out, plot=expr_v_expr, width=20, height=12, units="cm")
}

main(condition=snakemake@wildcards[["condition"]],
     control=snakemake@wildcards[["control"]],
     genic_path=snakemake@input[["genic"]],
     intra_path=snakemake@input[["intragenic"]],
     tsv_out=snakemake@output[["tsv"]],
     lfc_v_lfc_out=snakemake@output[["lfc_v_lfc"]],
     lfc_v_expr_out=snakemake@output[["lfc_v_expr"]],
     expr_v_expr_out=snakemake@output[["expr_v_expr"]],
     factor = snakemake@wildcards[["factor"]])

