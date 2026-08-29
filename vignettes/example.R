# -- Packages --------------------------------------------
## utility
library(tidyverse)
library(rlang)
library(devtools)

## example data
library(DESeq2)

## plotting
library(ggplot2)
library(ggrepel)

# -- Source functions ------------------------------------
# via source
source('R/utils.R')
source('R/volcano.R')

# or via devtools if interactive dev
# devtools::load_all()

# -- Generate example DEG data ---------------------------
dds = makeExampleDESeqDataSet(
    n      = 1000,
    m      = 8,
    betaSD = 1
)

dds = DESeq(dds)
res = results(dds) %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    drop_na()

# -- Plot volcano of example DEG data --------------------
plot_volcano(
    # dataframe or file of DEGs
    degs           = res,
    
    plot_version   = "v1",
    plot_name      = "volcano",
    plot_title     = "Condition A vs Condition B",

    # column headers that are in the data
    label_col      = "gene",
    log2FC_col     = "log2FoldChange",
    pval_col       = "padj",

    # conditions to be added to the plot filename
    condition_down = "ConditionA",
    condition_up   = "ConditionB",

    # annotation labels at the top of the plot
    anno_down      = paste0("Gene expression\nhigher in\n", "Condition A"),
    anno_up        = paste0("Gene expression\nhigher in\n", "Condition B"),

    # thresholds for intercept lines and point colors
    padj_thresh    = 0.05,
    log2fc_thresh  = 1,

    # colors for the points
    color_down     = "skyblue",
    color_up       = "orange",
    ns_color       = "gray50",

    # plot limits for x and y axes
    xlims          = c(-(max(res$log2FoldChange) * 1.5), max(res$log2FoldChange) * 1.5),
    ylims          = c(0, max(-log10(res$padj)) * 1.5),

    # y position of annotation labels
    anno_ypos      = max(-log10(res$padj)) * 1.5,

    # nudging of point labels (down) to attempt better autoplacement
    nudge_x_dn     = -0.5,
    nudge_y_dn     = 0.5,

    # nudging of point labels (up) to attempt better autoplacement
    nudge_x_up     = 0.5,
    nudge_y_up     = 0.5,

    # PDF and PNG dimensions (in inches)
    plot_height    = 4,
    plot_width     = 4
)
