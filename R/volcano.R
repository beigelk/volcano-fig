#' Generate a volcano plot from differential expression results
#'
#' Reads a tab/whitespace-delimited differential expression results file,
#' classifies genes as up- or down-regulated based on log2 fold change and
#' the \code{padj_thresh} significance threshold, logs the top 10 genes in
#' each direction via \code{\link{log_msg}}, and produces an annotated,
#' labeled volcano plot. The finished plot is displayed and saved to disk as
#' both PDF and PNG via \code{\link{plot_pdf_and_png}}.
#'
#' Genes with a p-value below \code{padj_thresh} and a positive log2 fold
#' change are colored/labeled with \code{name_up}; genes with a p-value
#' below \code{padj_thresh} and a negative log2 fold change are
#' colored/labeled with \code{name_down}; all other genes are marked as
#' non-significant (\code{"n.s."}). The points that are outlined and labeled
#' on the plot are restricted further to genes passing both
#' \code{padj_thresh} and \code{log2fc_thresh}, taking the top 10 most
#' extreme genes in each direction (by log2 fold change) via
#' \code{ggrepel::geom_label_repel}. Dashed threshold lines are drawn at
#' \code{padj_thresh} (horizontal) and at \eqn{\pm}\code{log2fc_thresh}
#' (vertical) for visual reference.
#'
#' @param degs Character. Path to the differential expression
#'   results file. Must be readable by \code{read.table} with
#'   \code{header = TRUE} and must contain a column \code{X} (gene
#'   identifiers), which is renamed to \code{gene}.
#' @param plot_version Character/numeric. Version label appended to the
#'   output filename, e.g. \code{"v1"}.
#' @param plot_dirname Character. Directory name used (together with
#'   \code{plot_version}, \code{name_down}, \code{name_up},
#'   and \code{plot_filebase}) to construct the output PDF/PNG filenames, in the
#'   form \code{<plot_dirname>/<plot_filebase>__<name_up>_v_<name_down>__volcano_<plot_version>.pdf} / \code{.png}.
#' @param plot_filebase Character. String used as the base/root for the 
#'   filename of the plot, e.g. a project name or other informative prefix for the filename.
#' @param plot_title Character. Title displayed on the plot and used in the
#'   log messages announcing the top gene tables.
#' @param label_col Character. Name of the column in \code{degs}
#'   used to label the top genes (e.g. gene symbol).
#' @param log2FC_col Character. Name of the column in \code{degs}
#'   containing log2 fold change values.
#' @param pval_col Character. Name of the column in \code{degs}
#'   containing the p-value or adjusted p-value used for significance
#'   filtering (against \code{padj_thresh}) and for the y-axis
#'   (\code{-log10(pval_col)}).
#' @param n_labels Character. Number of items to label from \code{degs}
#'   in each direction, used to label the top genes (e.g. gene symbol).
#' @param name_down Character. Value used to label/color genes with a
#'   significant, negative log2 fold change (i.e. the \code{color} grouping
#'   value, matched in \code{scale_color_manual}). Compare with
#'   \code{anno_down}, which is the human-readable text shown in the plot
#'   annotation.
#' @param name_up Character. Value used to label/color genes with a
#'   significant, positive log2 fold change (i.e. the \code{color} grouping
#'   value, matched in \code{scale_color_manual}). Compare with
#'   \code{anno_up}, which is the human-readable text shown in the plot
#'   annotation.
#' @param anno_down Character. Label describing the group/condition in which
#'   expression is higher for down-regulated (negative log2FC) genes; used in
#'   the down-regulated annotation text (colored with \code{color_down}) and
#'   in the log messages.
#' @param anno_up Character. Label describing the group/condition in which
#'   expression is higher for up-regulated (positive log2FC) genes; used in
#'   the up-regulated annotation text (colored with \code{color_up}) and in
#'   the log messages.
#' @param color_down Character. Color (name or hex code) used for the
#'   down-regulated group: point fill for the outlined top genes, the
#'   \code{name_down} entry in the color scale, and the down-regulated
#'   annotation text.
#' @param color_up Character. Color (name or hex code) used for the
#'   up-regulated group: point fill for the outlined top genes, the
#'   \code{name_up} entry in the color scale, and the up-regulated
#'   annotation text.
#' @param ns_color Character. Color (name or hex code) used for non-
#'   significant (\code{"n.s."}) points in the color scale.
#' @param padj_thresh Numeric. Adjusted p-value (or p-value) significance
#'   threshold used consistently throughout: to classify/color genes as
#'   up-/down-regulated, to filter the logged top-gene tables and labeled
#'   points, to filter which genes are eligible to be outlined as top
#'   down-/up-regulated points, and to draw the horizontal dashed threshold
#'   line at \code{-log10(padj_thresh)}.
#' @param log2fc_thresh Numeric. Log2 fold change magnitude threshold used to
#'   (a) draw the two vertical dashed threshold lines at
#'   \eqn{\pm}\code{log2fc_thresh} and (b) filter which genes are eligible to
#'   be outlined as top down-/up-regulated points.
#' @param xlims Numeric vector of length 2. x-axis (log2 fold change) limits.
#' @param ylims Numeric vector of length 2. y-axis (\code{-log10(p-value)})
#'   limits.
#' @param anno_ypos Numeric. y-coordinate at which the "higher in" annotation
#'   text blocks are placed.
#' @param nudge_x_up Numeric. Horizontal label-repel nudge applied to labels
#'   for the top up-regulated genes.
#' @param nudge_y_up Numeric. Vertical label-repel nudge applied to labels
#'   for the top up-regulated genes.
#' @param nudge_x_dn Numeric. Horizontal label-repel nudge applied to labels
#'   for the top down-regulated genes.
#' @param nudge_y_dn Numeric. Vertical label-repel nudge applied to labels
#'   for the top down-regulated genes.
#' @param plot_height Numeric. Height, in inches, of the saved PDF/PNG (passed
#'   through to \code{\link{plot_pdf_and_png}}).
#' @param plot_width Numeric. Width, in inches, of the saved PDF/PNG (passed
#'   through to \code{\link{plot_pdf_and_png}}).
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects:
#'   logging summary tables of the top genes via \code{\link{log_msg}},
#'   printing the assembled ggplot to the active graphics device, and writing
#'   \code{<plot_name>_<plot_version>.pdf} and \code{.png} to disk via
#'   \code{\link{plot_pdf_and_png}}.
#'
#' @details Requires \pkg{dplyr}, \pkg{ggplot2}, and \pkg{ggrepel} to be
#'   loaded (uses \code{\%>\%}, \code{case_when}, \code{sym}/\code{!!}
#'   tidy-eval, and \code{geom_label_repel}). All significance filtering
#'   (color classification, logged top-gene tables, outlined highlight
#'   points, and labeled points) consistently uses \code{padj_thresh}, so the
#'   highlighted/outlined genes and their labels will always agree.
#'
#' @importFrom dplyr %>% rename filter mutate case_when slice_min slice_max
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal coord_cartesian
#'   geom_hline geom_vline scale_color_manual xlab ylab theme element_text
#'   labs xlim ylim annotate
#' @importFrom ggrepel geom_label_repel
#'
#' @examples
#' \dontrun{
#' plot_volcano(
#'   degs   = "deg_results.txt",
#'   plot_version   = "v1",
#'   plot_dirname   = "output/project",
#'   plot_filebase  = "volcano",
#'   plot_title     = "Condition A vs B",
#'   label_col      = "gene",
#'   log2FC_col     = "log2FoldChange",
#'   pval_col       = "padj",
#'   n_labels       = 10,
#'   name_down      = "Condition A",
#'   name_up        = "Condition B",
#'   anno_down      = "Condition A",
#'   anno_up        = "Condition B",
#'   color_down     = "darkgreen",
#'   color_up       = "darkred",
#'   ns_color       = "gray20",
#'   padj_thresh    = 0.05,
#'   log2fc_thresh  = 1,
#'   xlims          = c(-5, 5),
#'   ylims          = c(0, 15),
#'   anno_ypos      = 12,
#'   nudge_x_up     = 0.5,
#'   nudge_y_up     = 0.5,
#'   nudge_x_dn     = -0.5,
#'   nudge_y_dn     = 0.5,
#'   plot_height    = 4,
#'   plot_width     = 4
#' )
#' }
#'
#' @export
plot_volcano = function(
    degs,
    plot_version,
    plot_dirname,
    plot_filebase,
    plot_title,
    label_col,
    log2FC_col,
    pval_col,
    n_labels,
    name_down,
    name_up,
    anno_down,
    anno_up,
    color_down,
    color_up,
    ns_color,
    padj_thresh,
    log2fc_thresh,
    xlims,
    ylims,
    anno_ypos,
    nudge_x_up,
    nudge_y_up,
    nudge_x_dn,
    nudge_y_dn,
    plot_height,
    plot_width
) {
    # Accept either a file path or an already-loaded data frame/tibble
    if (is.character(degs) && length(degs) == 1) {
        deg_file = read.table(
            file = degs,
            header = TRUE
        ) %>%
        dplyr::rename(gene = !!sym(label_col))
    } else if (is.data.frame(degs)) {
        deg_file = degs
    } else {
        log_msg(
            paste0(
                "`degs` must be either a single file path (character) or a ",
                "data frame/tibble of DEG results, not: ", class(degs)[1]
            ),
            level = 'ERROR'
        )
        stop()
    }
    
    vol_data = deg_file %>%
        # Filter NA/missing padj
        dplyr::filter(!is.na(!!sym(pval_col))) %>%
        mutate(
            color = case_when(
                !!sym(log2FC_col) > 0 & !!sym(pval_col) < padj_thresh ~ name_up,
                !!sym(log2FC_col) < 0 & !!sym(pval_col) < padj_thresh ~ name_down,
                .default = 'n.s.'
            )
        )
  
    log_msg(
        paste0(plot_title, ' top genes for ', anno_down),
        level = 'INFO')
    log_msg(
        vol_data %>%
            filter(!!sym(pval_col) < padj_thresh) %>%
            slice_min(!!sym(log2FC_col), n = n_labels),
        level = 'MSG'
    )
    
    log_msg(paste0(plot_title, ' top genes for ', anno_up))
    log_msg(
        vol_data %>%
            filter(!!sym(pval_col) < padj_thresh) %>%
            slice_max(!!sym(log2FC_col), n = n_labels),
        level = 'MSG'
    )
  
    ggvol = vol_data %>%
        ggplot(
            aes(
                x = !!sym(log2FC_col),
                y = -log10(!!sym(pval_col))
            )
        ) + 
        geom_point(
            aes(colour = color),
            size = 2,
            alpha = 0.9
        ) +
        # Points with outline to highlight the labeled top n genes down
        geom_point(
            data = vol_data %>%
                filter(!!sym(pval_col)   < padj_thresh) %>%
                filter(!!sym(log2FC_col) < -log2fc_thresh) %>%
                slice_min(!!sym(log2FC_col), n = n_labels),
            size = 2,
            shape = 21,
            colour = "black",
            fill = color_down,
            stroke = 1
        ) +

        # Points with outline to highlight the labeled top n genes up
        geom_point(
            data = vol_data %>%
                filter(!!sym(pval_col)   < padj_thresh) %>%
                filter(!!sym(log2FC_col) > log2fc_thresh) %>%
                slice_max(!!sym(log2FC_col), n = n_labels),
            size = 2,
            shape = 21,
            colour = "black",
            fill = color_up,
            stroke = 1
        ) +

        # Plot settings
        theme_minimal() +
        coord_cartesian(ylim = ylims, clip = "off") +

        # Threshold lines
        geom_hline(yintercept = -log10(padj_thresh), lty = 2) +
        geom_vline(xintercept = -log2fc_thresh, lty = 2) +
        geom_vline(xintercept = log2fc_thresh, lty = 2) +
        
        # Color scale for points
        scale_color_manual(
            breaks = c(name_down, name_up),
            values = rlang::list2(
                !!name_down := color_down,
                'n.s.' = ns_color,
                !!name_up := color_up
            )
        ) +
        
        # Title and axis labels
        labs(title = plot_title) +
        xlab(bquote(~log[2]~ Fold ~Change)) +
        ylab(bquote(~-log[10]~'(p-value)')) +
        
        # X and Y limits
        xlim(xlims) +
        ylim(ylims) +

        # Theme adjustments
        theme(
            text = element_text(size = 12),
            legend.position = 'none',
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            plot.title = element_text(hjust = 0.5)
        ) +

        # Annotation for the genes that are down (<0)
        annotate(
            geom = "text",
            x = log2fc_thresh + 0.5,
            y = anno_ypos,
            label = anno_up,
            color = color_up,
            fontface = "bold",
            size = 4,
            vjust = 1,
            hjust = 0,
            lineheight = 0.8
        ) +
        
        # Annotation for the genes that are up (>0)
        annotate(
            geom = "text",
            x = -(log2fc_thresh + 0.5),
            y = anno_ypos,
            label = anno_down,
            color = color_down,
            fontface = "bold",
            size = 4,
            vjust = 1,
            hjust = 1,
            lineheight = 0.8
        ) +        

        # Position nudging for genes that are down
        ggrepel::geom_label_repel(
            data = vol_data %>%
                filter(!!sym(pval_col) < padj_thresh) %>%
                slice_min(!!sym(log2FC_col), n = n_labels),
                aes(label = !!sym(label_col)),
            min.segment.length = 0,
            # arrow = arrow(length = unit(0.01, "npc")),
            size = 3,
            nudge_y = nudge_y_dn,
            nudge_x = nudge_x_dn,
            point.padding = 0.1,
            max.overlaps = 20
        ) +

        # Position nudging for genes that are up
        ggrepel::geom_label_repel(
            data = vol_data %>%
                filter(!!sym(pval_col) < padj_thresh) %>%
                slice_max(!!sym(log2FC_col), n = n_labels),
                aes(label = !!sym(label_col)),
            min.segment.length = 0,
            # arrow = arrow(length = unit(0.01, "npc")),
            size = 3,
            nudge_y = nudge_y_up,
            nudge_x = nudge_x_up,
            point.padding = 0.1,
            max.overlaps = 20
        )
    
    print(ggvol)
    
    dir.create(plot_dirname, recursive = TRUE, showWarnings = FALSE)
    plot_pdf_and_png(
        plot = ggvol,
        filename = file.path(plot_dirname, paste0(plot_filebase, paste0(name_up, '_v_', name_down, '__', plot_version))),
        width = plot_width,
        height = plot_height
    )

    log_msg('Completed!', level = 'DONE')
}


#' Generate a volcano plot from differential expression results (no sourcing of utils.R)
#'
#' Reads a tab/whitespace-delimited differential expression results file,
#' classifies genes as up- or down-regulated based on log2 fold change and
#' the \code{padj_thresh} significance threshold, logs the top \code{n_labels}
#' genes in each direction, and produces an annotated, labeled volcano plot. 
#' The finished plot is displayed and saved to disk as
#' both PDF and PNG via \code{pdf()} and \code{png()}.
#'
#' Genes with a p-value below \code{padj_thresh} and a positive log2 fold
#' change are colored/labeled with \code{name_up}; genes with a p-value
#' below \code{padj_thresh} and a negative log2 fold change are
#' colored/labeled with \code{name_down}; all other genes are marked as
#' non-significant (\code{"n.s."}). The points that are outlined and labeled
#' on the plot are restricted further to genes passing both
#' \code{padj_thresh} and \code{log2fc_thresh}, taking the top 10 most
#' extreme genes in each direction (by log2 fold change) via
#' \code{ggrepel::geom_label_repel}. Dashed threshold lines are drawn at
#' \code{padj_thresh} (horizontal) and at \eqn{\pm}\code{log2fc_thresh}
#' (vertical) for visual reference.
#'
#' @param degs Character. Path to the differential expression
#'   results file. Must be readable by \code{read.table} with
#'   \code{header = TRUE} and must contain a column \code{X} (gene
#'   identifiers), which is renamed to \code{gene}.
#' @param plot_version Character/numeric. Version label appended to the
#'   output filename, e.g. \code{"v1"}.
#' @param plot_dirname Character. Directory name used (together with
#'   \code{plot_version}, \code{name_down}, \code{name_up},
#'   and \code{plot_filebase}) to construct the output PDF/PNG filenames, in the
#'   form \code{<plot_dirname>/<plot_filebase>__<name_up>_v_<name_down>__volcano_<plot_version>.pdf} / \code{.png}.
#' @param plot_filebase Character. String used as the base/root for the 
#'   filename of the plot, e.g. a project name or other informative prefix for the filename.
#' @param plot_title Character. Title displayed on the plot and used in the
#'   log messages announcing the top gene tables.
#' @param label_col Character. Name of the column in \code{degs}
#'   used to label the top genes (e.g. gene symbol).
#' @param log2FC_col Character. Name of the column in \code{degs}
#'   containing log2 fold change values.
#' @param pval_col Character. Name of the column in \code{degs}
#'   containing the p-value or adjusted p-value used for significance
#'   filtering (against \code{padj_thresh}) and for the y-axis
#'   (\code{-log10(pval_col)}).
#' @param n_labels Character. Number of items to label from \code{degs}
#'   in each direction, used to label the top genes (e.g. gene symbol).
#' @param name_down Character. Value used to label/color genes with a
#'   significant, negative log2 fold change (i.e. the \code{color} grouping
#'   value, matched in \code{scale_color_manual}). Compare with
#'   \code{anno_down}, which is the human-readable text shown in the plot
#'   annotation.
#' @param name_up Character. Value used to label/color genes with a
#'   significant, positive log2 fold change (i.e. the \code{color} grouping
#'   value, matched in \code{scale_color_manual}). Compare with
#'   \code{anno_up}, which is the human-readable text shown in the plot
#'   annotation.
#' @param anno_down Character. Label describing the group/condition in which
#'   expression is higher for down-regulated (negative log2FC) genes; used in
#'   the down-regulated annotation text (colored with \code{color_down}) and
#'   in the log messages.
#' @param anno_up Character. Label describing the group/condition in which
#'   expression is higher for up-regulated (positive log2FC) genes; used in
#'   the up-regulated annotation text (colored with \code{color_up}) and in
#'   the log messages.
#' @param color_down Character. Color (name or hex code) used for the
#'   down-regulated group: point fill for the outlined top genes, the
#'   \code{name_down} entry in the color scale, and the down-regulated
#'   annotation text.
#' @param color_up Character. Color (name or hex code) used for the
#'   up-regulated group: point fill for the outlined top genes, the
#'   \code{name_up} entry in the color scale, and the up-regulated
#'   annotation text.
#' @param ns_color Character. Color (name or hex code) used for non-
#'   significant (\code{"n.s."}) points in the color scale.
#' @param padj_thresh Numeric. Adjusted p-value (or p-value) significance
#'   threshold used consistently throughout: to classify/color genes as
#'   up-/down-regulated, to filter the logged top-gene tables and labeled
#'   points, to filter which genes are eligible to be outlined as top
#'   down-/up-regulated points, and to draw the horizontal dashed threshold
#'   line at \code{-log10(padj_thresh)}.
#' @param log2fc_thresh Numeric. Log2 fold change magnitude threshold used to
#'   (a) draw the two vertical dashed threshold lines at
#'   \eqn{\pm}\code{log2fc_thresh} and (b) filter which genes are eligible to
#'   be outlined as top down-/up-regulated points.
#' @param xlims Numeric vector of length 2. x-axis (log2 fold change) limits.
#' @param ylims Numeric vector of length 2. y-axis (\code{-log10(p-value)})
#'   limits.
#' @param anno_ypos Numeric. y-coordinate at which the "higher in" annotation
#'   text blocks are placed.
#' @param nudge_x_up Numeric. Horizontal label-repel nudge applied to labels
#'   for the top up-regulated genes.
#' @param nudge_y_up Numeric. Vertical label-repel nudge applied to labels
#'   for the top up-regulated genes.
#' @param nudge_x_dn Numeric. Horizontal label-repel nudge applied to labels
#'   for the top down-regulated genes.
#' @param nudge_y_dn Numeric. Vertical label-repel nudge applied to labels
#'   for the top down-regulated genes.
#' @param plot_height Numeric. Height, in inches, of the saved PDF/PNG (passed
#'   through to \code{pdf()} and \code{png()}.
#' @param plot_width Numeric. Width, in inches, of the saved PDF/PNG (passed
#'   through to \code{pdf()} and \code{png()}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects:
#'   logging summary tables of the top genes, printing the 
#'   assembled ggplot to the active graphics device, and writing
#'   \code{<plot_name>_<plot_version>.pdf} and \code{.png} to disk via
#'   \code{pdf()} and \code{png()}.
#'
#' @details Requires \pkg{dplyr}, \pkg{ggplot2}, and \pkg{ggrepel} to be
#'   loaded (uses \code{\%>\%}, \code{case_when}, \code{sym}/\code{!!}
#'   tidy-eval, and \code{geom_label_repel}). All significance filtering
#'   (color classification, logged top-gene tables, outlined highlight
#'   points, and labeled points) consistently uses \code{padj_thresh}, so the
#'   highlighted/outlined genes and their labels will always agree.
#'
#' @importFrom dplyr %>% rename filter mutate case_when slice_min slice_max
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal coord_cartesian
#'   geom_hline geom_vline scale_color_manual xlab ylab theme element_text
#'   labs xlim ylim annotate
#' @importFrom ggrepel geom_label_repel
#'
#' @examples
#' \dontrun{
#' plot_volcano(
#'   degs   = "deg_results.txt",
#'   plot_version   = "v1",
#'   plot_dirname   = "output/project",
#'   plot_filebase  = "volcano",
#'   plot_title     = "Condition A vs B",
#'   label_col      = "gene",
#'   log2FC_col     = "log2FoldChange",
#'   pval_col       = "padj",
#'   n_labels       = 10,
#'   name_down      = "Condition A",
#'   name_up        = "Condition B",
#'   anno_down      = "Condition A",
#'   anno_up        = "Condition B",
#'   color_down     = "darkgreen",
#'   color_up       = "darkred",
#'   ns_color       = "gray20",
#'   padj_thresh    = 0.05,
#'   log2fc_thresh  = 1,
#'   xlims          = c(-5, 5),
#'   ylims          = c(0, 15),
#'   anno_ypos      = 12,
#'   nudge_x_up     = 0.5,
#'   nudge_y_up     = 0.5,
#'   nudge_x_dn     = -0.5,
#'   nudge_y_dn     = 0.5,
#'   plot_height    = 4,
#'   plot_width     = 4
#' )
#' }
#'
#' @export
plot_volcano_basic = function(
    degs,
    plot_version,
    plot_dirname,
    plot_filebase,
    plot_title,
    label_col,
    log2FC_col,
    pval_col,
    n_labels,
    name_down,
    name_up,
    anno_down,
    anno_up,
    color_down,
    color_up,
    ns_color,
    padj_thresh,
    log2fc_thresh,
    xlims,
    ylims,
    anno_ypos,
    nudge_x_up,
    nudge_y_up,
    nudge_x_dn,
    nudge_y_dn,
    plot_height,
    plot_width
) {
    # Accept either a file path or an already-loaded data frame/tibble
    if (is.character(degs) && length(degs) == 1) {
        deg_file = read.table(
            file = degs,
            header = TRUE
        ) %>%
        dplyr::rename(gene = !!sym(label_col))
    } else if (is.data.frame(degs)) {
        deg_file = degs
    } else {
        stop(
            paste0(
                "`degs` must be either a single file path (character) or a ",
                "data frame/tibble of DEG results, not: ", class(degs)[1]
            )
        )
    }
    
    vol_data = deg_file %>%
        # Filter NA/missing padj
        dplyr::filter(!is.na(!!sym(pval_col))) %>%
        mutate(
            color = case_when(
                !!sym(log2FC_col) > 0 & !!sym(pval_col) < padj_thresh ~ name_up,
                !!sym(log2FC_col) < 0 & !!sym(pval_col) < padj_thresh ~ name_down,
                .default = 'n.s.'
            )
        )

    message(paste0(plot_title, '\n  ', gsub('\n', ' ', anno_down)))
    print(
        vol_data %>%
            filter(!!sym(pval_col) < padj_thresh) %>%
            slice_min(!!sym(log2FC_col), n = n_labels)
    )
    
    message(paste0(plot_title, '\n  ', gsub('\n', ' ', anno_up)))
    print(
        vol_data %>%
            filter(!!sym(pval_col) < padj_thresh) %>%
            slice_max(!!sym(log2FC_col), n = n_labels)
    )
  
    ggvol = vol_data %>%
        ggplot(
            aes(
                x = !!sym(log2FC_col),
                y = -log10(!!sym(pval_col))
            )
        ) + 
        geom_point(
            aes(colour = color),
            size = 2,
            alpha = 0.9
        ) +
        # Points with outline to highlight the labeled top n genes down
        geom_point(
            data = vol_data %>%
                filter(!!sym(pval_col)   < padj_thresh) %>%
                filter(!!sym(log2FC_col) < -log2fc_thresh) %>%
                slice_min(!!sym(log2FC_col), n = n_labels),
            size = 2,
            shape = 21,
            colour = "black",
            fill = color_down,
            stroke = 1
        ) +

        # Points with outline to highlight the labeled top n genes up
        geom_point(
            data = vol_data %>%
                filter(!!sym(pval_col)   < padj_thresh) %>%
                filter(!!sym(log2FC_col) > log2fc_thresh) %>%
                slice_max(!!sym(log2FC_col), n = n_labels),
            size = 2,
            shape = 21,
            colour = "black",
            fill = color_up,
            stroke = 1
        ) +

        # Plot settings
        theme_minimal() +
        coord_cartesian(ylim = ylims, clip = "off") +

        # Threshold lines
        geom_hline(yintercept = -log10(padj_thresh), lty = 2) +
        geom_vline(xintercept = -log2fc_thresh, lty = 2) +
        geom_vline(xintercept = log2fc_thresh, lty = 2) +
        
        # Color scale for points
        scale_color_manual(
            breaks = c(name_down, name_up),
            values = rlang::list2(
                !!name_down := color_down,
                'n.s.' = ns_color,
                !!name_up := color_up
            )
        ) +
        
        # Title and axis labels
        labs(title = plot_title) +
        xlab(bquote(~log[2]~ Fold ~Change)) +
        ylab(bquote(~-log[10]~'(p-value)')) +
        
        # X and Y limits
        xlim(xlims) +
        ylim(ylims) +

        # Theme adjustments
        theme(
            text = element_text(size = 12),
            legend.position = 'none',
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            plot.title = element_text(hjust = 1)
        ) +

        # Annotation for the genes that are down (<0)
        annotate(
            geom = "text",
            x = log2fc_thresh + 0.5,
            y = anno_ypos,
            label = anno_up,
            color = color_up,
            fontface = "bold",
            size = 4,
            vjust = 1,
            hjust = 0,
            lineheight = 0.8
        ) +
        
        # Annotation for the genes that are up (>0)
        annotate(
            geom = "text",
            x = -(log2fc_thresh + 0.5),
            y = anno_ypos,
            label = anno_down,
            color = color_down,
            fontface = "bold",
            size = 4,
            vjust = 1,
            hjust = 1,
            lineheight = 0.8
        ) +        

        # Position nudging for genes that are down
        ggrepel::geom_label_repel(
            data = vol_data %>%
                filter(!!sym(pval_col) < padj_thresh) %>%
                slice_min(!!sym(log2FC_col), n = n_labels),
                aes(label = !!sym(label_col)),
            min.segment.length = 0,
            # arrow = arrow(length = unit(0.01, "npc")),
            size = 3,
            nudge_y = nudge_y_dn,
            nudge_x = nudge_x_dn,
            point.padding = 0.1,
            max.overlaps = 20
        ) +

        # Position nudging for genes that are up
        ggrepel::geom_label_repel(
            data = vol_data %>%
                filter(!!sym(pval_col) < padj_thresh) %>%
                slice_max(!!sym(log2FC_col), n = n_labels),
                aes(label = !!sym(label_col)),
            min.segment.length = 0,
            # arrow = arrow(length = unit(0.01, "npc")),
            size = 3,
            nudge_y = nudge_y_up,
            nudge_x = nudge_x_up,
            point.padding = 0.1,
            max.overlaps = 20
        )
    
    print(ggvol)
    
    dir.create(plot_dirname, recursive = TRUE, showWarnings = FALSE)
    filename = file.path(plot_dirname, paste0(plot_filebase, paste0(name_up, '_v_', name_down, '__', plot_version)))

    pdf_file = paste0(filename, '.pdf')
    grDevices::pdf(pdf_file, width = plot_width, height = plot_height)
    print(ggvol)
    dev.off()

    message(sprintf("Saved plot to %s", pdf_file))

    # Write PNG
    png_file = paste0(filename, '.png')
    grDevices::png(png_file, width = plot_width, height = plot_height, res = 600, units = 'in')
    print(ggvol)
    dev.off()

    message('Completed!')
}
