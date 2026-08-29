#' Save a plot as both PDF and PNG
#'
#' Writes the same plot object to disk twice: once as a vector PDF and once
#' as a 600 dpi PNG, using a shared filename stem. Each save is logged via
#' \code{\link{log_msg}}, and if either device write fails the error is
#' logged and re-thrown so the calling code can handle or halt on it.
#'
#' @param plot A plot object (e.g. a ggplot object) to be printed and saved.
#' @param filename Character. File path/stem (without extension) used for
#'   both outputs. Results in \code{<filename>.pdf} and \code{<filename>.png}.
#' @param height Numeric. Height of the output image, in inches.
#' @param width Numeric. Width of the output image, in inches.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effects of
#'   writing \code{<filename>.pdf} and \code{<filename>.png} to disk and
#'   logging progress/errors via \code{\link{log_msg}}.
#'
#' @details If either the PDF or PNG device fails to write (e.g. an invalid
#'   \code{filename} path), the underlying error is logged at the
#'   \code{"ERROR"} level via \code{\link{log_msg}} and then re-raised with
#'   \code{stop()}, halting execution.
#'
#' @examples
#' \dontrun{
#' p = ggplot2::ggplot(mtcars, ggplot2::aes(mpg, hp)) + ggplot2::geom_point()
#' plot_pdf_and_png(p, "mtcars_scatter", height = 4, width = 6)
#' }
#'
#' @export
plot_pdf_and_png = function(plot, filename, height, width) {
  tryCatch({
    # Write PDF
    pdf_file = paste0(filename, '.pdf')
    grDevices::pdf(pdf_file, width = width, height = height)
    print(plot)
    dev.off()

    log_msg(sprintf("Saved plot to %s", pdf_file), level = "PLOT")

    # Write PNG
    png_file = paste0(filename, '.png')
    grDevices::png(png_file, width = width, height = height, res = 600, units = 'in')
    print(plot)
    dev.off()

    log_msg(sprintf("Saved plot to %s", png_file), level = "PLOT")
  },
    error = function(e) {
      log_msg(sprintf("Could not write plot: %s", filename), level = "ERROR")
      log_msg(e$message, level = "ERROR")
      stop(e)
    }
  )
}

#' Log a message with a timestamp and level
#'
#' Prints a formatted, color-coded log message to the console, prefixed with
#' a timestamp and a right-aligned log level tag. Colors are applied via ANSI
#' escape codes, so output is colorized in terminals that support ANSI (RStudio
#' console, most Unix terminals); consoles without ANSI support will show the
#' raw escape codes.
#'
#' @param msg The message to log. Typically a single character string, but
#'   any object can be passed (e.g. a data frame or tibble): non-scalar-
#'   character input is first rendered with \code{\link[utils]{capture.output}}
#'   and collapsed into one multi-line string, so it is logged as a single
#'   block under one timestamp/level prefix rather than one line per element.
#'   See \strong{Details}.
#' @param level Character. The log level/tag to display, e.g. \code{"INFO"},
#'   \code{"START"}, \code{"PLOT"}, \code{"WARN"}, \code{"ERROR"},
#'   \code{"DEBUG"}, \code{"MSG"}, or \code{"DONE"}. Matched case-insensitively
#'   (internally converted to uppercase). Levels not in the predefined color
#'   map fall back to plain white. Default is \code{"INFO"}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   printing a formatted line (or block) to the console.
#'
#' @details Timestamps are formatted as \code{"%Y-%m-%d %H:%M:%S"} using the
#'   system time (\code{\link[base]{Sys.time}}). This function has no
#'   external dependencies beyond base R.
#'
#'   Because the underlying \code{\link[base]{sprintf}} call is vectorized,
#'   passing a multi-element object (such as a multi-row data frame) as
#'   \code{msg} directly would otherwise produce one timestamped line per
#'   element/row rather than a single readable message. To avoid this,
#'   any \code{msg} that is not a length-1 character string is first passed
#'   through \code{print()} and \code{\link[utils]{capture.output}}, and the
#'   resulting lines are joined with \code{"\\n"} into a single string before
#'   formatting — so e.g. a data frame prints as one aligned table under a
#'   single timestamp/level prefix, matching what you'd see from
#'   \code{print()} at the console.
#'
#' @examples
#' log_msg("Starting analysis")
#' log_msg("File not found", level = "ERROR")
#' log_msg(head(mtcars), level = "DEBUG")
#'
#' @importFrom utils capture.output
#' @export
log_msg = function(msg, level = "INFO") {
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # ANSI color codes
  colors = c(
    "START" = "\033[33m",
    "INFO"  = "\033[97m",
    "PLOT"  = "\033[92m",
    "WARN"  = "\033[38;5;208m",
    "ERROR" = "\033[31m",
    "DEBUG" = "\033[34m",
    "MSG"   = "\033[90m",
    "DONE"  = "\033[32m"
  )
  
  # reset color
  reset = "\033[0m"  

  level = toupper(level)
  level_fmt = sprintf("%5s", level)
  
  # Pick color for the level, default to white
  color = colors[[level]]
  if (is.null(color)) color = "\033[37m"

  # Collapse non-scalar-character messages (e.g. data frames/tibbles) into
  # a single formatted block, so sprintf() below treats the whole thing as
  # ONE %s substitution instead of recycling the format string across every
  # element (which is what produces one timestamp per row).
  if (!is.character(msg) || length(msg) != 1) {
    msg = paste(utils::capture.output(print(msg)), collapse = "\n")
  }
  
  # Print colored message
  cat(sprintf("%s[%s] [%s] %s%s\n", color, timestamp, level_fmt, msg, reset))
}
