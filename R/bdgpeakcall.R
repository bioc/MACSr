#' bdgpeakcall
#'
#' Call peaks from bedGraph output. Note: All regions on the same
#' chromosome in the bedGraph file should be continuous so only
#' bedGraph files from MACS3 are accpetable.
#'
#' @param ifile MACS score in bedGraph. REQUIRED.
#' @param cutoff Cutoff depends on which method you used for score
#'     track. If the file contains pvalue scores from MACS3, score 5
#'     means pvalue 1e-5. DEFAULT: 5", default = 5.
#' @param minlen minimum length of peak, better to set it as d
#'     value. DEFAULT: 200", default = 200.
#' @param maxgap maximum gap between significant points in a peak,
#'     better to set it as tag size. DEFAULT: 30", default = 30.
#' @param call_summits If set, MACS will use a more sophisticated
#'     approach to find all summits in each enriched peak region
#'     DEFAULT: False",default=False.
#' @param cutoff_analysis While set, bdgpeakcall will analyze number
#'     or total length of peaks that can be called by different cutoff
#'     then output a summary table to help user decide a better
#'     cutoff. Note, minlen and maxgap may affect the
#'     results. DEFAULT: False", default = False.
#' @param trackline Tells MACS not to include trackline with bedGraph
#'     files. The trackline is required by UCSC.
#' @param outputfile The output file.
#' @param outdir The output directory.
#' @param log Whether to capture logs.
#' @param verbose Set verbose level of runtime message. 0: only show
#'     critical message, 1: show additional warning message, 2: show
#'     process information, 3: show debug messages. DEFAULT:2
#' @return `macsList` object.
#' @export
#' @examples
#' \donttest{
#' eh <- ExperimentHub::ExperimentHub()
#' CHIP <- eh[["EH4558"]]
#' CTRL <- eh[["EH4563"]]
#' p1 <- pileup(CHIP, outdir = tempdir(),
#'              outputfile = "pileup_ChIP_bed.bdg", format = "BED")
#' p2 <- pileup(CTRL, outdir = tempdir(),
#'              outputfile = "pileup_CTRL_bed.bdg", format = "BED")
#' c1 <- bdgcmp(p1$outputs, p2$outputs, outdir = tempdir(),
#'              oprefix = "bdgcmp", pseudocount = 1, method = "FE")
#' bdgpeakcall(c1$outputs, cutoff = 2,
#'             outdir = tempdir(), outputfile = "bdgpeakcall")
#' }
bdgpeakcall <- function(ifile, cutoff = 5, minlen = 200L, maxgap = 30L,
                        call_summits = FALSE, cutoff_analysis = FALSE,
                        trackline = TRUE, outdir = ".",
                        outputfile = character(),
                        log = TRUE, verbose = 2L){
    cl <- basiliskStart(env_macs)
    ifile = normalizePath(ifile)
    on.exit(basiliskStop(cl))
    res <- basiliskRun(cl, function(.namespace, outdir){
        opts <- .namespace()$Namespace(ifile = ifile,
                                       cutoff = cutoff,
                                       minlen = minlen,
                                       maxgap = maxgap,
                                       call_summits = call_summits,
                                       cutoff_analysis = cutoff_analysis,
                                       trackline = trackline,
                                       outdir = outdir,
                                       ofile = outputfile,
                                       verbose = verbose)
        .bdgpeakcall <- reticulate::import("MACS3.Commands.bdgpeakcall_cmd")
        if(log){
            reticulate::py_capture_output(.bdgpeakcall$run(opts))
        }else{
            .bdgpeakcall$run(opts)
        }
    }, .namespace = .namespace, outdir = outdir)
    if(log){
        message(res)
    }
    ofile <- file.path(outdir, outputfile)
    args <- as.list(match.call())
    macsList(arguments = args, outputs = ofile, log = res)
}
