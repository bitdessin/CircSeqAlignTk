# nocov start

#' @importFrom stats setNames
.dir_roots <- function() {
    ws_dpath <- getwd()
    home_dpath <- ifelse(Sys.info()[['sysname']] == 'Windows',
                         Sys.getenv('USERPROFILE'),
                         normalizePath(path.expand('~')))
    volume_dpath <- ifelse(Sys.info()[['sysname']] == 'Windows', 'C:/', '/')
    
    return(setNames(
        c(ws_dpath, home_dpath, volume_dpath),
        c(getwd(), normalizePath(home_dpath), normalizePath(volume_dpath))
    ))
}


#' @importFrom shinyjs enable disable
.refresh_buttons_status <- function(enable_button = NULL, rv = NULL) {
    if (is.null(enable_button)) {
        enable_button <- FALSE
        if (!is.null(rv$inputs$fastq) && !is.null(rv$inputs$refseq))
            if (rv$proc_status == 'idle')
                enable_button <- TRUE
    }
    if (enable_button) {
        enable('run_qc')
        enable('run_align')
    } else {
        disable('run_qc')
        disable('run_align')
    }
}


.parse_args <- function(args, arg_prefix) {
    arg_prefix <- paste0('^', arg_prefix, '__')
    str_params <- ''
    for (i in grep(arg_prefix, names(args))) {
        arg_name <- sub(arg_prefix, '', names(args)[i])
        arg_val <- args[[names(args)[i]]]
        if (arg_name == 'params') {
            str_params <- paste0(str_params, ' ', arg_val)
        } else {
            str_params <- paste0(str_params, paste0('--', arg_name, ' ', arg_val), ' ')
        }
    }
    if (arg_prefix == '^qc__')
        str_params <- paste0(str_params,
            ' --minlength ', args[['_qc__read_length_range']][1],
            ' --maxlength ', args[['_qc__read_length_range']][2])
    if (nchar(gsub(' ', '', str_params)) == 0)
        str_params <- NULL
    return(str_params)
}



#' @importFrom tools file_path_sans_ext file_ext
.basename <- function(fpath) {
    fname <- ifelse(file_ext(fpath) %in% c('gz', 'gzip', 'bz2', 'bzip2'),
                file_path_sans_ext(file_path_sans_ext(basename(fpath))),
                file_path_sans_ext(basename(fpath)))
    return(fname)
}


#' @importFrom tools file_ext
.is_compressed_file <- function(fpath) {
    if (file.exists(fpath))
        if (file_ext(fpath) %in% c('gz', 'gzip', 'bz2', 'bzip2'))
            return(TRUE)
    return(FALSE)
}


#' @importFrom R.utils gunzip bunzip2
#' @importFrom tools file_ext
.decompress_file <- function(input_fpath, output_fpath) {
    if (.is_compressed_file(input_fpath)) {
        if (file_ext(input_fpath) %in% c('gz', 'gzip')) {
            gunzip(input_fpath, output_fpath, overwrite = TRUE, remove = FALSE)
        } else if (file_ext(input_fpath) %in% c('bz2', 'bzip2')) {
            bunzip2(input_fpath, output_fpath, overwrite = TRUE, remove = FALSE)
        }
    } else {
        output_fpath <- input_fpath
    }
    return(output_fpath)
}


#' @importFrom R.utils gzip bzip2
#' @importFrom tools file_ext
.compress_file <- function(input_fpath, output_fpath) {
    if (!.is_compressed_file(input_fpath)) {
        if (file_ext(output_fpath) %in% c('gz', 'gzip')) {
            gzip(input_fpath, output_fpath, overwrite = TRUE, remove = TRUE)
        } else if (file_ext(output_fpath) %in% c('bz2', 'bzip2')) {
            bzip2(input_fpath, output_fpath, overwrite = TRUE, remove = TRUE)
        }
    } else {
        output_fpath <- input_fpath
    }
    return(output_fpath)
}


#' @importFrom Rbowtie2 remove_adapters
.run_qc <- function(fastq_fpath, ws_dpath, params = NULL) {
    fastq_fname <- .basename(fastq_fpath)
    clean_fastq_fpath <- file.path(ws_dpath,
                                   paste0(fastq_fname, '.clean.fastq'))
    
    # decompress file to avoid errors in `remove_adaptors` function
    fastq_fpath_ <- .decompress_file(fastq_fpath, tempfile())
    # qc
    remove_adapters(file1 = fastq_fpath_,
                    adapter1 = params$`_qc__adapter1`,
                    adapter2 = NULL,
                    output1 = clean_fastq_fpath,
                    basename = file.path(ws_dpath, 'AdapterRemoval.log'),
                    overwrite = TRUE,
                    .parse_args(params, 'qc'))
    # compressed file if the input is compressed
    if (.is_compressed_file(fastq_fpath)) {
        clean_fastq_fpath <- .compress_file(clean_fastq_fpath,
                        paste0(clean_fastq_fpath, '.', file_ext(fastq_fpath)))
    }
    return(clean_fastq_fpath)
}


.run_alignment <- function(fastq_fpath, refseq_fpath, ws_dpath, params = NULL) {
    ref_index <- build_index(
        input = refseq_fpath,
        output = file.path(ws_dpath, 'ref_index'),
        n_threads = params$`_aln__threads`,
        aligner = tolower(params$`_aln__aligner`)
    )
    aln <- align_reads(
            input = fastq_fpath,
            index = ref_index,
            output = file.path(ws_dpath, 'align_results'),
            n_threads = params$`_aln__threads`,
            n_mismatch = params$`_aln__nmismatches`,
            aligner = tolower(params$`_aln__aligner`),
            add_args = .parse_args(params, 'aln')
    )
    return(calc_coverage(aln))
}

#' @importFrom utils write.table
#' @importFrom ggplot2 ggsave
.save_results <- function(ws_dpath, alncov) {
    output_dpath <- file.path(ws_dpath, 'results')
    dir.create(output_dpath, showWarnings = FALSE)
    write.table(alncov@forward,
                file.path(output_dpath, 'coverage_forward.txt'),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
    write.table(alncov@reverse,
                file.path(output_dpath, 'coverage_reverse.txt'),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
    fdata <- plot(alncov)
    ggsave(file.path(output_dpath, 'coverage_figure.png'),
           fdata,
           width = 10, height = 5, dpi = 300)
}


    

#' @importFrom parallel detectCores
#' @importFrom shiny fluidPage titlePanel wellPanel fluidRow column
#' @importFrom shiny HTML div h1 h2 h3 p hr icon
#' @importFrom shiny textInput sliderInput selectInput actionButton
#' @importFrom shiny verbatimTextOutput
#' @importFrom shiny tableOutput tabsetPanel tabPanel
#' @importFrom shinyFiles shinyDirButton shinyFilesButton
#' @importFrom shinyjs useShinyjs disabled hidden extendShinyjs
#' @importFrom plotly plotlyOutput
#' @importFrom htmltools tags
.ui <- fluidPage(
    useShinyjs(),
    extendShinyjs(text = 'shinyjs.closeWindow = function() { window.close(); }',
                  functions = c("closeWindow")),
    tags$head(tags$script('
        window.onbeforeunload = function() {
            return "Changes you made may not be saved.";
        }')),
    
    fluidRow(
        column(width = 12,
            h1('CircSeqAlignTk',
               style = 'margin-bottom: 1.8rem;'),
            p(paste(
                'CircSeqAlignTk',
                paste0('v', as.character(packageVersion('CircSeqAlignTk'))),
                'is designed as an end-to-end analysis tool,',
                'from alignment to visualisation of small RNA-Seq data of',
                'viroids',
                '(246-401 nt; single-stranded circular non-coding RNAs).'),
              style = 'margin-bottom: 2rem;')
        )
    ),
    
    wellPanel(
        fluidRow(
            column(width = 3,
                h2(
                    'Inputs/Outputs',
                   style = 'margin-top: 0'
                ),
                p(paste('Select a FASTQ of small RNA-seq reads file',
                        'and a FASTA of reference sequence,',
                        'and set the working directory.'))
            ),
            column(width = 9,
                fluidRow(
                    column(width = 4,
                           shinyFilesButton(
                               'fastq_fpath',
                               label = 'FASTQ',
                               title = 'FASTQ',
                               multiple = FALSE,
                               icon = icon('file'),
                               style = 'width:100%'
                           )
                    ),
                    column(width = 8,
                           verbatimTextOutput('fastq_fpath')
                    )
                ),
                
                fluidRow(
                    column(width = 4,
                           shinyFilesButton(
                               'refseq_fpath',
                               label = 'FASTA',
                               title = 'FASTA (viroid reference sequence)',
                               multiple = FALSE,
                               icon = icon('file'),
                               style = 'width:100%'
                           )
                    ),
                    column(width = 8,
                           verbatimTextOutput('refseq_fpath')
                    )
                ),
                
                fluidRow(
                    column(width = 4,
                           shinyDirButton(
                               'workspace_dpath',
                               label = 'Working Directory',
                               title = 'Working Directory',
                               multiple = FALSE,
                               icon = icon('folder'),
                               style = 'width:100%'
                           )
                    ),
                    column(width = 8,
                           verbatimTextOutput('workspace_dpath')
                    )
                ),
            )
        ),
    ),
    
    wellPanel(
        fluidRow(
            column(width = 3,
                h2(
                    'FASTQ QC',
                    style = 'margin-top: 0;'
                ),
                p(HTML(paste('Quality control for FASTQ using',
                '<a href="https://doi.org/10.1186/s13104-016-1900-2"',
                'target="_blank">AdapterRemoval v2</a>',
                'packaged in',
                '<a href="https://www.bioconductor.org/packages/release/bioc/html/Rbowtie2.html"',
                'target="_blank">RBowtie2</a>',
                'package.'))),
                disabled(actionButton(
                    'run_qc',
                    'Run QC',
                    width = '100%',
                    icon = icon('play'),
                    onclick = '$("#run_qc").html("RUNNING ...");'
                ))
            ),
            column(width = 9,
                textInput(
                    '_qc__adapter1',
                    'Adapter sequence expected to be found in reads:',
                    value = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                    placeholder = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
                    width = '100%'
                ),
                sliderInput(
                    '_qc__read_length_range',
                    'Range of read length (reads outside of this range are discarded):',
                    min = 10,
                    max = 60,
                    value = c(21, 24),
                    step = TRUE,
                    width = '100%'
                ),
                sliderInput(
                    'qc__minquality',
                    'Minimum Phred score (bases at 5\'/3\'-end with quality scores less than this are trimmed):',
                    min = 0,
                    max = 40,
                    value = 30,
                    step = TRUE,
                    width = '100%'
                ),
                sliderInput(
                    'qc__maxns',
                    'Maximum number of Ns allowed in a read:',
                    min = 0,
                    max = 5,
                    value = 1,
                    step = TRUE,
                    width = '100%'
                ),
                sliderInput(
                    'qc__threads',
                    'Threads:',
                    min = 0,
                    max = detectCores(),
                    step = TRUE,
                    value =  as.integer(detectCores() / 2)
                ),
                textInput(
                    'qc__params',
                    'Additional arguments to be passed on to the AdapterRemoval:',
                    value = '',
                    placeholder = '--seed 123',
                    width = '100%'
                ),
            )
        ),
        
    ),
    
    wellPanel(
        fluidRow(
            column(width = 3,
                h2(
                    'Alignment',
                    style = 'margin-top: 0'
                ),
                p(HTML('Align reads to reference sequence using <a href="https://bioconductor.org/packages/release/bioc/html/Rbowtie2.html" target="_blank">Bowtie2</a> or <a href="https://bioconductor.org/packages/release/bioc/html/Rhisat2.html" target="_blank">HISAT2</a>.')),
                disabled(actionButton(
                    'run_align',
                    'Run Alignment',
                    width = '100%',
                    icon = icon('play'),
                    onclick = '$("#run_align").html("RUNNING ...");'
                ))
            ),
            column(width = 9,
                selectInput(
                    '_aln__aligner',
                    'Alignment tool',
                    c('Bowtie2', 'HISAT2')
                ),
                sliderInput(
                    '_aln__nmismatches',
                    'Number of mismatches allowed in alignment:',
                    min = 0,
                    max = 5,
                    step = TRUE,
                    value = 0
                ),
                sliderInput(
                    '_aln__threads',
                    'Threads:',
                    min = 0,
                    max = detectCores(),
                    step = TRUE,
                    value = as.integer(detectCores() / 2)
                ),
                textInput(
                    'aln__params',
                    'Additional arguments to be passed on to the alignment tool:',
                    value = '',
                    placeholder = '--seed 123 --quiet',
                    width = '100%'
                ),
            )
        )
    ),
    
    shinyjs::hidden(wellPanel(
        id = 'alignment_results',
        fluidRow(
            column(width = 12,
                h2('Alignment Result'),
                tabsetPanel(
                    type = 'tabs',
                    tabPanel('Figure',
                             h3('Alignment Coverage'),
                             p('A plot showing the alignment coverage of the reads to the reference sequence. The upper and lower directions of the y-axis represent the alignment coverage of the reads aligned in the forward and reverse strands, respectively.'),
                             plotlyOutput('alignment_summary_figure')),
                    tabPanel('Forward Coverage',
                             h3('Coverage of Forward Reads'),
                             p('A matrix containing the alignment coverage of the forward strand reads.'),
                             tableOutput("forward_alncov")),
                    tabPanel('Reversed Coverage',
                             h3('Coverage of Reversed Reads'),
                             p('A matrix containing the alignment coverage of the reverse strand reads.'),
                             tableOutput("reverse_alncov"))
                )
            )
        )
    )),

    fluidRow(
        column(width = 12,
            div(
                actionButton(
                    'shutdown_app',
                    'Shutdown App',
                    icon = icon("power-off"),
                    style = 'color: #C51605; background-color: #FFF5E0; border-color: #d04848;'),
                style = 'text-align: right;'
            )
        )
    ),
    
    fluidRow(
        column(width = 12,
            div(
                hr(),
                style = 'text-align: center;margin-bottom: 5rem;',
                p(paste0('CircSeqAlignTk v',
                         as.character(packageVersion('CircSeqAlignTk')))),
            )
        )
    )
    
)


#' @importFrom shiny reactive reactiveValues observe observeEvent
#' @importFrom shiny renderText renderPrint renderTable stopApp updateActionButton
#' @importFrom shinyFiles shinyDirChoose shinyFileChoose parseDirPath parseFilePaths
#' @importFrom shinyjs showElement js
#' @importFrom plotly ggplotly renderPlotly
.server <- function(input, output, session) {
    
    rv <- reactiveValues(
        app_path = as.character(normalizePath(getwd())),
        workspace = NULL,
        inputs = list(fastq = NULL, refseq = NULL),
        outputs = list(
            clean_fastq = NULL,
            ref_index = NULL,
            aln = NULL,
            alncov = NULL
        ),
        proc_status = 'idle'
    )
    
    # set workspace
    shinyDirChoose(
        input,
        id = 'workspace_dpath',
        roots = .dir_roots(),
        allowDirCreate = TRUE
    )
    workspace_dpath <- reactive({
        if (all(c('path', 'root') %in% names(input$workspace_dpath))) {
            rv$workspace <- as.character(parseDirPath(.dir_roots(),
                               input$workspace_dpath))
        } else {
            rv$workspace <- rv$app_path
        }
        return(rv$workspace)
    })
    output$workspace_dpath <- renderText({
        workspace_dpath()
    })
    
    # set file inputs
    shinyFileChoose(input, 'fastq_fpath', roots = .dir_roots())
    fastq_fpath <- reactive({
        if (all(c('root') %in% names(input$fastq_fpath))) {
            fpath <- parseFilePaths(.dir_roots(),
                           input$fastq_fpath)
            rv$inputs$fastq <- fpath$datapath[[1]]
        } else {
            rv$inputs$fastq <- NULL
        }
        rv$inputs$fastq
    })
    output$fastq_fpath <- renderText({
        if (is.null(fastq_fpath())) {
            'FASTQ (small RNAs sequenced from viroid-infected plants)'
        } else {
            fastq_fpath()
        }
    })
    
    shinyFileChoose(input, 'refseq_fpath', roots = .dir_roots())
    refseq_fpath <- reactive({
        if (all(c('root') %in% names(input$refseq_fpath))) {
            fpath <- parseFilePaths(.dir_roots(),
                           input$refseq_fpath)
            rv$inputs$refseq <- fpath$datapath[[1]]
        } else {
            rv$inputs$refseq <- NULL
        }
        rv$inputs$refseq
    })
    output$refseq_fpath <- renderText({
        if (is.null(refseq_fpath())) {
            'FASTA (reference sequence of a viroid)'
        } else {
            refseq_fpath()
        }
    })
    
    # check inputs and enable buttons
    observe({
        .refresh_buttons_status(NULL, rv)
    })

    # QC
    observeEvent(input$run_qc, {
        .refresh_buttons_status(FALSE, rv)
        rv$proc_status <- 'qc'
        rv$outputs$clean_fastq <- .run_qc(rv$inputs$fastq,
                                          rv$workspace,
                                          input)
        rv$proc_status <- 'idle'
        updateActionButton(session, 'run_qc', label = 'Run QC', icon = icon('play'))
        .refresh_buttons_status(TRUE, rv)
    })
    
    # Quantification
    observeEvent(input$run_align, {
        .refresh_buttons_status(FALSE, rv)
        rv$proc_status <- 'align'
        if (is.null(rv$outputs$clean_fastq))
            rv$outputs$clean_fastq <- rv$inputs$fastq
        rv$outputs$alncov <- .run_alignment(rv$outputs$clean_fastq,
                                rv$inputs$refseq,
                                rv$workspace,
                                input)
        output$forward_alncov <- renderTable({
            data.frame(rv$outputs$alncov@forward)
            })
        output$reverse_alncov <- renderTable({
            data.frame(rv$outputs$alncov@reverse)
            })
        output$alignment_summary_figure <- renderPlotly({
            ggplotly(plot(rv$outputs$alncov))
        })
        .save_results(rv$workspace, rv$outputs$alncov)
        showElement('alignment_results')
        rv$proc_status <- 'idle'
        updateActionButton(session, 'run_align', label = 'Run Alignment', icon = icon('play'))
        .refresh_buttons_status(TRUE, rv)
    })
    

    # shutdown
    observeEvent(input$shutdown_app, {
        stopApp()
    })
    
}


#' Build a GUI application of CircSeqAlignTk
#'
#' Build a graphical user interface (GUI) application for CircSeqAlignTk
#' using Shiny package.
#' 
#' The CircSeqAlignTk graphical user interface (GUI) application
#' is built using the Shiny package. Users need to install the Shiny package
#' and associated packages (shinyFiles, shinyjs) to use this application.
#' Additionally, the installation of RBowtie2 and Rhisat2 is required
#' to use the full functionality of this application.
#' 
#' @param ... Arguments to be passed to \code{\link[shiny]{shinyApp}}.
#' @return A Shiny application object.
#' @examples
#' \dontrun{
#' library(shiny)
#' library(CircSeqAlignTk)
#' 
#' app <- build_app()
#' runApp(app)
#' }
#' @importFrom shiny shinyApp runApp
#' @export
build_app <- function(...) {
    useShinyjs()
    app <- shinyApp(ui = .ui, server = .server, ...)
}

# nocov end
