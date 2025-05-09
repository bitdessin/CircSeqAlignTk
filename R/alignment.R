#' Build indexes of reference sequences for alignment
#'
#' This function internally calls Bowtie2 or HISAT2
#' to build indexes of reference sequences for alignment preparation.
#'
#' This function generates two types of reference sequences
#' from a genome and indexes them in preparation for alignment.
#' The type 1 reference sequence is identical to
#' the sequence provided by the `input` argument.
#' The type 2 reference sequence is generated
#' by restoring the type 1 reference sequence to a circular RNA
#' and opening the circle at the position opposite to that of type 1.
#' The type 1 and 2 reference sequences are then saved as FASTA format files,
#' `refseq.t1.fa` and `refseq.t2.fa`, respectively,
#' under the directory specified by the `output` argument.
#' Next, the function builds indexes for `refseq.t1.fa` and `refseq.t2.fa`.
#'
#' Two alignment tools (Bowtie2 and HISAT2) can be specified for
#' building indexes through the aligner argument.
#' This function first attempts
#' to call the specified alignment tool
#' installed on the operation system directly;
#' however, if the tool is not installed,
#' then the function attempts to call
#' \code{\link[Rbowtie2]{bowtie2_build}} or
#' \code{\link[Rhisat2]{hisat2_build}} functions
#' implemented in the Rbowtie2 or Rhisat2 packages for indexing.
#'
#' @param input A path to a FASTA format file containing a reference sequence
#'     of a genome for indexing.
#' @param output A path to a directory for saving the reference sequences and
#'     indexes.
#' @param n_threads Number of threads to use for aligning reads.
#' @param overwrite Overwrite the existing files if `TRUE`.
#' @param aligner A string to specify the alignment for indexing.
#' @param add_args A string of additional arguments to be passed on to the
#'     alignment tool directly (e.g., `--quiet`).
#' @return A \code{\link{CircSeqAlignTkRefIndex-class}} object.
#' @examples
#' output_dpath <- tempdir()
#'
#' genome_seq <- system.file(package="CircSeqAlignTk", "extdata", "FR851463.fa")
#' ref_index <- build_index(input = genome_seq,
#'                          output = file.path(output_dpath, "index"))
#' @seealso \code{\link{CircSeqAlignTkRefIndex-class}}
#' @importFrom methods new
#' @importFrom Biostrings DNAStringSet writeXStringSet
#' @export
build_index <- function(input, output = NULL, n_threads = 1, overwrite = TRUE,
                        aligner = c('hisat2', 'bowtie2'), add_args = NULL) {
    aligner <- match.arg(aligner)
    check_files(input)
    if (is.null(output)) stop('Set `output` to specify a directory',
                              'to save index files.')

    dir.create(output, showWarnings = FALSE, recursive = TRUE)

    # generate metadata object to store sequence and index path
    ref_fasta <- read_single_fasta(input)
    ref <- new('CircSeqAlignTkRefIndex',
                name = ref_fasta$name,
                seq = ref_fasta$seq,
                length = ref_fasta$length,
                fasta = c(file.path(output, 'refseq.t1.fa'),
                          file.path(output, 'refseq.t2.fa')),
                index = c(file.path(output, 'refseq.t1'),
                          file.path(output, 'refseq.t2')),
                cut_loc = floor(ref_fasta$length / 2))

    # generate type 1 and type 2 reference
    check_overwrite(c(ref@fasta, ref@index), overwrite,
                    ext = c('bt2', 'bt2l', 'ht2', 'ht2l'))

    seq_1 <- DNAStringSet(ref@seq)
    names(seq_1) <- paste0(ref@name, '_t1')
    writeXStringSet(seq_1, ref@fasta[1], append = FALSE, format = 'fasta')

    seq_2 <- DNAStringSet(paste0(substr(ref@seq, ref@cut_loc, ref@length),
                          substr(ref@seq, 1, ref@cut_loc - 1)))
    names(seq_2) <- paste0(ref@name, '_t2')
    writeXStringSet(seq_2, ref@fasta[2], append = FALSE, format = 'fasta')

    # create indexes
    .build_index(aligner, add_args, n_threads,
                 ref@fasta[1], ref@index[1], overwrite)
    .build_index(aligner, add_args, n_threads,
                 ref@fasta[2], ref@index[2], overwrite)

    ref
}


#' @importFrom Rbowtie2 bowtie2_build
#' @importFrom Rhisat2  hisat2_build
.build_index <- function(aligner, add_args, n_threads, input, output,
                         overwrite) {
    if (is.null(add_args)) add_args <- ''
    add_args <- paste(add_args, ifelse(aligner == 'bowtie2', '--threads', '-p'),
                      n_threads)
    aligner <- paste0(aligner, '-build')

    if (check_cmd(aligner)) {
        message(aligner, ' found on the system (', Sys.which(aligner), '). ',
                'Use the system installed ', aligner, ' to build index.')
        run_status <- system2(aligner,
                              args = c(add_args, '-f',
                                       shQuote(input), shQuote(output)))
        if (run_status != 0) stop(aligner, 'could not be executed.')
    } else {
        if (aligner == 'bowtie2-build') {
            suppressWarnings(suppressMessages(bowtie2_build(references = input,
                    bt2Index = output, overwrite = overwrite, add_args)))
        } else if (aligner == 'hisat2-build') {
            hisat2_args <- argparse_hisat2(add_args)
            hisat2_args$references <- input
            hisat2_args$outdir <- dirname(output)
            hisat2_args$prefix <- basename(output)
            hisat2_args$force <- overwrite
            do.call('hisat2_build', hisat2_args)
        }
    }
    invisible(output)
}





#' Filter sequence reads in a FASTQ file by length
#'
#' This function removes sequence reads with lengths outside
#' the specified range from the FASTQ file.
#'
#' Studies on small RNA-seq data from viroid-infected plants
#' have mostly focused on reads with lengths ranging from 21 nt to 24 nt.
#' This function is intended to be used to
#' remove sequence reads with lengths outside the specified range.
#' The default range is 21-24 nt,
#' which can be changed through the `read_lengths` argument.
#'
#' Note that, if filtering by read length has already been performed
#' during the quality control process,
#' there is no need to use this function.
#'
#' @param input A path to a FASTQ file targeted for filtering.
#' @param output A path to save the filtered reads in FASTQ format.
#' @param read_lengths A series of integers to specify read length.
#'                     Reads other than the length specified will be excluded
#'                     during alignment.
#' @param overwrite Overwrite the existing files if \code{TRUE}.
#' @return A path to the filtered FASTQ file.
#' @examples
#' output_dpath <- tempdir()
#'
#' fq <- system.file(package="CircSeqAlignTk", "extdata", "srna.fq.gz")
#' output_fq <- file.path(output_dpath, "sran.filtered.fq.gz")
#' filter_reads(fq, output_fq, seq(21, 24))
#' @importFrom BiocGenerics width
#' @importFrom ShortRead FastqStreamer yield writeFastq
#' @export
filter_reads <- function(input, output,
                         read_lengths = seq(21, 24), overwrite = TRUE) {
    check_overwrite(output, overwrite)
    remove_files(output)

    c_mode <- ifelse (file_ext(output) %in% c('gz', 'gzip'), TRUE, FALSE)

    fh <- FastqStreamer(input, n = 1e6)
    while (length(fq <- yield(fh))) {
        if (!is.null(read_lengths)) fq <- fq[width(fq) %in% read_lengths]
        if (length(fq) > 0) writeFastq(fq, output, mode = 'a', full = FALSE,
                                       compress = c_mode)
    }
    close(fh)
    remove_objects(fq)
    invisible(output)
}




#' Align sequence reads to a genome sequence
#'
#' This function aligns sequence reads in a FASTQ file
#' to the reference sequences of a genome.
#'
#' This function aligns sequence reads in a FASTQ format file in two stages:
#' (i) aligning reads to the type 1 reference sequence (i.e., `refseq.t1.fa`)
#' and (ii) collecting the unaligned reads
#' and aligning them with the type 2 reference (i.e., `refseq.t2.fa`).
#' The alignment results are saved as BAM format files
#' in the specified directory with the suffixes `*.t1.bam` and `*.t2.bam`.
#' The original alignment results may contain mismatches.
#' Hence, filtering is performed to remove the alignment with mismatches
#' over the specified value from the BAM format file.
#' The filtered results of the `*.t1.bam` and `*.t2.bam` are saved as
#' `*.clean.t1.bam` and `*.clean.t2.bam`, respectively.
#'
#' Two alignment tools (Bowtie2 and HISAT2) can be specified
#' for building indexes through the aligner argument.
#' This function first attempts to call the specified alignment tool
#' installed on the operation system directly;
#' however, if the tool is not installed,
#' then the function attempts to call
#' \code{\link[Rbowtie2]{bowtie2_build}} or \code{\link[Rhisat2]{hisat2_build}}
#' functions implemented in the Rbowtie2 or Rhisat2 packages for alignment.
#'
#' @param input A path to a FASTQ format file for alignment.
#' @param index A \code{\link{CircSeqAlignTkRefIndex-class}} object generated by
#'     the \code{\link{build_index}} function.
#' @param output A path to a directory for saving the intermediate and final
#'     results of alignment.
#' @param n_threads Number of threads to use for aligning reads.
#' @param n_mismatch Number of allowed mismatches in alignment.
#' @param overwrite Overwrite the existing files if \code{TRUE}.
#' @param aligner A string to specify the alignment is for alignment.
#' @param add_args A string of additional arguments to be passed on to the
#'     alignment tool directly.
#'     For example, \code{-N 0 -L 22}, \code{--no-spliced-alignment -k 10}, etc.
#' @return A \code{\link{CircSeqAlignTkAlign-class}} object.
#' @seealso \code{\link{CircSeqAlignTkAlign-class}}
#' @examples
#' output_dpath <- tempdir()
#'
#' genome_seq <- system.file(package="CircSeqAlignTk", "extdata", "FR851463.fa")
#' fq <- system.file(package="CircSeqAlignTk", "extdata", "srna.fq.gz")
#'
#' ref_index <- build_index(input = genome_seq,
#'                          output = file.path(output_dpath, 'index'))
#' aln <- align_reads(input = fq, index = ref_index,
#'                    output = file.path(output_dpath, 'align_results'))
#'
#' slot(aln, 'stats')
#' @importFrom methods new slot
#' @export
align_reads <- function(input, index, output,
                        n_threads = 1, n_mismatch = 1,
                        overwrite = TRUE,
                        aligner = c('hisat2', 'bowtie2'),
                        add_args = NULL) {
    aligner <- match.arg(aligner)

    # check FASTA
    check_files(input)

    # preparation (create output directory for storing alignment results)
    dir.create(output, showWarnings = FALSE, recursive = TRUE)
    output_prefix <- file.path(output, filename(input))

    # craete a object to store the file path
    aln <- new('CircSeqAlignTkAlign',
               input_fastq = input,
               fastq = c(input, paste0(output_prefix, '.t2.fq.gz')),
               bam = paste0(output_prefix,
                            c('.t1.bam', '.t2.bam')),
               clean_bam = paste0(output_prefix,
                                  c('.clean.t1.bam', '.clean.t2.bam')),
               stats = data.frame(),
               reference = index)
    check_overwrite(c(aln@fastq, aln@bam, aln@clean_bam), overwrite)
    remove_files(c(aln@fastq[2], aln@bam, aln@clean_bam))

    # two-stage alignment
    message('Started align sequence reads onto the type 1 reference sequence ...')
    .align_reads(aligner, add_args, index@index[1], aln@bam[1], aln@fastq[1],
                 n_threads, overwrite)
    message('Collecting the unaligned sequence reads ...')
    .save_unmapped_reads(aln@bam[1], aln@fastq[2])
    message('Started align sequence reads onto the type 2 reference sequence ...')
    .align_reads(aligner, add_args, index@index[2], aln@bam[2], aln@fastq[2],
                 n_threads, overwrite)
    message('Finilize the process ...')
    .clean_bam(aln@bam, aln@clean_bam, n_mismatch)

    # check mapping results
    aln@stats <- .bam2stats(c(aln@bam[1], aln@bam[2],
                              aln@clean_bam[1], aln@clean_bam[2]))

    aln
}


#' @importFrom tools file_path_sans_ext
#' @importFrom Rbowtie2 bowtie2
#' @importFrom Rhisat2 hisat2
#' @importFrom Rsamtools sortBam indexBam asBam
.align_reads <- function(aligner, add_args,
                         index, output, input, n_threads, overwrite) {
    dest <- NULL
    output_sam <- paste0(tempfile(), '.sam')
    
    if (is.null(add_args)) add_args <- ''

    if (check_cmd(aligner)) {
        message(aligner, ' found on the system (', Sys.which(aligner), '). ',
                'Use the system installed ', aligner, ' to align reads.\n')
        # Bowtie2 and HISAT2 has the same usages
        add_args <- paste(add_args, '--threads', n_threads,
                          '-x', shQuote(index), '-U', input,
                          '-S', shQuote(output_sam))
        run_status <- system2(aligner, args = add_args, wait = TRUE)
        if (run_status != 0) stop(aligner, 'could not be executed.\n')
    } else {
        if (aligner == 'bowtie2') {
            add_args <- paste(paste('--threads', n_threads), ' ', add_args)
            suppressMessages(suppressWarnings(bowtie2(bt2Index = index,
                    samOutput = output_sam, seq1 = input, seq2 = NULL,
                    add_args, overwrite = overwrite)))
        } else if (aligner == 'hisat2') {
            decompressed_input <- unzip_fastq(input)
            hisat2_args <- argparse_hisat2(add_args)
            hisat2_args$sequences <- decompressed_input
            hisat2_args$index <- index
            hisat2_args$outfile <- output_sam
            hisat2_args$type <- 'single'
            hisat2_args$force <- overwrite
            hisat2_args$threads <- n_threads
            
            do.call('hisat2', hisat2_args)
            remove_files(decompressed_input)
        }
    }
    
    dest <- asBam(output_sam)
    sortBam(dest, file_path_sans_ext(output))
    indexBam(output)

    remove_files(c(output_sam, dest))

    invisible(output)
}


#' @importFrom Biostrings BStringSet
#' @importFrom ShortRead ShortReadQ
#' @importFrom Rsamtools ScanBamParam scanBam
.save_unmapped_reads <- function(bam_file, fastq_file) {
    scan_bam_params <- ScanBamParam(what = c('qname', 'seq', 'qual', 'flag'))
    bam <- scanBam(bam_file, param = scan_bam_params)[[1]]
    unmapped_reads <- ShortReadQ(bam$seq, bam$qual, BStringSet(bam$qname))
    unmapped_reads <- unmapped_reads[bam$flag == 4]
    remove_files(fastq_file)
    writeFastq(unmapped_reads, fastq_file, mode = 'w', full = FALSE,
               compress = TRUE)
}


#' @importFrom tools file_path_sans_ext
#' @importFrom Rsamtools ScanBamParam filterBam sortBam indexBam BamFile
#' @importFrom S4Vectors FilterRules
.clean_bam <- function(bam_files, clean_bam_files, n_mismatch) {
    scan_bam_params <- ScanBamParam(tag = c('NM'))
    filter_bam_params <- FilterRules(list(nm = function(x) x$NM <= n_mismatch))

    .check_bam <- function(bam_file) {
        f <- scanBam(bam_file, param = ScanBamParam(what = c('flag')))[[1]]$flag
        c(length(f), sum(f != 4))
    }

    dest <- c(NA, NA)
    for (i in seq(bam_files)) {
        bam_reads <- .check_bam(bam_files[i])
        if (bam_reads[2] > 0) {
            temp_file <- paste0(tempfile(), '.bam')
            bam_contents <- BamFile(bam_files[i], yieldSize = bam_reads[1])
            dest[i] <- filterBam(bam_contents, temp_file,
                                 param = scan_bam_params,
                                 filter = filter_bam_params)
            dest[i] <- sortBam(temp_file,
                               file_path_sans_ext(clean_bam_files[i]))
            indexBam(clean_bam_files[i])
            remove_files(temp_file)
        }
    }

    dest
}


#' @importFrom Rsamtools ScanBamParam scanBam
.bam2stats <- function(bam_files) {
    bam_contents <- vector('list', length = length(bam_files))
    scan_bam_params <- ScanBamParam(what = c('flag'))
    for (i in seq(bam_contents)) {
        if (file.exists(bam_files[i])) {
            bam_contents[[i]] <- scanBam(bam_files[i],
                                         aram = scan_bam_params)[[1]]
        }
    }

    stats <- data.frame(
        n_reads = vapply(bam_contents,
                    function(x) {ifelse(is.null(x), NA, length(x[[1]]))},
                    numeric(1)),
        aligned_fwd = vapply(bam_contents,
                    function(x) {ifelse(is.null(x), NA, sum(x$flag == 0))},
                    numeric(1)),
        aligned_rev = vapply(bam_contents,
                    function(x) {ifelse(is.null(x), NA, sum(x$flag == 16))},
                    numeric(1)),
        unaligned = vapply(bam_contents,
                    function(x) {ifelse(is.null(x), NA, sum(x$flag == 4))},
                    numeric(1)),
        unsorted_reads = vapply(bam_contents,
                    function(x) {ifelse(is.null(x), NA,
                                        sum(!(x$flag %in% c(0, 4, 16))))},
                    numeric(1))
    )
    rownames(stats) <- basename(bam_files)
    stats
}





argparse_hisat2 <- function (args_string) {
    arg_list <- NULL
    if (!(is.null(args_string) && args_string == '')) {
        arg_nv <- strsplit(args_string, ' -')[[1]]
        arg_nv <- arg_nv[arg_nv != '']
        arg_nv <- regmatches(arg_nv, regexpr(' ', arg_nv), invert = TRUE)
        arg_n <- gsub('^\\-+', '', vapply(arg_nv, function(x) x[1], ''))
        arg_v <- vapply(arg_nv, function(x) x[2], '')
        names(arg_v) <- arg_n
        arg_list <- as.list(arg_v)
        for (i in seq(arg_list)) {
            if (is.na(arg_list[[i]])||arg_list[[i]] == '') arg_list[[i]] <- TRUE
        }
    }
    arg_list
}


