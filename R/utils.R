abspath <- function(fpath) {
    file.path(normalizePath(dirname(fpath)), basename(fpath))
}


#' Get the slot contents from a formal class
#' 
#' This function returns the slot contents from a formal class.
#' It is convenient to use `@` when accessing the contents of a slot,
#' however, using `@` will generate warnings
#' during the unit tests under software development.
#' This function was created to avoid that warning.
#' Users do not have to use this function.
#' 
#' @param object An object from a formally defined class.
#' @param name The name of the slot.
#' @return The contents of the specified slot from the given object.
#' @examples
#' output_dpath <- tempdir()
#' sim <- generate_reads(output = file.path(output_dpath, 'sample1.fq.gz'))
#' head(get_slot_contents(sim, 'peak'))
#' @importFrom methods slot
#' @export
get_slot_contents <- function(object, name) {
    slot(object, name)
}


#' @importFrom tools file_path_sans_ext file_ext
filename <- function(fpath) {
    fnames <- NULL
    for (i in seq(fpath)) {
        fname <- basename(fpath[i])
        if (file_ext(fname) %in% c('gz', 'gzip')) {
            fname <- file_path_sans_ext(fname)
        }
        if (file_ext(fname) %in% c('fa', 'fasta', 'fq', 'fastq')) {
            fname <- file_path_sans_ext(fname)
        }
        fnames <- c(fnames, fname)
    }
    fnames
}


check_files <- function(fpath) {
    for (i in seq(fpath)) {
        if (!file.exists(fpath[i])) {
            stop('File "', fpath[i], '" not found. ',
                 'Check the file path and run the function again.\n')
        }
    }
}


remove_files <- function(fpath) {
    for (i in seq(fpath)) {
        if (!is.na(fpath[i]) && file.exists(fpath[i])) file.remove(fpath[i])
    }
}


#' @importFrom tools file_ext
check_overwrite <- function(fpath, overwrite, ext = NULL) {
    existing_files <- NULL
    for (i in seq(fpath)) {
        if (file.exists(fpath[i])) {
            existing_files <- c(existing_files, abspath(fpath[i]))
        }
        if (!is.null(ext)) {
            for (f in list.files(dirname(fpath[i]), full.names = TRUE)) {
                m <- length(grep(basename(fpath[i]),
                                 file_path_sans_ext(basename(f))))
                if ((m > 0) && (file_ext(f) %in% ext)) {
                    existing_files <- c(existing_files, abspath(f))
                }
            }
        }
    }
    if (length(existing_files) > 0) {
        msg <- c('The following files are found.\n')
        for (i in seq(existing_files)) {
            msg <- c(msg, '  - ', existing_files[i], '\n')
        }
        if (overwrite) {
            warning(msg, 'These files will be overwrite.\n')
        } else {
            stop(msg, 'Move these files to other location',
                 ' or use `overwrite = TRUE` argument to overwrite.\n')
        }
    }
}


remove_objects <- function(...) {
    for (x in list(...)) {
        rm(x)
    }
    gc(verbose = FALSE, reset = TRUE)
}


check_cmd <- function(cmd) {
    if (file.exists(Sys.which(cmd))) {
        return (TRUE)
    } else {
        return (FALSE)
    }
}


unzip_fastq <- function(input, output = NULL) {
    if (is.null(output)) {
        output <- paste0(tempfile(), '.fq')
    }
    fh <- FastqStreamer(input, n = 1e6)
    while (length(fq <- yield(fh))) {
        writeFastq(fq, output, mode = 'a', full = FALSE, compress = FALSE)
    }
    close(fh)
    invisible(output)
}


#' @importFrom Biostrings readDNAStringSet
read_single_fasta <- function(fpath) {
    seq <- readDNAStringSet(fpath, 'fasta')
    if (length(seq) > 1) {
        stop('The sequence should only contain a single sequence.',
             'Multiple sequences are found.\n')
    }
    list(name = names(seq)[1], seq = toString(seq[[1]]), length = width(seq)[1])
}


#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate arrange
reshape_coverage_df <- function(x, read_strand) {
    data.frame(position = seq(1, nrow(x)), x) %>%
        pivot_longer(!.data$position,
                     names_to = 'read_length',
                     values_to = 'coverage') %>%
        mutate(read_length = factor(sub('L', '', .data$read_length),
                                    levels = sub('L', '',
                                        sort(unique(.data$read_length))))) %>%
        arrange(.data$read_length) %>%
        mutate(strand = read_strand) %>%
        as.data.frame()
}

