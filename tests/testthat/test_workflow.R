test_that('workflow with Bowtie2', {
    setwd(tempdir())

    fa <- system.file(package='CircSeqAlignTk', 'extdata', 'FR851463.fa')
    fq <- system.file(package='CircSeqAlignTk', 'extdata', 'srna.fq.gz')
    # index
    ref_index <- build_index(input = fa, output = 'tmp_bt2index')
    # run with default arguments
    aln <- align_reads(input = fq, index = ref_index, output = 'tmp_bt2output')
    # run with additional arguments
    aln <- align_reads(input = fq, index = ref_index, output = 'tmp_bt2output',
                       aligner = 'bowtie2', add_args = '-N 0 -L 22')
})


test_that('workflow with HISAT2', {
    setwd(tempdir())

    fa <- system.file(package='CircSeqAlignTk', 'extdata', 'FR851463.fa')
    fq <- system.file(package='CircSeqAlignTk', 'extdata', 'srna.fq.gz')
    # index
    ref_index <- build_index(input = fa,
                             output = 'tmp_ht2index', aligner = 'hisat2')
    # run with default arguments
    aln <- align_reads(input = fq, index = ref_index, output = 'tmp_ht2output',
                       aligner = 'hisat2')
    # run with additional arguments
    aln <- align_reads(input = fq, index = ref_index, output = 'tmp_ht2output',
                       aligner = 'hisat2',
                       add_args = '--no-spliced-alignment -k 10')
})



test_that('workflow - no reads aligned on edges', {
    setwd(tempdir())

    peaks <- data.frame(mean = c(150, 250),
                        std = c(2, 5),
                        strand = c('+', '-'),
                        prob = c(0.6, 0.4))
    sim <- generate_reads(output = 'tmp_sim_reads_1.fq.gz',
                          adapter = NA,
                          peaks = peaks)

    fa <- system.file(package='CircSeqAlignTk', 'extdata', 'FR851463.fa')
    ref_index <- build_index(input = fa, output = 'tmp_index')
    aln <- align_reads(input = 'tmp_sim_reads_1.fq.gz',
                       index = ref_index, output = 'tmp_sim_output_1')
    alncov <- calc_coverage(aln)
})




test_that('workflow - reads aligned on edges', {
    setwd(tempdir())

    peaks <- data.frame(mean = c(5, 150, 250, 354),
                        std = c(5, 5, 5, 5),
                        strand = c('-', '+', '-', '+'),
                        prob = c(0.3, 0.2, 0.1, 0.4))
    sim <- generate_reads(output = 'tmp_sim_reads_2.fq.gz',
                          adapter = NA,
                          peaks = peaks)

    fa <- system.file(package='CircSeqAlignTk', 'extdata', 'FR851463.fa')
    ref_index <- build_index(input = fa, output = 'tmp_index')
    aln <- align_reads(input = 'tmp_sim_reads_2.fq.gz',
                       index = ref_index, output = 'tmp_sim_output_2')
    alncov <- calc_coverage(aln)
})



test_that('workflow - reads without adapter and mismatches correctly aligned', {
    setwd(tempdir())

    viroid_seq <- system.file(package="CircSeqAlignTk", "extdata", "FR851463.fa")
    ref_index <- build_index(input = viroid_seq, output = 'tmp_index')

    fwd_rmse <- rev_rmse <- rep(NA, 10)
    for (i in seq(fwd_rmse)) {
        # prepare file names and directory to store the simulation results
        dir.create(paste0('tmp_tries_', i))
        syn_fq <- paste0('tmp_tries_', i, '/synthetic_reads.fq.gz')
        align_result <- paste0('tmp_tries_', i, '/align_results')
        fig_coverage <- paste0('tmp_tries_', i, '/alin_coverage.png')

        # generate synthetic reads
        set.seed(i)
        sim <- generate_reads(output = syn_fq, mismatch_prob = 0, adapter = NA)

        # alignment
        aln <- align_reads(input = syn_fq,
                           index = ref_index,
                           output = align_result)

        # calculate alignment coverage
        alncov <- calc_coverage(aln)

        # calculate RMSE
        fwd_pred <- slot(alncov, 'forward')
        fwd_true <- slot(slot(sim, 'coverage'), 'forward')
        fwd_rmse[i] <- sqrt(sum((fwd_pred - fwd_true) ^ 2) / length(fwd_true))
        rev_pred <- slot(alncov, 'reverse')
        rev_true <- slot(slot(sim, 'coverage'), 'reverse')
        rev_rmse[i] <- sqrt(sum((rev_pred - rev_true) ^ 2) / length(rev_true))
    }

    expect_equal(sum(fwd_rmse), 0)
    expect_equal(sum(rev_rmse), 0)
})
