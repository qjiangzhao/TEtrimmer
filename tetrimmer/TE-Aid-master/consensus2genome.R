#! /bin/Rscript
#######################################################################################
### consensus2genome - v2 - Clement Goubert (2020) - goubert.clement@gmail.com      ###
### ------------------------------------------------------------------------------- ###
### This R function blast a TE consensus against a reference genome and then plots  ###
### the genomic fragments found relative to the consensus sequence                  ###
### see https://github.com/clemgoub/consensus2genome for the full documentation     ###
### to use, copy and paste the following code into a R console                      ###
### USAGE: consensus2genome(query, db , ...)                                        ###
### query: path to query (fasta file)                                               ###
### db: path to blast db (blast formated nucleotide database)                       ###
#######################################################################################

consensus2genome <- function(query = NULL,
                             db = NULL,
                             evalue = 10e-8,
                             FL_thresh = 0.9,
                             alpha = 0.5,
                             full_alpha = 1,
                             auto_y = TRUE,
                             bins = NULL,
                             output = NULL,
                             v_x_line_1 = 0,
                             v_x_line_2 = 0) {
  
  if (is.null(query)) { print('query not specified') }
  if (is.null(db)) { print('db not specified') }
  
  # Perform the blast
  blast <- read.table(
    text = system(
      paste(
        "blastn -max_target_seqs 10000 -query", query, "-db", db, "-evalue", evalue,
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue qcovhsp' | sed 's/#/-/g'"
      ),
      intern = TRUE
    )
  )
  
  # Write BLAST results to the output directory
  output_filepath <- file.path(output, "blastn.txt")
  write.table(blast, file = output_filepath, quote = FALSE, row.names = FALSE)
  
  # TE consensus size
  cons_len <- as.numeric(system(paste(bins, "/getlength.sh ", query, sep = ""), intern = TRUE))
  print(paste("consensus length: ", cons_len, "bp", sep = " "))
  
  if (v_x_line_1 > 0 && v_x_line_2 > 0) {
    real_cons_len <- abs(as.numeric(v_x_line_2) - as.numeric(v_x_line_1))
  } else {
    real_cons_len <- cons_len
  }
  
  # List half full of length fragments
  half_full <- blast[0.5 * as.numeric(real_cons_len) <= abs(blast$V7 - blast$V8) & 
                       abs(blast$V7 - blast$V8) < FL_thresh * as.numeric(real_cons_len), ]
  
  # List of almost full length fragments
  full <- blast[abs(blast$V7 - blast$V8) >= FL_thresh * as.numeric(real_cons_len), ]
  
  # List blast hit for grey color
  blast_grey <- blast[abs(blast$V7 - blast$V8) < 0.5 * as.numeric(real_cons_len), ]
  
  # Graph
  if (auto_y == TRUE) {
    plot(
      range(0, cons_len),
      range(0, max(100 - blast$V3)),
      type = "n",
      main = paste(
        "TE: ", as.character(blast[1, 1]),
        "\n size: ", as.character(cons_len),
        "bp; fragments: ", as.character(length(blast$V1)),
        "; full length: ", as.character(length(full$V1)),
        " (>=", as.character(as.numeric(FL_thresh) * cons_len), "bp)",
        sep = ""
      ),
      cex.main = 2,
      xlab = "TE consensus (bp)",
      ylab = "divergence to consensus (%)",
      xaxs = "i",
      yaxs = "i",
      cex.lab = 2,
      cex.axis = 1.5
    )
    
    # Add vertical blue lines
    if (v_x_line_1 > 0) { abline(v = v_x_line_1, col = "blue", lwd = 2) }
    if (v_x_line_2 > 0) { abline(v = v_x_line_2, col = "blue", lwd = 2) }
    
    for (i in 1:length(blast$V1)) {
      segments(blast_grey$V7[i], 100 - blast_grey$V3[i], blast_grey$V8[i], 100 - blast_grey$V3[i], col = rgb(0, 0, 0, alpha = alpha))
    }
    
    for (i in 1:length(blast$V1)) {
      segments(half_full$V7[i], 100 - half_full$V3[i], half_full$V8[i], 100 - half_full$V3[i], col = rgb(0, 0.7, 0.1, alpha = 0.8), lwd = 1.1)
    }
    
    for (i in 1:length(blast$V1)) {
      segments(full$V7[i], 100 - full$V3[i], full$V8[i], 100 - full$V3[i], col = rgb(1, 0, 0, alpha = full_alpha), lwd = 1.3)
    }
    
  } else {
    plot(
      range(0, cons_len),
      range(0, auto_y),
      type = "n",
      main = paste(
        "TE: ", as.character(blast[1, 1]),
        "\n size: ", as.character(cons_len),
        "bp; fragments: ", as.character(length(blast$V1)),
        "; full length: ", as.character(length(full$V1)),
        " (>=", as.character(as.numeric(FL_thresh) * cons_len), "bp)",
        sep = ""
      ),
      cex.main = 2,
      xlab = "TE consensus (bp)",
      ylab = "divergence to consensus (%)",
      cex.lab = 2,
      cex.axis = 1.5
    )
    
    for (i in 1:length(blast$V1)) {
      segments(blast$V7[i], 100 - blast$V3[i], blast$V8[i], 100 - blast$V3[i], col = rgb(0, 0, 0, alpha = alpha))
    }
    
    for (i in 1:length(blast$V1)) {
      segments(full$V7[i], 100 - full$V3[i], full$V8[i], 100 - full$V3[i], col = rgb(1, 0, 0, alpha = full_alpha), lwd = 1.5)
    }
  }
  
  # --- Optimized Coverage Calculation ---
  cons_len_num <- as.numeric(cons_len)
  cov_vector <- rep(0, cons_len_num)
  
  for (i in 1:nrow(blast)) {
    start <- max(1, min(blast$V7[i], blast$V8[i]))
    end <- min(cons_len_num, max(blast$V7[i], blast$V8[i]))
    
    if (start <= end) {
      cov_vector[start:end] <- cov_vector[start:end] + 1
    }
  }
  
  # --- Updated removator function ---
  removator <- function(covV) {
    plot(
      covV,
      type = "l",
      main = "",
      xaxs = "i",
      yaxs = "i",
      xlab = "TE consensus genomic coverage plot (bp)",
      ylab = "coverage (bp)",
      cex.lab = 2,
      cex.axis = 1.5,
      ylim = c(0, max(covV, na.rm = TRUE) * 1.1)
    )
    
    if (as.numeric(v_x_line_1) > 0) { abline(v = as.numeric(v_x_line_1), col = "blue", lwd = 2) }
    if (as.numeric(v_x_line_2) > 0) { abline(v = as.numeric(v_x_line_2), col = "blue", lwd = 2) }
  }
  
  # --- Execution ---
  removator(cov_vector)
}