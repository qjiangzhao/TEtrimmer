# self sequence db is being made by the shell script // or you need to make it to use in R
blastdotplot <- function(query = NULL,
                         db = NULL,
                         blast = NULL,
                         os = NULL,
                         tables = NULL,
                         output = NULL) {
  
  # Run the selfblast
  bl <- read.table(
    text = system(
      paste(
        "blastn -max_target_seqs 10000 -query", query,
        "-db", db, "-evalue 0.05 -outfmt 6 -word_size 11",
        "-gapopen 5 -gapextend 2 -reward 2 -penalty -3 | cut -f 1,7-10 | sed 's/#/-/g'"
      ),
      intern = TRUE
    )
  )
  
  # Order from left to right
  bl <- bl[order(bl$V2, decreasing = FALSE), ]

  # ONLY take the top 17 hits
  if (nrow(bl) > 17) {
    bl <- bl[1:17, ]
  }
  
  # Test if there are ORF detected; will later store in orf if TRUE
  if (exists("blast") && file.exists(as.character(blast))) {
    test <- try(read.table(as.character(blast)), silent = TRUE)
  } else {
    test <- FALSE
  }
  
  if (class(test) == "data.frame") {
    orfs <- read.table(as.character(blast))
  } else {
    print("no orf to plot...")
    orfs <- suppressWarnings(as.data.frame(0))
    names(orfs) <- "V1"
    # Note: V2 is intentionally missing here to trigger the safe 'has_data' check
  }
  
  ###########################
  ## dot-plot (bottom left)##
  ###########################
  
  plot(
    x = 1, type = "n",
    xlim = c(0, bl$V3[1]), ylim = c(0, bl$V3[1]),
    col = "white", main = "",
    ylab = "TE consensus self dotplot (bp)",
    xlab = "TE consensus self dotplot (bp)",
    cex.lab = 2, cex.axis = 1.5, xaxs = "i", yaxs = "i"
  )
  
  for (i in 1:length(bl$V1)) {
    if (bl$V5[i] > bl$V4[i]) {
      segments(x0 = bl$V2[i], x1 = bl$V3[i], y0 = bl$V4[i], y1 = bl$V5[i], col = "black", lwd = 1.5)
    } else {
      segments(x0 = bl$V2[i], x1 = bl$V3[i], y0 = bl$V4[i], y1 = bl$V5[i], col = "#009E73", lwd = 1.5)
    }
  } 
  
  #####################################
  ## Annotation graph (bottom right) ##
  #####################################
  
  plot_te_structure <- function(orfs, bl) {
    # FIX: Ensure V2 exists to prevent "argument is of length zero" error
    has_data <- is.data.frame(orfs) && nrow(orfs) > 0 && !is.null(orfs$V2)
    
    # Constants for visual consistency
    row_height <- 1
    bar_thickness <- 0.25 
    
    # Logic to prevent "fat" boxes: force a minimum of 15 rows in plot space
    n_hits <- if (has_data) nrow(orfs) else 0
    view_limit <- max(12, n_hits) 
    
    plot(
      x = 1, type = "n",
      xlim = c(0, bl$V3[1]),
      ylim = c(0, view_limit + 1),
      col = "white", yaxt = "n", main = "",
      xlab = "TE consensus structure and protein hits (bp)",
      ylab = "", cex.lab = 1.5, cex.axis = 1.2, xaxs = "i", yaxs = "i"
    )
    
    if (!has_data) {
      text("No TE domain detected", x = bl$V3[1] / 2, y = view_limit / 2, cex = 2)
    } else {
      for (i in 1:nrow(orfs)) {
        # Scaling: place hits from the top down
        y_pos <- (view_limit - i + 1) * row_height
        
        orf_color <- if (orfs$V1[i] < orfs$V2[i]) "black" else "red"
        
        # Outer Frame
        rect(
          xleft = orfs$V1[i], xright = orfs$V2[i],
          ybottom = y_pos - (bar_thickness / 2),
          ytop = y_pos + (bar_thickness / 2),
          lwd = 1.2, border = orf_color
        )
        
        # Protein Hit Fill
        rect(
          xleft = orfs$V5[i], xright = orfs$V6[i],
          ybottom = y_pos - (bar_thickness / 2.5),
          ytop = y_pos + (bar_thickness / 2.5),
          col = paste0("#", orfs$V9[i]), border = NA
        )
        
        # Label
        text(
          paste0("(", orfs$V8[i], ") ", orfs$V3[i], " ", orfs$V4[i]),
          x = (orfs$V1[i] + orfs$V2[i]) / 2,
          y = y_pos + (bar_thickness / 1.5),
          pos = 3, cex = 1.1
        )
      }
    }
  }
  
  plot_te_structure(orfs, bl)
  
  if (tables == "TRUE") {
    names(bl) <- c("TE", "from.1", "to.1", "from.2", "to.2")
    write.table(bl, file = paste0(output, "/TE.self-blast.txt"), quote = FALSE, row.names = FALSE)
  }
}