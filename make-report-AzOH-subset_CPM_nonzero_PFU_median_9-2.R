# R script to process campaign and LiGA data - GL 2024-11-21


# MODULE 1 - SETUP and INITIALIZATION

library(dplyr)
library(tidyverse)
library(knitr)
library(rmarkdown)
library(edgeR)
library(scales)
library(grid)
library(gridExtra)
library(gtable)
library(readr)

campaign.filename.pattern <- "*UEA-campaign.csv"
identified.seq.table.name <- "./identified-sequence-table-Jan.txt"

LiGA.data.dir   <- "LiGA-data"
dict_dir        <- "Dictionaries"
order_table_dir <- "Order-tables"
campaign_dir    <- "Campaigns"
output.dir      <- "Outputs"
output.filename <- "Outputs/output-plot.pdf"
font.family     <- "Helvetica"   # other likely options "Arial", "Helvetica", "Arimo"
                    # must be installed on system
                    # mainfont option in header of report.Rmd should match 

# names of the internal standards
is.SDBs <- c("SDB267", "SDB268", "SDB271", "SDB272", "SDB273",
			 "SDB274", "SDB275", "SDB276", "SDB277", "SDB278", "SDB279", "SDB280")



if (!dir.exists(output.dir)){
    dir.create(output.dir, recursive = TRUE)
}
pdf(file = output.filename, width = 7.1, height = 7.1/1.618 * 1, onefile = TRUE)


# MODULE 2 - UTILITY FUNCTIONS

# extract the description fields from the LiGA datafile names
# inputs: a dataframe and the index of a column of the filenames
# results are additional columns prepended to the result frame

parse_LiGA_filenames <- function(files, index) {
    date_regex <- "(2[0-9]{3})(0[1-9]|1[012])(0[1-9]|[12][0-9]|3[01])"
    chem_regex <- "-([0-9]+)([A-Z]{2,})([a-z]{2,})([A-Z]{2,})-([A-Z]+)"
    parts <- str_match(basename(as.character(files[[index]])), 
                       paste0(date_regex, chem_regex))
    parsed_names <- data.frame(parts[,2:9], files, stringsAsFactors = FALSE)
    colnames(parsed_names)[1:8] <- c("Year", "Month", "Day", "Library",
                                    "Modification", "Target", "Postprocess", "ID")
    parsed_names
}

# given the root (or more) of the data filename and a directory:
# find the file's full name, then read it into a data.frame
# handles both space- and tab-separated files and gzip compression
read.LiGA.data <- function(filename, dir = "."){
    # search for the given filename part and obtain the full filename
    # obtain the full list of files and filter with grepl - avoids a regex
    
    filename.list <- list.files(path = dir)
    matches.index <- which(grepl(pattern = trimws(filename),
                                 filename.list, fixed = TRUE))
    if(length(matches.index) == 0){
        stop(paste0("No file matching '", filename, "' found in '", dir, "'."))
    } else if(length(matches.index) > 1){
        warning(paste0("Found multiple files, using '",
                       filename.list[matches.index[1]], "'."))
    }
    full.name <- file.path(dir, filename.list[matches.index[1]])
    print(paste0('Reading data from file: ', full.name))

    # handle case of a zip file - use first matching filename in ZIP directory
    if(grepl("\\.zip$", full.name)){
        zip.directory <- unzip(full.name, list = TRUE)
        matches.index <- which(grepl(pattern = trimws(filename),
                                     zip.directory$Name, fixed = TRUE))
        zip.name <- zip.directory$Name[matches.index[1]]
        zip.length = zip.directory$Length[matches.index[1]]
        # decompress and copy file from ZIP archive to temporary file
        temp.file.name <- tempfile()
        source.connection <- unz(full.name, zip.name, "rb")
        writeBin(readBin(source.connection, "raw", zip.length, 1),
                 temp.file.name, 1, useBytes = TRUE)
        close(source.connection)
        # switch the data source to the temporary file
        full.name <- temp.file.name
    }

    # read in start of the file as a list of lines and find header line
    file.header.lines <- readLines(full.name, n = 100)
    header.line <- which(grepl(" Mod Nuc AA ", file.header.lines))[1]

    # read in a data file - try reading a tab separated table first
    file.data <- tryCatch(
      read.table(full.name, header = TRUE, skip = header.line - 1,
                 sep = "\t", quote = "\"", stringsAsFactors = FALSE),
      # silence errors - may not be a tab-separated table
      error = function(condition) {NULL}
    )
    
    # if there is only one column, retry as space separated instead
    if(is.null(file.data) || NCOL(file.data) == 1){
      file.data <- tryCatch(
        read.table(full.name, header = TRUE, skip = header.line - 1,
                   sep = " ", quote = "\"", stringsAsFactors = FALSE),
        # silence errors - may not be a space-separated table
        error = function(condition) {NULL}
      )
    }

    # read.table has failed - try manually parsing lines instead
    # works for space separated tables when the Mod column (unquoted) contains spaces
    if(is.null(file.data)){
      # read in the file as a list of lines
      file.lines <- readLines(full.name)
      
      # find the indices of the header line and the Mod column within the header
      header.line <- which(grepl("Nuc", file.lines))[1]
      mod.column <- which(grepl("Mod", strsplit(file.lines[header.line], " ")[[1]]))
      header <- strsplit(file.lines[header.line], " ")[[1]]
      header.columns <- length(header)
      post.mod.columns <- header.columns - mod.column
      
   print(header.line)
   print(full.name)

      # split each line using spaces
      # keep the mod - 1 columns at the left and the post.mod.columns at the right
      # merge (paste) the mod column(s) back into one column using spaces,
      resplit.lines <- lapply(strsplit(file.lines[(header.line + 2):length(file.lines)], " "),
                              function(x){len <- length(x);
                              c(x[1:(mod.column - 1)],
                                paste(x[mod.column:(len-post.mod.columns)], collapse = " "),
                                x[(len-post.mod.columns+1):len])})
      
      file.data <- data.frame(row.names = 1:length(resplit.lines))
      for(i in 1:header.columns){
        data.column <- unlist(lapply(resplit.lines, function(x){x[i]}))
        conversion <- suppressWarnings(as.numeric(data.column))
        if(!any(is.na(conversion))){  # if column is all numeric convert it
          data.column <- conversion
        }
        file.data$data <- data.column
        colnames(file.data)[i] <- header[i]
      }
    }
    
    # as a sanity check, ensure that there is a Nuc (nucleotide sequence) column
    if(!"Nuc" %in% colnames(file.data)){
       warning(paste0("Expected column Nuc not found.  File '", full.name,
                      "' may not be read correctly"))
    }

    # if a temporary file has been used - delete it
    if(exists("temp.file.name")){
      unlink(temp.file.name)
    }
    
    # return the data.frame read by read.table()
    file.data
}


### read in the dictionary file, keep (non-blank) SDB and Axis.name pairs
### remove leading and trailing periods in the column names, duplicate SDB entries
read_dictionary_tables <- function(dictionary_basename, order_table_basename,
        label_column = "Axis.name", known_SDB_list, dict_dir = "dictionaries",
        order_table_dir = "figure data tables/axes") {

    # list possible dictionaries - if more than one is found, use the first
    dict_filenames <- list.files(path = dict_dir,
                        pattern = paste0(dictionary_basename, ".txt"))
    if(length(dict_filenames) < 1){
        stop("Dictionary not found: ", dictionary_basename)
    }

    # read in the choosen file
    dictionary <- read.delim(file.path(dict_dir, dict_filenames[1]),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # clean-up by stripping leading and trailing periods from column names
    colnames(dictionary) <- sub("^\\.*(.*?)\\.*$", "\\1", colnames(dictionary))

    # for simplicity drop the columns we are not using
    dictionary <- dictionary[,colnames(dictionary) %in% c("SDB", label_column, "Alphanum")]

    # standardize the name of the label column
    colnames(dictionary)[colnames(dictionary) == label_column] <- "Axis.name"

    # handle duplicated SDB entries - keep only one copy and label it "mix"
    dupl_dictionary_SDBs <- duplicated(dictionary$SDB)
    dictionary <- dictionary[!dupl_dictionary_SDBs,]
    dictionary$Axis.name[dictionary$SDB %in% dupl_dictionary_SDBs] = "mix"

    # drop unknown SDBs and blank entries from the dictionary
    dictionary <- dictionary[dictionary$SDB %in% known_SDB_list,]
    dictionary <- dictionary[dictionary$Axis.name != "" & dictionary$SDB != "",]

    # check that there is still at least one entry and three columns
    if(length(dictionary) != 3 | length(dictionary$SDB) < 1){
       stop(paste0("The dictionary file appears to be empty or lacks the SDB, Alphanum, and/or Label columns: ",
           dict_filenames[1], " ", input_columns$Label[i,]))
    }

    # read in the order table - drop extra columns and duplicate Alphanum entries
    order_table <- read.delim(file.path(order_table_dir, paste0(order_table_basename, ".csv")), 
                              header = TRUE, sep = ",", stringsAsFactors = FALSE)
    colnames(order_table) <- sub("^\\.*(.*?)\\.*$", "\\1", colnames(order_table))
    order_table <- order_table[,colnames(order_table) %in% c("Alphanum", "Order")]

    # remove all copies of duplicate Alphanum entries - multiples propagate through joins
    # TODO: could rescue repeated Alphanum's where all entries agree on the same order
    dupl_order_table_Alphanum <- order_table$Alphanum[duplicated(order_table$Alphanum)]
    order_table <- order_table[!order_table$Alphanum %in% dupl_order_table_Alphanum,]

    if(length(order_table) != 2 | length(order_table$Order) < 1){
        stop(paste0("The Order Table file appears to be empty or lacks the Order column: ",
                    order_table_basename))
    }

    # join the order table to the dictionary and reduce it to
    # the SDB, Axis Label, and Order columns
    dictionary <- dictionary %>% left_join(order_table, by = "Alphanum")
    dictionary <- dictionary[,colnames(dictionary) %in% c("SDB", "Axis.name", "Order")]
    
    # make NA in the Order a special position at the end of the list
    # and any duplicated SDBs a special order position at the front of the list
    dictionary$Order[is.na(dictionary$Order)] <- max(dictionary$Order, na.rm = TRUE) + 1000
    dictionary$Order[dictionary$SDB %in% dupl_dictionary_SDBs] <- min(dictionary$Order) - 1

    # sort the dictionary based on the Order column
    dictionary <- dictionary[order(dictionary$Order, nchar(dictionary$SDB), dictionary$SDB),]

    # return the dictionary table

    dictionary
}

# MODULE 3 - DATA PROCESSING

### read dataframe listing datafiles, columns to select, and matching dictionary
### return a single merged table of the requested columns
read_and_merge_columns <- function(input_columns, LiGA.data.dir = "LiGA-data",
        dict_dir = "dictionaries", order_table_dir = "order-tables/axes") {

    # verify the inputs
    if(length(colnames(input_columns)[colnames(input_columns) %in%
               c("Filename", "Columns", "Dictionary", "Labels", "OrderTable")]) != 5){
        stop("One or more required columns (Filename, Columns, Dictionary, Labels, OrderTable) not found.")
    }
    # TODO check that columns are numeric values

    if(length(unique(input_columns$OrderTable)) != 1){
        stop("Multiple OrderTables in input.")
    }

    # TODO: what if filenames do not meet the expected pattern?
    parsed_names <- parse_LiGA_filenames(input_columns, "Filename")
    parsed_names$root <- unite(unite(unite(parsed_names, Date, Year:Day, sep = ""),
                       root, Library:Postprocess, sep = ""), 
                       filename_root, Date:ID, sep = "-")$filename_root

    # if option columns are not present create minimal placeholders
    if(!"ExcludedSDBs" %in% colnames(parsed_names)){
       parsed_names$ExcludedSDBs <- ""
    }
    if(!"Type" %in% colnames(parsed_names)){
       parsed_names$Type <- ""
    }
    if(!"Reference" %in% colnames(parsed_names)){
        parsed_names$Reference <- ifelse(parsed_names$Type == "Test", "Control", "")
    }
    if(!"Standards" %in% colnames(parsed_names)){
       parsed_names$Standards <- ""			# TODO generate an approopriate placeholder
    }


    # read in the table of pre-labelled SDB sequences
    identified_seqs <- read.table(identified.seq.table.name,
                                  header = TRUE, stringsAsFactors = FALSE)


    # create data structures to store the results
    merged_data <- data.frame(SDB = character(), Barcode = character(),
        stringsAsFactors = FALSE)
    column_types <- c("Barcode", "Barcode")
    type.list <- unique(as.character(parsed_names$Type))
    data.source <- data.frame(Filename = c("", ""), Column = c(NA, NA), stringsAsFactors = FALSE)
    data.types <- setNames(data.frame(matrix(FALSE, nrow = 2, ncol = length(type.list)), stringsAsFactors = FALSE), type.list)

    # loop over the listed files, merging columns into the output data frame
    for(i in 1:length(parsed_names$Filename)){
        # read in the file - LiGA sequence count table
        file_data <- read.LiGA.data(parsed_names$Filename[i], dir = LiGA.data.dir) 
        # THIS VARIABLE STORES ALL THE DATA FROM THE TEXT FILE, NOT ASSGINED YET TO A KNOWN SDB

        # discard sequences that are not 96 base in length
        file_data <- file_data[nchar(file_data$Nuc) == 96,]
        # discard the sum total line (Nuc == XX) - can be recalculated later
        file_data <- file_data[file_data$Nuc != "XX",]

        # match each sequence with known SDBs - keep the SDB name and counts
        file_data <- full_join(file_data, identified_seqs, by = "Nuc")
        noncount.columns <- colnames(file_data) %in%
                c("Distance", "index", "mindex", "Primer", "Nuc", "Mod", "AA")
        file_data <- file_data[,!noncount.columns] # exclude the columns above, not needed anymore
        
        # handle any missing values - zero out in count columns
        # replace any NAs in the count data with zeros - restore NA SDBs
        sdb_names <- file_data$SDB
        file_data[is.na(file_data)] <- 0
        file_data$SDB <- sdb_names
        file_data <- file_data %>% group_by(SDB) %>% summarize_all(sum)

        # THIS VARIABLE STORES ALL THE DATA FROM THE TEXT FILE and ASSIGN THEM TO KNOWN SDBS, NOT DICT YET
        # UNKNNOWN SDBS are automatically named "NA" in the SDB COLUMN because there was no match
        # CONTAIS THE SDB NAME AND THE NUMBER OF READS FOR EACH REPLICATE
        
        
        
        # read in the dictionary file and order table
        dictionary <- read_dictionary_tables(parsed_names$Dictionary[i], parsed_names$OrderTable[i],
                          parsed_names$Label[i], known_SDB_list = identified_seqs$SDB,
                          dict_dir = dict_dir, order_table_dir = order_table_dir)
        # THIS IS JUST TO LOAD THE DICIONARY YOU WANT TO LOOK AT IT
        # CONTAINS THE SDB NUMBER, THE GLYCAN NAME AND THE ORDER
      

        
        
        # use Axis.name labels where the SDBs match, fallback to the SDB if unknown
        # and drop data for any excluded cases
        file_data <- full_join(file_data, dictionary, by = "SDB")
        file_data$Barcode <- file_data$Axis.name
        unknown_barcode <- is.na(file_data$Axis.name)
        file_data$Barcode[unknown_barcode] <- "" #file_data$SDB[unknown_barcode]
        file_data <- file_data[!file_data$SDB %in% parsed_names$ExcludedSDBs[[i]],]
        # THIS FINALLY LOADS FILE_DATA WITH ALL SDBS, UNKNOWN SDBS and NA
        # CONTAINS SDB NAME, READ COUNT, AXIS NAME AND ORDER

        #file_data <- file_data[order(file_data$Order, nchar(file_data$SDB), file_data$SDB, na.last = TRUE),]
        
        # list the indices of columns containing the read counts 
        # MAYBE IS HERE WHERE I NEED TO ADD THE OTHER COLUMN
        # exclude columns named "Mod", "Nuc", "AA", etc.
        data_columns <- (1:NCOL(file_data))[!colnames(file_data) %in% c("SDB",
            "Distance", "index", "mindex", "Primer", "Nuc", "Mod", "AA",
            "Axis.name", "Order", "Barcode")]
        
   
        # choose the requested subset of the data columns
        active_columns <- data_columns[as.vector(parsed_names$Columns[[i]])]
        for(col in 1:length(active_columns)){
            old.row <- data.source$Filename == parsed_names$Filename[i] &
                       data.source$Column == active_columns[col]
            if(!any(old.row)){
                # if not already read, append the requested column/row to the count data
                data.source <- rbind(data.source, data.frame(
                        Filename = parsed_names$Filename[i], Column = active_columns[col]))
                data.types <- rbind(data.types, data.types[1,])
            }
            # now that the column is added to data.types mark this entry
             data.types[data.source$Filename == parsed_names$Filename[i] &
                       data.source$Column == active_columns[col], as.character(parsed_names$Type[i])] <- TRUE
            # read in the column if new
            if(!any(old.row)){
                column_types <- c(column_types,
                                  as.character(parsed_names$Type[i]))
                if(active_columns[col] > ncol(file_data)){
                   stop("Requested column not in data file. Column:",
                        active_columns[col])
                }
                # new column name: filename root followed by the column number
                # column number as in the campaign file
                column_name <- paste0(parsed_names$root[i], 
                                      parsed_names$Columns[[i]][col])
                colnames(file_data)[active_columns[col]] <- column_name

                # extract a data_frame with the barcode and one column
                # of read counts, summing read counts over barcodes
                data_column <- file_data %>% group_by(SDB, Barcode) %>%
                           summarize(count_sum = sum(!!as.name(column_name)),
                                     .groups = 'drop')
                # note: added experimental .groups option to avoid warning

                # fix the column name and merge it to the result frame
                colnames(data_column)[3] <- column_name
                merged_data <- full_join(merged_data, data_column,
                                         by = c("SDB", "Barcode"))

                # replace any NAs in the count data with zeros - restore NA SDBs
                sdb_names <- merged_data$SDB
                merged_data[is.na(merged_data)] <- 0
                merged_data$SDB <- sdb_names
                
            }
        }
    }
    
    
    ordering_frame <- merged_data %>% left_join(dictionary,
                                        by = c("SDB", "Barcode" = "Axis.name"))
    merged_data$Order <- ordering_frame$Order
    data.types <- rbind(data.types, data.types[1,])
    ordering <- order(ordering_frame$Order,		# primary sort by provided Order
                      ordering_frame$Barcode == "",		# empty labels at the end
                      sub("-\\[[0-9<]*]$", "",  ordering_frame$Barcode),   # break ties by glycan (alphabetic)
                      suppressWarnings(as.integer(
                          sub(".*-\\[<?0*([1-9][0-9]*)\\]$", "\\1", ordering_frame$Barcode))),
                      ordering_frame$Barcode,		# finally just use the raw string
                      nchar(ordering_frame$SDB),		# then use the SDB
                      ordering_frame$SDB, na.last = TRUE)


#       colnames(order_table) <- sub("^\\.*(.*?)\\.*$", "\\1", colnames(order_table))
#ddd_frame <<- ordering_frame
    merged_data <- merged_data[ordering,]
    # CONTAINS DATA FROM TEST AND REFERENCE WITH SDB, READ COUNTS AND ORDER COLUMN
    
#    print(merged_data)
#    print("merged")
    
    
    # print(ordering_frame[ordering, c("Order", "SDB", "Barcode")])

    # finally extract a table of requested two-way comparisons
    contrasts <- unique(parsed_names[parsed_names$Reference != "", c("Type","Reference")])
    
    list(table = merged_data, Contrasts = contrasts, Column.types = data.types, dict = dictionary)
}

make.barplot <- function(output_df, labels.on, basename, font.family = font.family) {
    relative.error <- sqrt(output_df$Test_SD ** 2    / output_df$Test_CPM ** 2 +
                           output_df$Control_SD ** 2 / output_df$Control_CPM ** 2)
    plot_data <- data.frame(
                Label = reorder(factor(paste0(output_df$Barcode, " ", output_df$SDB)), output_df$rank),
                Order = output_df$Order,
                Fold.change = 2 ** output_df$logFC,
                QValue = p.adjust(output_df$PValue, "BH"),
                Noise = (2 ** output_df$logFC) * (1 + relative.error),
                Valid = output_df$Valid
                )
    #plot_data <- plot_data[!is.na(plot_data$QValue),]
    error.bar.data <- plot_data[!is.na(plot_data$Noise),]  # avoid plotting NAs
    significant.data.pos <- plot_data[plot_data$Fold.change > 1 & plot_data$QValue <= 0.05 & plot_data$Valid,]
    significant.data.neg <- plot_data[plot_data$Fold.change < 1 & plot_data$QValue <= 0.05 & plot_data$Valid,]

    bar_plot <- ggplot(data = plot_data,
                     aes(x = Label, y = Noise, fill = Valid)) +
            theme_light() + theme(text = element_text(family = font.family, size = 7),
                                  axis.title = element_text(family = font.family, size = 7),
                                  panel.grid.major.x = element_blank(),
                                  panel.grid.minor.y = element_blank(),
                                  legend.position = "top") +
            ggtitle(unique(paste0(filename_root, ": ", difference.frame$Type, " vs. ", difference.frame$Reference))) +
            labs(x = "Glycan", y = "Fold Change") +
            #scale_y_sqrt(limits = c(0, 1250), expand = expand_scale(mult = c(0, 0.0))) +
            #scale_y_sqrt(expand = expand_scale(mult = c(0, 0.08))) +
            #scale_y_continuous(trans = "identity",
            #                   # limits = NA, #c(0, 15.0),
            #                   expand = expand_scale(mult = c(0, 0.08))) +
            #scale_y_log10() +
            scale_y_sqrt(breaks = c(1, 20, 50, 100, 250, 500, 1500)) +
            #scale_y_continuous(limits = c(0, 100)) +
            geom_hline(yintercept = 0.01) + 
            geom_col(aes(y = Fold.change), 
                     #position = "dodge", 
                     position = position_dodge2(padding = 0.1, width = 0.7),
                     linewidth = 0.0, width = 0.5) +
            geom_errorbar(data = error.bar.data,	# top of error bars
                          #aes(ymin = Noise, ymax = Noise, width = 2)) +
                         aes(ymin = Noise, ymax = Noise, width = 0.8  #, color = "lightblue"
                             )) +
      
            geom_linerange(data = error.bar.data,	# vertical of error bar
                           aes(ymin = Fold.change, ymax = Noise       #, color = "lightblue"
                               ))
    
    if(nrow(significant.data.pos) > 0){
        bar_plot <- bar_plot + geom_text(label = "\u2022",  # "\u2022"= round dot;  "."= square;  "*"= *.
                                         hjust = 0.5,    # Adjust this value to move the dot left or right
                                         vjust = -0.05,  # Adjust this value to move the dot higher
                                         size = 3,       # Adjust this value to change the dot size
                                         data = significant.data.pos, aes(y = Noise))
    }
    #if(nrow(significant.data.neg) > 0){
    #    bar_plot <- bar_plot + geom_text(label = "*", hjust = 0.5, colour = "red",
    #                                     data = significant.data.neg, aes(y = Noise))
    #}

    bar_plot <- bar_plot + theme(axis.text.y  = element_text(face = "bold", size = 7)) +
            scale_fill_manual(values = c("#888888", "#000000"),
                     #labels = c("FALSE", "TRUE"),
                     #guide = guide_legend(reverse = TRUE)
                     )
            theme(axis.text.y  = element_text(size = 7),
                  axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3))

    if(labels.on){
        bar_plot <- bar_plot + theme(axis.text.x  = element_text(
                   face = "bold", angle = 90, hjust = 1, vjust = 0.5, size = 4))
    } else {
        bar_plot <- bar_plot + theme(axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank(),
                                     axis.title.x = element_blank())
    }

    bar_plot = bar_plot   #, output_df = output_df)
}

 


# routine to run the process on the 
scriptDir <- getwd()
input_filename_list <- list.files(path = campaign_dir, pattern = campaign.filename.pattern)



data_list <- list()
diff.tests <- list()
for(filenumber in 1:length(input_filename_list)){
    data_list <- list()
    diff.tests <- list()
    # drop the extension (if any) from the filename to get the root used for naming new outputs
    # note: special handling of indexing at length == 1 - does not work for all dot filenames
    filename_root <- unlist(strsplit(basename(input_filename_list[filenumber]), split = "[.]"))
    filename_root <- paste(filename_root[1:length(filename_root)-1], collapse = ".")

    # read the campaign file
    print(paste0("Processing campaign file: ", filename_root))
    testing <- read.table(file.path(campaign_dir, input_filename_list[filenumber]),
                          sep = ",", quote = "\"", header = TRUE, stringsAsFactors = FALSE)


    # check for the expected columns and skip any file without
    if(length(unique(colnames(testing)[colnames(testing) %in%
               c("Type", "Filename", "Columns", "Dictionary", "Labels", "OrderTable")])) != 6){
   		stop("One or more required columns (Type, Filename, Columns, Dictionary, Labels, OrderTable) not found.")
    }

    # create the directory for output if not already existing - parents too
    # remove spaces and underscores as per http://yihui.org/2018/03/space-pain
    compare_dir_name <- paste0(output.dir, "/", gsub("[ _]", "-", filename_root), "-compare-output")
    if (!dir.exists(compare_dir_name)){
        dir.create(compare_dir_name, recursive = TRUE)
    }

    testing$Type <- factor(testing$Type)
    testing$Columns <- lapply(strsplit(as.character(testing$Columns), ","), as.integer)

    # R converts blank strings into NA on input so reverse this before extracting the list
    testing$ExcludedSDBs[is.na(testing$ExcludedSDBs)] <- ""
    testing$ExcludedSDBs <- lapply(strsplit(testing$ExcludedSDBs, ","), trimws)

    testing_data <- read_and_merge_columns(testing, LiGA.data.dir = LiGA.data.dir,
            dict_dir = "dictionaries", order_table_dir = order_table_dir)
    # RETURNS A LIST OF ALL SDBS (NOT UNKNOWN), AXIS.name and ORDER
    # UNKNOWN ARE EXCLUDED
    #NO READ COUNTS HERE


    # all the read counts have now been read in - write them out in the count-table
	sum.table <- data.frame(Glycan = testing_data$table$Barcode,
                       SDB = testing_data$table$SDB,
                       testing_data$table[,3:(ncol(testing_data$table)-1)],
                       check.names = FALSE)

	# delete unkown SDBs
	###sum.table <- sum.table[sum.table$Glycan != "",]
	
	write.table(sum.table,
        file = paste0(compare_dir_name, "/", "count-table-", filename_root, ".tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE)

    
    # check that at least one comparison is requested
    if(NROW(testing_data$Contrasts) < 1) {
        print("No Test - Control contrasts available.")
    }

    is.indices <- which(testing_data$table$SDB %in% is.SDBs)

    # loop over three subsets of the data
    for(case in 3:3){							# or only do one case instead
        if(case == 1){
            analysis_case <- "all"
            active_rows <- TRUE  #!is.na(testing_data$sequence)
        } 
        if(case == 2){
            analysis_case <- "knownSDBs"
            active_rows <- !is.na(testing_data$table$SDB)
        } 
        if(case == 3){
            analysis_case <- "knownGlycans"
            active_rows <- testing_data$table$Barcode != ""
            # The code belows analyzes the presence of barcode with 0 reads
        } 
  
      
        # use edgeR, apply TMM and calculate the mean & st. dev. for each condition

        # generate a mask to select data rows/columns  
        count.mask <- apply(testing_data$Column.types, 1, any)

        # generate a placeholder group - select the first if multiple groups
        a.group <- apply(testing_data$Column.types, 1, function(x){head(names(which(x)), n=1)})

        dge <- DGEList(counts = testing_data$table[active_rows,count.mask],
                      group = factor(unlist(a.group[count.mask])))
        dge <- calcNormFactors(dge, method = "TMM")
        ### print("TMM over all cases")
        ### print(dge$samples)
        ### dge.normalized <- cpm(dge, normalized.lib.sizes = TRUE)

        # is.SDBs - names of the internal standards - set at top of file

        # normalize a "library" of just the internal standards
        is.rows <- (testing_data$table$SDB %in% is.SDBs) & active_rows
        is.dge <- DGEList(counts = testing_data$table[is.rows,count.mask],
                          group = factor(unlist(a.group[count.mask])))
        is.dge$samples$lib.size <- dge$samples$lib.size
        is.dge <- calcNormFactors(is.dge, method = "RLE")  # I think RLE is most appropriate uses all rows

        print("Per Sample Sums of Internal Standards:")
        print(colSums(testing_data$table[is.rows,count.mask]))
        print("Sample Normalizations:")
        print(is.dge$samples)
 
        # copy the normalization factors to the whole library and calculate CPM
        dge$samples$norm.factors <- is.dge$samples$norm.factors
        dge.normalized <- cpm(dge, normalized.lib.sizes = TRUE)
 

                # extract the subset of the test and control data
        testing_data$Contrasts$Type      <- as.character(testing_data$Contrasts$Type)
        testing_data$Contrasts$Reference <- as.character(testing_data$Contrasts$Reference)
        group.names <- data.frame(Group = unique(c(testing_data$Contrasts$Type,
                                                   testing_data$Contrasts$Reference)), stringsAsFactors = FALSE)
        group.names$Name <- make.names(unique(group.names$Group), unique = TRUE, allow_ = FALSE)
        contrast.table <- left_join(left_join(testing_data$Contrasts, group.names, by = c("Type" = "Group")),
                                                                      group.names, by = c("Reference" = "Group"))
        contrast.table$String <- paste0(contrast.table$Name.x, "-", contrast.table$Name.y)

        print(contrast.table)
        
        #TODO: Would it be better to estimate dispersion from ALL the data?
        for(pair in 1:nrow(contrast.table)){
            # extract masks selecting which datasets are used in the current contrast        
            type.data       <- testing_data$Column.types[,contrast.table$Type[pair]]
            compare.to.data <- testing_data$Column.types[,contrast.table$Reference[pair]]

            # need omit any datasets on both sides of the comparison    
            both.data <- type.data & compare.to.data
            if(any(both.data)){
                stop("Dataset(s) used on both sides of comparison:", colnames(testing_data$table)[both.data])
            }             

            # combine the indexing masks and extract the group names
            subset_indices <- (type.data | compare.to.data) & !both.data
            subset_group <- ifelse(testing_data$Column.types[subset_indices,contrast.table$Type[pair]],
                    contrast.table$Type[pair], contrast.table$Reference[pair])

            # edgeR object for the subset, take norm.factors from whole TMM
            dge_subset <- DGEList(counts = testing_data$table[active_rows,subset_indices],
                                 norm.factors = dge$samples$norm.factors[subset_indices[count.mask]],
                                 group = subset_group)

            # identify the low-count data (rows) in the subset
            Valid <- filterByExpr(dge_subset)
            #dge_subset <- dge_subset[Valid, , keep.lib.sizes = FALSE]

            # the following is from the tutorial / Jessica's script, needs to be further modified
            # build the model matrix and fix-up the output names (meet R's variable name rules)
            names <- left_join(data.frame(Group = as.character(dge_subset$samples$group),
                                 stringsAsFactors = FALSE), group.names, by = "Group")$Name
            design <- model.matrix(~ 0 + names, data = dge_subset$samples)
            colnames(design) <- gsub("^names", "", colnames(design))
            
            # estimate dispersion (by one method or another)
            tmp <- try(estimateDisp(dge_subset, design, robust = TRUE), TRUE)
            if(inherits(tmp, "try-error")){
                dge_subset <- estimateCommonDisp(dge_subset)
                dge_subset <- estimateTagwiseDisp(dge_subset)
            } else {
                # print("Using estimateDisp")
                dge_subset <- estimateDisp(dge_subset, design, robust = TRUE)
                ###plotMDS(dge_subset, top = dim(dge_subset$counts)[1], col = colgroup)
            }
            
            # now use GLM to model the glycan dispersion and test for differences
            fit <- glmQLFit(dge_subset, design, robust = TRUE)
            contrasts <- makeContrasts(contrasts = contrast.table$String[pair], levels = design)
			qlfi <- glmQLFTest(fit, contrast = contrasts)
            
            # build a data.frame of output information...
            difference.frame <- data.frame(
                    rank = 1:sum(active_rows),
                    type.rank = 1 + length(diff.tests),
                    testing_data$table[active_rows, c("SDB", "Barcode", "Order")],
                    contrast.table[pair,1:2])

            # ...include the TMM normalized mean and st. deviation...
            # to calculate relative error for the error bars
            difference.frame[["Test_CPM"]]    <- apply(dge.normalized[,type.data[count.mask]], 1, mean)
            difference.frame[["Test_SD"]]     <- apply(dge.normalized[,type.data[count.mask]], 1, sd)
            difference.frame[["Control_CPM"]] <- apply(dge.normalized[,compare.to.data[count.mask]], 1, mean)
            difference.frame[["Control_SD"]]  <- apply(dge.normalized[,compare.to.data[count.mask]], 1, sd)

            # join the quasi-liklihood F-test results back the appropriate IDs
            difference.frame <- left_join(difference.frame,
                data.frame(testing_data$table[active_rows, c("SDB", "Barcode", "Order")],
                    Test.CPM = exp(fit$coefficients[,contrast.table$Name.x[pair]])*1e6,
                    Control.CPM = exp(fit$coefficients[,contrast.table$Name.y[pair]])*1e6,
                    qlfi$table, QValue = p.adjust(qlfi$table$PValue, "BH"), Valid),
                by = c("SDB", "Barcode", "Order"))

            # NAs are inserted for low-count cases excluded from the testing
            #difference.frame$PValue[!Valid] <- NA

            bar.plot <- make.barplot(difference.frame, TRUE, filename_root,
                                     font.family = font.family)
            diff.tests[[length(diff.tests)+1]] <- list(frame = difference.frame,
                      plot = bar.plot,
                      data1.set = type.data,
                      data2.set = compare.to.data)

         }
    }   ### end of subset case loop - still in campaign file loop


	# arrange the difference tests into one large table for plotting
	long.data <- do.call(rbind, lapply(diff.tests, function(x){x$frame}))
	#long.data$PValue <- ifelse(long.data$Valid, long.data$PValue, NA)
	long.data$QValue <- p.adjust(long.data$PValue, method = "BH")

	long.data$Glycan <- sub("-\\[[0-9<]*]$", "",  long.data$Barcode)
	long.data$Density <- suppressWarnings(as.integer(
    	                      sub(".*-\\[<?0*([1-9][0-9]*)\\]$", "\\1", long.data$Barcode)))

	permutation <- order(long.data$Order,			# primary sort by provided Order
    	                 long.data$Barcode == "",           # empty labels at the end
        	             long.data$Glycan,		        # group each glcan's entries (alphabetic)
            	         long.data$Density,
                	     #long.data$Barcode,	                # finally just use the raw string
                    	 nchar(long.data$SDB),		# resolve ties by the SDB, but put
                     	long.data$SDB, na.last = TRUE)	# shorter ones first (numeric order)

	long.data <- long.data[permutation,]
	x.position.frame <- data.frame(unique(long.data[,c("Order", "Barcode", "SDB")]), stringsAsFactors = FALSE)
	x.position.frame$xpos = 1:nrow(x.position.frame)
	long.data <- left_join(long.data, x.position.frame, by = c("SDB", "Barcode", "Order"))

	long.data <- left_join(long.data, data.frame(type.rank = c(1:7), ypos=c(7, 4, 3, 2, 5, 6, 1)), by = "type.rank")


	long.data$Label <- paste0(sub("'", "\u2032", long.data$Barcode), " ", long.data$SDB)
	long.data$Label <- reorder(factor(long.data$Label), long.data$xpos)

	significant.rows.1 <- long.data$Valid & long.data$QValue <= 0.05 & log2(2) <= long.data$logFC & long.data$logFC < log2(3.8)
	significant.rows.2 <- long.data$Valid & long.data$QValue <= 0.05 & log2(3.8) <= long.data$logFC





	#comparison-table
	compare.table <- long.data[order(long.data$PValue), c("Barcode", "SDB", "Type", "Reference",
    	         "Test.CPM", "Test_SD", "Control.CPM", "Control_SD", "logFC", "F", "PValue", "QValue")]
	colnames(compare.table)[3:12] <- c("Test", "Reference", "Test (CPM)", "Test (SD)",
    	                               "Reference (CPM)", "Ref (SD)", "log2FC", "F", "p-value", "q-value")

	# Merge columns from sum.table based on matching 'SDB' values
	compare.table <- cbind(compare.table, sum.table[match(compare.table$SDB, sum.table$SDB), 3:NCOL(sum.table)])
	# Sort rows alphabetically by 'SDB' (or another column if desired)
	compare.table <- compare.table[order(compare.table$SDB), ]




	# experimental test - using the row numbers to control the order
	Order <- as.integer(row.names(compare.table))

	# Add the "Order" column and the "NonZeroReference" to compare.table

	reference_column <- grep("Reference \\(CPM\\)", colnames(compare.table))  # Adjust "Control" to match your actual column names

	compare.table$NonZero <- compare.table[, reference_column] >= 50
	compare.table$log2FC[!compare.table$NonZero] <- NA
	###compare.table$NonZero <- compare.table[, "Reference (CPM)"] >= 50
	compare.table$Order <- Order



output_path <- paste0(compare_dir_name, "/", "comparison-table-", sub(".[^.]*$", "", input_filename_list[filenumber]), ".tsv")
write.table(compare.table, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)


#data-table
row.idx <- match(testing_data$table$SDB, compare.table$SDB)
row.flags <- !is.na(row.idx)
ratmir.table <- cbind(testing_data$table[row.flags,-ncol(testing_data$table)],
                      compare.table[row.idx[row.flags], c("Test (CPM)", "Test (SD)", "Reference (CPM)", "Ref (SD)", "log2FC", "q-value")])

colnames(ratmir.table)[ncol(ratmir.table)-1] <- "FC"

ratmir.table$FC <- 2**ratmir.table$FC
write.table(ratmir.table,
        file = paste0(compare_dir_name, "/", "data-table-", sub(".[^.]*$", "", input_filename_list[filenumber]), ".tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE)


# Read data from the tsv. file generated above
data <- read_tsv(output_path)
print(data)





#input_test =10^6
cat("Enter PFU values separated by whitespace (e.g., 5.06E+02 6.48E+02 8.98E+02):\n")
input_string <- readline()
values <- as.numeric(strsplit(input_string, "\\s+")[[1]])
input_test <- prod(values)^(1 / length(values))
total_FC <- sum(2^data$log2FC, na.rm = TRUE)  # Summing values in the 'FC' column



data$PFU_column <- (2^data$log2FC / total_FC) * input_test

#filtered_data <- data %>%  #if i use this code I can generate the plot without spaces 
  #filter(Nonzero == TRUE)




# print(total_sum)


# Ensure unique order and create a data frame
plot_data <- data %>%
  arrange(Order) %>% # Arrange data by Order column
  mutate(xlabel = factor(Barcode, levels = unique(Barcode))) # Convert Barcode to a factor with the correct levels


# Display options to the user
cat("Select the glycan types you want to plot:\n")
cat("1 - GalNAc/Gal-containing glycans\n")
cat("2 - Blood Groups\n")
cat("3 - Mannose-containing glycans\n")
cat("4 - N-glycans\n")
cat("5 - Sia-containing glycans\n")
cat("0 - Plot all glycans\n")


selected_group <- as.integer(readline("Enter your choice: "))


glycan_groups <- list(
  "1" = c("a-GalNAc-COOH-[135]", "a-GalNAc-COOH-[540]", "a-GalNAc-COOH-[1350]",
          "b-GalNAc-Az-[10]", "b-GalNAc-Az-[20]", "b-GalNAc-Az-[50]", "b-GalNAc-Az-[100]",
          "b-GalNAc-Az-[500]", "b-GalNAc-Az-[1000]", "b-GalNAc-COOH-[10]", "b-GalNAc-COOH-[20]",
          "b-GalNAc-COOH-[50]", "b-GalNAc-COOH-[135]", "b-GalNAc-COOH-[540]", "b-GalNAc-COOH-[1350]",
          "Tri-GalNAc-COOH-[100]", "Tri-GalNAc-COOH-[500]", "Tri-GalNAc-COOH-[1000]",
          "GN-[590]", "GN-[1050]", "Lac-peg4-[1080]", "Lac-[1030]", "Lac-[1240]",
          "LacNAc, LN-[970]", "LacNAc, LN-[1110]", "Lec-[680]", "Galili-tri-[1000]",
          "Pk-[860]", "Gala3-type1-[350]", "B2 tri-[350]", "P1 tri-[620]", "LacDiNAc-[50]",
          "LNT-2-[430]", "GNLN-[810]", "3'GN type1-[860]", "LNnT-[240]", "Globoside-P-[730]",
          "Globoside-P-[1030]", "P1 tetra-[970]", "P1 penta-[620]", "P1x penta-[620]",
          "Tri-LN-[380]", "Di-N3-[970]", "AzOH", "blank"),
  "2" = c("H type 1-[700]", "H-type-1-COOH-[135]", "H-type-1-COOH-[1350]", "H type 2-[540]",
          "H-type-2-COOH-[135]", "H-type-2-COOH-[540]", "H-type-2-COOH-[1350]", "H-type-3-COOH-[100]",
          "H-type-3-COOH-[500]", "H-type-3-COOH-[1000]", "H-type-4-COOH-[100]", "H-type-4-COOH-[500]",
          "H-type-4-COOH-[1000]", "H-type-5-COOH-[100]", "H-type-5-COOH-[500]", "H-type-5-COOH-[1000]",
          "H type 6-[950]", "H-type-6-COOH-[100]", "H-type-6-COOH-[500]", "H-type-6-COOH-[1000]",
          "A type 1-[700]", "A-type-1-COOH-[135]", "A-type-1-COOH-[540]", "A-type-1-COOH-[1350]",
          "A type 2-[920]", "A-type-2-COOH-[135]", "A-type-2-COOH-[540]", "A-type-2-COOH-[1350]",
          "A-type-3-COOH[135]", "A-type-3-COOH[540]", "A-type-3-COOH[1350]", "A-type-4-COOH-[100]",
          "A-type-4-COOH-[500]", "A-type-4-COOH-[1000]", "A-type-5-COOH-[135]", "A-type-5-COOH-[540]",
          "A-type-5-COOH-[1350]", "A type 6-[590]", "A-type-6-COOH-[135]", "A-type-6-COOH-[540]",
          "A-type-6-COOH-[1350]", "B type 1-[620]", "B-type-1-COOH-[135]", "B-type-1-COOH-[540]",
          "B-type-1-COOH-[1350]", "B type 2-[970]", "B-type-2-COOH-[540]", "B-type-3-COOH-[100]",
          "B-type-3-COOH-[500]", "B-type-3-COOH-[1000]", "B-type-4-COOH-[135]", "B-type-4-COOH-[540]",
          "B-type-4-COOH-[1350]", "B-type-5-COOH-[540]", "B-type-5-COOH-[1350]", "B type 6-[920]",
          "B-type-6-COOH-[135]", "B-type-6-COOH-[540]", "B-type-6-COOH-[1350]", "2'F-B type 2-[220]",
          "2'F-B type 2-[760]", "H2-[430]", "H3-[190]", "Tri-AN3-[1080]", "LeA-[950]", "Lex-[810]",
          "Gala3Lex-[620]", "Di-Lex-[410]", "Ley-Lex-[350]", "LeALex-[350]", "Lex-LeA-[410]",
          "Lec-LeX-[570]", "Tri-Lex-[430]", "Ley-Di-Lex-[510]","AzOH","blank"),
  "3" = c("aMan-[840]", "Mana-sp-[20]", "Mana-sp-[50]", "Mana-sp-[150]", "Mana-sp-[500]", "Mana-sp-[700]",
          "Mana-C4-COOH-[25]", "Mana-C4-COOH-[50]", "Mana-C4-COOH-[150]", "Mana-C4-COOH-[500]",
          "Mana-PEG4-COOH-[25]", "Mana-PEG4-COOH-[50]", "Mana-PEG4-COOH-[150]", "Mana-PEG4-COOH-[500]",
          "L-Mana-C8-COOH-[25]", "L-Mana-C8-COOH-[50]", "L-Mana-C8-COOH-[150]", "L-Mana-C8-COOH-[500]",
          "Man4-C5-N3-[20]", "Man4-C5-N3-[50]", "Man4-C5-N3-[150]", "Man4-C5-N3-[500]", "Man4-C5-N3-[700]",
          "Man4-P6-N3-[25]", "Man4-P6-N3-[50]", "Man4-P6-N3-[150]", "Man4-P6-N3-[500]", "Man4-P9-N3-[25]",
          "Man4-P9-N3-[50]", "Man4-P9-N3-[150]", "Man4-P9-N3-[500]", "Man4-P12-N3-[25]", "Man4-P12-N3-[50]",
          "Man4-P12-N3-[150]", "Man4-P12-N3-[500]", "(Man)3-[410]", "(Man)3-[1300]", "(Man)3-[1730]","AzOH","blank"),
  "4" = c("11-[50]", "11-[150]", "11-[500]", "11-[750]", "11-[780]", "11-[1000]",
          "6-[50]", "6-[150]", "6-[500]", "6-[750]", "6-[810]", "6-[1000]",
          "10-[50]", "10-[150]", "10-[500]", "10-[540]", "10-[750]", "10-[1000]",
          "9-[50]", "9-[160]", "9-[510]", "9-[730]", "9-[950]", "9-[970]",
          "8-[50]", "8-[140]", "8-[490]", "8-[760]", "8-[1000]",
          "7-[50]", "7-[140]", "7-[510]", "7-[760]", "7-[970]","AzOH","blank"),
  "5" = c("GM1-[460]", "GM1-[1190]", "GM2-[190]", "GM2-[460]", "GM2-[590]", "GM2-[1190]",
          "CT|Sda-[350]", "GD1a-[110]", "3'S-Di-LeA-[160]", "3'SLeA-Lex-[570]", "3'S-Tri-LeX-[160]",
          "3'SLec-[430]", "GM3-[190]", "GM3-[540]", "GM3-[1190]", "3'SLDN-[460]", "3'S-Di-Lec-[270]",
          "3'SLecLN-[860]", "3'SLN-Lec-[410]", "3'S-Di-LN-[1350]", "3'STri-LN-[160]", "3'SLec (Gc)-[350]",
          "3'SL (Gc)-[490]", "3'SLN (Gc)-[570]", "3'-KDNLec-[760]", "3'KDNLN-[510]", "6'SL-[620]",
          "6'SL-[760]", "6'SLN-[350]", "6'SLN (Gc)-[430]", "6'SLDN-[620]", "6'SLN-Lec-[300]",
          "6'S-Di-LN-[300]", "GD2-[460]", "GD3-[320]", "GT2-[160]", "GT3-[510]", "GQ2-[140]", "(Galf)4-[<8]","AzOH", "blank")
)



plot_data_filtered <- plot_data %>%
  {
    if (selected_group == 0) {
      . # Plot all glycans
    } else if (selected_group %in% c(1, 2, 3, 4, 5)) {
      filter(., Barcode %in% glycan_groups[[as.character(selected_group)]])
    } else {
      stop("Invalid input. Please enter a correct number.")
    }
  } %>%
  mutate(
    bar_height = ifelse(NonZero, 2 ^ `log2FC`, NA) # Set to NA where NonZero is FALSE
  )

### try to print nice labels
plot_data_filtered$Decoration <- sub("-\\[[^]]*\\]$", "", plot_data_filtered$Barcode)
plot_data_filtered$Density <- sub(".*\\[([^]]*)\\]|(.*[^]]())$", "\\1", plot_data_filtered$Barcode)
plot_data_filtered[!plot_data_filtered$NonZero,]$log2FC <- NA

data_by_density <- plot_data_filtered %>% 
						group_by(xlabel) %>% 
						summarize(Density = first(Density),
								  log2FC = median(log2FC, na.rm = TRUE),
								 ) %>%
						mutate(bar_height = 2 ^ log2FC)

# estimate some likely breaks on the PFU axis
rough_breaks <- ((0:10 / 10) * sqrt((50 * 1.05) * input_test / total_FC)) ** 2
cleaner_breaks <- 10 ** unique(round(log10(rough_breaks)))
cleaner_breaks <- c(0, 100, 500, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000)

# Plot the bar chart with error bars
#g.bar <- ggplot(plot_data_filtered, aes(x = xlabel, y = bar_height)) + # changed here to remove bars with problematic sdbs

plot.max <- 50 * 1.05   # limit plot to 50 + 5% expansion
g.bar <- ggplot(data_by_density, aes(x = xlabel, y = pmin(bar_height, plot.max))) + # changed here to remove bars with problematic sdbs
  geom_rect(mapping = aes(xmin = min(as.numeric(xlabel)) - 1,
                          xmax = max(as.numeric(xlabel)) + 1,
                          ymin = bar_height[xlabel == "blank"],
                          ymax = bar_height[xlabel == "AzOH"]),
            fill = "firebrick2") +
  geom_bar(stat = 'identity', width = 0.7, position = position_dodge2(padding = 0.1, width = 1.5), fill = "black") +
  #scale_y_continuous(trans = "log10", limits = c(1, NA), breaks = c(1, 2, 5, 10, 20, 50, 100, 500, 5000, 20000)) +
  scale_x_discrete(expand = expansion(add = 0.0)) +
  #scale_x_discrete(drop = FALSE) +   # leave space for missing columns
           scale_y_sqrt(breaks = c(0.1, 1, 2, 5, 10, 20, 30, 40, 50, 100, 250, 500, 1000, 2000, 4000, 8000, 20000, 40000),
               minor_breaks = c(0.5,3,4,6,7,8,9,15,25,35,45),
               limits = c(NA, plot.max),  # Capping the fold change at 50
               name = "Fold Change",  # Primary y-axis label
               expand = expansion(mult = c(0, 0)),
               guide=guide_axis(minor.ticks = TRUE),
               sec.axis = sec_axis(
  									~ . * (input_test / total_FC),  # Transform FC to PFU proportionally
    								name = "PFU",
    								breaks = cleaner_breaks)) +  # Secondary y-axis label) +
  labs(title = '', x = 'Glycans', y = 'Fold Change') +
  theme_minimal() +
  theme(
    text = element_text(family = font.family, size = 7), # Set font to Helvetica and size to 7
    #plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.title.x = element_blank(),  #text(size = 7), # Set x-axis label font size
    axis.title.y = element_text(size = 7, colour = "black"), # Set y-axis label font size
    axis.text.x = element_blank(),  #text(size = 5, angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels 90 degrees
    #axis.ticks.length.x = unit(2, "mm"), # Customize x-axis tick length
    axis.ticks.x = element_blank(), #element_line(linewidth = 0.5), # Customize x-axis tick appearance
    axis.ticks.y = element_line(linewidth = 0.5), # Customize y-axis tick appearance
    axis.ticks.length.y = unit(2, "mm"), # Customize y-axis tick length
    axis.text.y = element_text(size = 7, color = "black"), # Set y-axis text font size
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 0, b = 0),  # remove top and bottom margins
    axis.line.y.right = element_blank(),# Remove minor grid lines
    legend.position = "none", # Ensure the legend is removed
    panel.border = element_rect(colour = "#888888", fill = NA, linewidth = 1) # Add border around the plot
  )

# Determine the range of Fold Change
   fold_change_values <- 2 ^ plot_data_filtered$`log2FC`
   min_fc <- min(fold_change_values)
   max_fc <- max(fold_change_values)
    
### try to print nice labels
plot_data_labels <- plot_data_filtered %>% 
						group_by(Decoration) %>%
						summarize(xmin = min(as.numeric(xlabel)), xmax = max(as.numeric(xlabel))) %>%
						arrange(xmin)
plot_data_labels$xpos <- (plot_data_labels$xmin + plot_data_labels$xmax) / 2						




# Customize this value to unify maximum value
   max_fc <- 200    
   HeatMap_threshold <- 200

# Adjust the heatmap with correct scaling
         plot_data_filtered <- plot_data_filtered %>% group_by(xlabel) %>% mutate(ylabel = row_number() - 0.5, ycount = n() )
         plot_data_filtered$NonZero <- compare.table$NonZero[match(plot_data_filtered$GlycanName, compare.table$GlycanName)]
         plot_data_filtered <- subset(plot_data_filtered, NonZero)
         g.heatmap <- ggplot(data_by_density, aes(x = xlabel, y = 0.0, height = 1.0, fill = sqrt(2 ^ `log2FC`))) +
           #geom_tile(mapping = aes(x = xlabel, y = ylabel/ycount, height = 1/ycount, colour = "skyblue"),
           #          data = plot_data_filtered[!plot_data_filtered$NonZero,]) +
           #geom_tile() + 
           #geom_tile(mapping = aes(x = xlabel, y = -0.5, height = 0.5, fill = sqrt(2 ^ log2FC)), data = data_by_density) +
           geom_tile(mapping = aes(
             x = xlabel,
             y = -0.5,
             height = 0.5,
             fill = ifelse(NonZero, sqrt(2 ^ log2FC), NA)
           ), data = data_by_density) +
           geom_text(mapping = aes(x = xlabel, y = -1, label = Density, angle = 90, fontface = "plain", hjust = 1, size = 5.2), 
           			data = data_by_density, size.unit = "pt", show.legend = FALSE, inherit.aes = FALSE) +
           ylim(-1.7, 0) +
           scale_x_discrete(breaks = plot_data_labels$xpos, expand = expansion(add = 1.0)) +
           #scale_x_discrete(labels = data_by_density$Density, drop = FALSE) +   # leave space for missing columns
           scale_fill_gradientn(
            colours = c("#ffffc0", "#c80000"),
            ###colors = c("#fcfbf1", "#C80000"), # Define color transitions. #220908, #fcfbf1,"blue", "yellow", "#FF4000"
            values = sqrt(c(0, (HeatMap_threshold - min_fc) / (max_fc - min_fc), 1)), # Normalize transitions
            limits = sqrt(c(min_fc, max_fc)), # Define range of the scale
            na.value = "skyblue",
            breaks = sqrt(c(HeatMap_threshold, 1, 10, 20, 50, 100, 250, 500, 1000, 2000, max_fc)), # Add new breakpoints
            labels = c(
              #as.character(round(min_fc, 1)), 
               as.character(HeatMap_threshold), 
               "1", "10", "20", "50", "100", "250", "500", "1000", "2000", 
               as.character(round(max_fc, 1))
             ), # Add labels corresponding to the breaks
            guide = guide_colorbar(
               ticks = TRUE, 
              ticks.linewidth = 0.5, 
               ticks.colour = "black", 
              label.position = "bottom", 
              barwidth = 15, 
             barheight = 0.8
            )
           ) +
       labs(x = NULL, y = NULL, fill = 'Fold Change') + # Label the legend
       #geom_vline(mapping = aes(xintercept = xmin - 0.5, colour = "grey70", linewidth = 0.25), data = plot_data_labels) +
       geom_vline(mapping = aes(xintercept = xmin - 0.5), data = plot_data_labels, show.legend = FALSE) +
       theme_minimal() +
       theme(
         #axis.ticks.x = element_line(linewidth = 0.5),
         axis.ticks.x = element_blank(),
         axis.text.x = element_text(size = 5.2, angle = 90, hjust = 1, vjust = 0.5),
         axis.text.y = element_blank(), # Hide y-axis for heatmap
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
         panel.grid = element_blank(),
         plot.margin = margin(t = 0, b = 0),  # remove top and bottom margins
         legend.position = "bottom",
         legend.direction = "horizontal", # Ensure the legend is horizontal
        legend.title = element_text(size = 7, family = font.family), # Customize legend title font size and family
        legend.text = element_text(size = 7, family = font.family)
       )


g.labels <- ggplot(data = plot_data_filtered, aes(x = log2FC, y = log2FC)) +
				scale_x_continuous(breaks = as.numeric(plot_data_labels$xpos), 
								   minor_breaks = sort(as.numeric(plot_data_labels$xmin))[-1] - 0.5,
								   labels = plot_data_labels$Decoration, 
								   limits = c(min(as.numeric(plot_data_filtered$xlabel)), max(as.numeric(plot_data_filtered$xlabel))),
								   expand = expansion(add = 0.6),
								   guide=guide_axis(minor.ticks = TRUE)
								  ) +
				labs(x = 'Glycans', y = NULL) +
				theme_minimal() +
                theme(
					axis.text.x = element_text(size = 7.0, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
					axis.title.x = element_text(size = 7),
					axis.text.y = element_blank(),
					#axis.ticks.x = element_line(linewidth = 0.5),
					#axis.minor.ticks.x.bottom = element_line(linewidth = 0.5),
					#axis.minor.ticks.length = rel(1.0),
					#axis.minor.ticks.length.x.bottom = unit(5, "pt"),
					#panel.grid = element_blank(),
					plot.margin = margin(t = 0)
					)



# Stack the plots vertically, aligning the x-axes
g.bar.grobs <- ggplotGrob(g.bar)
g.heatmap.grobs <- ggplotGrob(g.heatmap)
g.labels.grobs <- ggplotGrob(g.labels)

# exchange guides between plots
idx.heatmap <- which(g.heatmap.grobs$layout$name == "guide-box-bottom")
idx.labels  <- which(g.labels.grobs$layout$name == "guide-box-bottom")

top.row.heatmap <- g.heatmap.grobs$layout$t[idx.heatmap]
top.row.labels  <- g.labels.grobs$layout$t[idx.labels]

# swap the gap between the x-labels and the guides
temp <- g.labels.grobs$heights[top.row.labels - 1]
g.labels.grobs$heights[top.row.labels - 1] <- g.heatmap.grobs$heights[top.row.heatmap - 1]
g.heatmap.grobs$heights[top.row.heatmap - 1] <- temp

# swap heights of the guides
temp <- g.labels.grobs$heights[top.row.labels]
g.labels.grobs$heights[top.row.labels] <- g.heatmap.grobs$heights[top.row.heatmap]
g.heatmap.grobs$heights[top.row.heatmap] <- temp

# swap the guides
temp <- g.labels.grobs$grobs[[idx.labels]]
g.labels.grobs$grobs[[idx.labels]] <- g.heatmap.grobs$grobs[[idx.heatmap]]
g.heatmap.grobs$grobs[[idx.heatmap]] <- temp



g.combined <- rbind(g.bar.grobs, g.heatmap.grobs, g.labels.grobs)
g.combined$widths <- unit.pmax(g.bar.grobs$widths, g.heatmap.grobs$widths)
g.combined$heights[9] <- unit(10.0, "null")  # make the height of the bar chart panel large relative to that of the heatmap
g.combined$heights[41] <- unit(0.0, "null")  # make the height of the bar chart panel large relative to that of the heatmap

# Save the combined plot as EPS
ggsave(paste0(compare_dir_name, "2.pdf"), plot = g.combined, width = 15, height = 7)



###
###
###
###
        plots <- list()
        for(i in 1:length(rownames(dge$samples))){
             row.name <- rownames(dge$samples)[i]
             group.name <- dge$samples$group[i]
             raw.column <- compare.table[row.name][,1]   # second index strips colname
             if(group.name == "Test"){
                 ref.column <- compare.table$"Test (CPM)"
             } else {
                 ref.column <- compare.table$"Reference (CPM)"
             }
             raw.mask <- raw.column > 0    # remove zero counts for log-scale
             subplot <- ggplot(data.frame(x = ref.column[raw.mask], 
                                          y = raw.column[raw.mask])) +
                            geom_point(aes(x = x, y = y), size = 1) + 
                            scale_x_log10() + 
                            scale_y_log10() +
                            labs(x = row.name, y = "") +
                            theme(text = element_text(family = font.family, 
                                                      size = 7))
             #print(subplot)
             plots[[i]] <- ggplotGrob(subplot)
        }
        ###do.call(grid.arrange, plots)
        for(i in 1:ceiling(length(plots)/12)){
            do.call(grid.arrange, c(plots[(12*i-11):min(12*i,length(plots))],
                                    ncol = 4, nrow = 3))
        }






        # use knitr to produce the PDF report
        ### render("report.Rmd", output_dir = compare_dir_name, 
        ###                     intermediates_dir = compare_dir_name,
        ###	###output_format = "pdf_document",
        ###       params = list(testing_data = testing_data,
        ###                     difference.tests = diff.tests,
        ###                     main.table = long.data))
        #                     bar.plot = bar.plot[[1]]))

}



draw.heatmap <- function(long.data, aspect.ratio, horiz.size){

    # compress range of data to 1+ and set new colour scheme
    colour_steps2 <- c("white", "#E8E8B4", "#C08080", "#744646", "#1D0000")
    colour_intervals <- c(0, 0.1, 0.568, 0.856, 1)

    colour_steps2 <- c(colour_steps2, "#500000", "#500000",
                              "#680000", "#680000", "#800000", "#800000")
    colour_intervals <- c(colour_intervals * (sqrt(40)-1)/(sqrt(100)-1),
                                             (sqrt(40)-1)/(sqrt(100)-1)+0.0001,
             (sqrt(60)-1)/(sqrt(100)-1), (sqrt(60)-1)/(sqrt(100)-1)+0.0001,
             (sqrt(80)-1)/(sqrt(100)-1), (sqrt(80)-1)/(sqrt(100)-1)+0.0001, 1)

    colour_steps2 <- c(colour_steps2, "#500050", "#500050",
                              "#480080", "#480080", "#2000A0", "#2000A0")
    colour_intervals <- c(colour_intervals * (sqrt(40)-1)/(sqrt(70)-1),
                                        (sqrt(40)-1)/(sqrt(70)-1)+0.0001,
             (sqrt(50)-1)/(sqrt(70)-1), (sqrt(50)-1)/(sqrt(70)-1)+0.0001,
             (sqrt(60)-1)/(sqrt(70)-1), (sqrt(60)-1)/(sqrt(70)-1)+0.0001, 1)

    ### matplotlib.colors.rgb2hex(c) for c in matplotlib.cm.get_cmap(name = 'viridis').colors]
    colour_steps2 <- c('#440154', '#440256', '#450457', '#450559', 
            '#46075a', '#46085c', '#460a5d', '#460b5e', '#470d60', '#470e61', 
            '#471063', '#471164', '#471365', '#481467', '#481668', '#481769', 
            '#48186a', '#481a6c', '#481b6d', '#481c6e', '#481d6f', '#481f70', 
            '#482071', '#482173', '#482374', '#482475', '#482576', '#482677', 
            '#482878', '#482979', '#472a7a', '#472c7a', '#472d7b', '#472e7c', 
            '#472f7d', '#46307e', '#46327e', '#46337f', '#463480', '#453581', 
            '#453781', '#453882', '#443983', '#443a83', '#443b84', '#433d84', 
            '#433e85', '#423f85', '#424086', '#424186', '#414287', '#414487', 
            '#404588', '#404688', '#3f4788', '#3f4889', '#3e4989', '#3e4a89', 
            '#3e4c8a', '#3d4d8a', '#3d4e8a', '#3c4f8a', '#3c508b', '#3b518b', 
            '#3b528b', '#3a538b', '#3a548c', '#39558c', '#39568c', '#38588c', 
            '#38598c', '#375a8c', '#375b8d', '#365c8d', '#365d8d', '#355e8d', 
            '#355f8d', '#34608d', '#34618d', '#33628d', '#33638d', '#32648e', 
            '#32658e', '#31668e', '#31678e', '#31688e', '#30698e', '#306a8e', 
            '#2f6b8e', '#2f6c8e', '#2e6d8e', '#2e6e8e', '#2e6f8e', '#2d708e', 
            '#2d718e', '#2c718e', '#2c728e', '#2c738e', '#2b748e', '#2b758e', 
            '#2a768e', '#2a778e', '#2a788e', '#29798e', '#297a8e', '#297b8e', 
            '#287c8e', '#287d8e', '#277e8e', '#277f8e', '#27808e', '#26818e', 
            '#26828e', '#26828e', '#25838e', '#25848e', '#25858e', '#24868e', 
            '#24878e', '#23888e', '#23898e', '#238a8d', '#228b8d', '#228c8d', 
            '#228d8d', '#218e8d', '#218f8d', '#21908d', '#21918c', '#20928c', 
            '#20928c', '#20938c', '#1f948c', '#1f958b', '#1f968b', '#1f978b', 
            '#1f988b', '#1f998a', '#1f9a8a', '#1e9b8a', '#1e9c89', '#1e9d89', 
            '#1f9e89', '#1f9f88', '#1fa088', '#1fa188', '#1fa187', '#1fa287', 
            '#20a386', '#20a486', '#21a585', '#21a685', '#22a785', '#22a884', 
            '#23a983', '#24aa83', '#25ab82', '#25ac82', '#26ad81', '#27ad81', 
            '#28ae80', '#29af7f', '#2ab07f', '#2cb17e', '#2db27d', '#2eb37c', 
            '#2fb47c', '#31b57b', '#32b67a', '#34b679', '#35b779', '#37b878', 
            '#38b977', '#3aba76', '#3bbb75', '#3dbc74', '#3fbc73', '#40bd72', 
            '#42be71', '#44bf70', '#46c06f', '#48c16e', '#4ac16d', '#4cc26c', 
            '#4ec36b', '#50c46a', '#52c569', '#54c568', '#56c667', '#58c765', 
            '#5ac864', '#5cc863', '#5ec962', '#60ca60', '#63cb5f', '#65cb5e', 
            '#67cc5c', '#69cd5b', '#6ccd5a', '#6ece58', '#70cf57', '#73d056', 
            '#75d054', '#77d153', '#7ad151', '#7cd250', '#7fd34e', '#81d34d', 
            '#84d44b', '#86d549', '#89d548', '#8bd646', '#8ed645', '#90d743', 
            '#93d741', '#95d840', '#98d83e', '#9bd93c', '#9dd93b', '#a0da39', 
            '#a2da37', '#a5db36', '#a8db34', '#aadc32', '#addc30', '#b0dd2f', 
            '#b2dd2d', '#b5de2b', '#b8de29', '#bade28', '#bddf26', '#c0df25', 
            '#c2df23', '#c5e021', '#c8e020', '#cae11f', '#cde11d', '#d0e11c', 
            '#d2e21b', '#d5e21a', '#d8e219', '#dae319', '#dde318', '#dfe318', 
            '#e2e418', '#e5e419', '#e7e419', '#eae51a', '#ece51b', '#efe51c', 
            '#f1e51d', '#f4e61e', '#f6e620', '#f8e621', '#fbe723', '#fde725')
    colour_steps2 <- colour_steps2[54:256]
    colour_intervals <- 0:(length(colour_steps2) - 1) / (length(colour_steps2) - 1)


    long.data$Type <- factor(long.data$Type, levels = rev(unique(long.data$Type)))


    ### # throw away the Naive comparisons and recalculate the FDR
    ### long.data <- long.data[long.data$Reference != "Naive",]
    ### long.data$QValue <- p.adjust(long.data$PValue, method = "BH")

    ### long.data$logFC[long.data$logFC < 0] <- 0
    ### long.data$logFC[long.data$logFC > log2(30)] <- log2(30)

    heatmap <- ggplot(long.data, #[long.data$Label %in% unique(long.data$Label[long.data$Valid]),],
                  aes(x = Label, y = Type)) +
      theme_light() + theme(text = element_text(family = font.family),
                            panel.grid.major.x = element_blank(),
                            panel.grid.major.y = element_blank()) +
      coord_fixed(aspect.ratio) +
      geom_tile(aes(fill = 2^logFC), colour = "gray70", linewidth = 0.2343) +
      #geom_hline(yintercept = 0.5 + 2 * (0:11), colour = "black", size = 0.7) + 
      scale_fill_gradientn(colours = colour_steps2, values = colour_intervals,
                       limits = c(0, 75), trans = "sqrt",
                       ###limits = c(0.1, 65), trans = "log",
                       #breaks = c(0:8, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100),
                       #labels = as.character(c(0:8, 2:7 * 5, 4:9 * 10, " 100+")),
                       breaks = c(0.1, 0.5, 1:10, 3:15*5),
                       labels = as.character(c(0.1, 0.5, 1:10, 3:15 * 5)),
                       na.value = "white", guide = guide_colorbar(nbin = 1000)) +
      ###geom_text(label = "*", hjust = 0.5, vjust = 0.80, color = "#808080",
      ###          data = long.data[significant.rows.2,]) +
      ###geom_text(label = "\u25e6", hjust = 0.5, vjust = 0.40, color = "#808080",
      ###          data = long.data[significant.rows.1,]) +

      labs(x = "Glycan", y =  "Binding Target", aes(hjust = 0)) +
      theme(axis.text = element_text(size = 7.0, colour = "black"),
            axis.title = element_text(size = 7.0, colour = "black"), 
            axis.text.x  = element_text(hjust = 1, vjust = 0.5, angle = 90, size = horiz.size),
            legend.position = "bottom",
            legend.title = element_blank(), # element_text(size = 10),
            legend.text = element_text(size = 7.0),
            legend.key.size = unit(0.3, "cm"),
            legend.key.width = unit(2.5, "cm"))
            
      heatmap
}

###long.data$logFC <- pmax(log(0.1), long.data$logFC)   # clamp data for heatmaps

#heatmap <- draw.heatmap(long.data, 3, 3.8)
heatmap <- draw.heatmap(long.data, 1, 7.0)
print(heatmap)

xpos <- long.data$xpos
#heatmap <- draw.heatmap(long.data[xpos <= 34 | (127 <= xpos & xpos <= 146),], 1, 7)
heatmap <- draw.heatmap(long.data[xpos <= 21 | (34 <= xpos & xpos <= 40),], 1, 7)
print(heatmap)

##heatmap <- draw.heatmap(long.data[(35 <= xpos & xpos <= 92) | 127 <= xpos,], 1.62, 6.0)
heatmap <- draw.heatmap(long.data[22 <= xpos,], 1, 7.0)
print(heatmap)

# end the PDF file generation
dev.off()
