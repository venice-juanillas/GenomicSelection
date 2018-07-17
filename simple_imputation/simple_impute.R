
## optparse needed for taking care of command line arguments
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="hapmap file", metavar="gobii_hapmap_file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file name", metavar="output_file_name"),
  make_option(c("-i", "--imputeMethod"), type="character", default="mean", 
              help="Imputation method: mean or mode", metavar="impute_method"),
  make_option(c("-n", "--numSkip"), type="integer", default="36", 
              help="how many rows to skip", metavar="rows_to_skip")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## supply at least the input file 
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

## read the hapmap file
geno<-read.delim(opt$f,sep="\t",quote="\"",na.strings="NA", skip=opt$n)

## store the metadata headers in a different variable
metadata <- read.delim(opt$f, nrows=opt$n)

## write the metadata headers to the file
write.table(metadata, file = opt$o, append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".",  row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"))

## subset the matrix.
## imputation will start at column 11
M <- geno[,-1:-11]

## impute with population mean or mode 
if(opt$i=="mean"){
  
  M <- apply(M, 1, function(x){ y = mean(x, na.rm = TRUE); x[which(is.na(x))] <- y; x })
  print("missing values are imputed with population mean\n")
}

if(opt$i =="mode"){
  modeX <- function(x){
    md <- as.numeric(names(table(x))[which.max(table(x))])
    return(md)
  }
  M <- apply(M, 1, function(x){ y = modeX(x); x[which(is.na(x))] <- y; x })
  print("missing values are imputed with population mode\n")
}

## reconstruct data frame after imputation.
## transpose since prior imputation step have changed the orientation of the matrix
geno <- data.frame (geno[,1:11], t(M), check.names = FALSE)
geno1 <- geno[,-ncol(geno)]

## write the output into the out_file. 
## will make append = TRUE since we need to append this table right after the prior metadata var
write.table(geno1, file = opt$o, append = TRUE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"))

########################################## end of simple imputation ##############################################
