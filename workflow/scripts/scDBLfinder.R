



library(scDblFinder)
library(DropletUtils)



filtered_matrix <- snakemake@input[["filtered_matrix"]]
barcodeList <- snakemake@input[["barcodes"]]
barcodesclass <- snakemake@output[["barcodesclass"]]
scDblFinderpath <- snakemake@params[["scDblFinderpath"]]



dir.create(scDblFinderpath, showWarnings = FALSE)


## __2. Read data__


sce <- read10xCounts(filtered_matrix)
barcodeMask <-  read.csv(barcodeList, header=FALSE)$V1
sce <- sce[,sce$Barcode %in% barcodeMask]


## __3. scDBLfinder__


sce <- scDblFinder(sce)

calls <- DataFrame(scDBLfinder=sce$scDblFinder.class, row.names=sce$Barcode)



## __4. Writing result__

write.table(calls, barcodesclass, sep="\t", col.names=F, quote=F)
