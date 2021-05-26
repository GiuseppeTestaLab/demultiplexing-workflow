
library(scDblFinder)
library(DropletUtils)



filtered_matrix <- snakemake@params[["filtered_matrix"]]
barcodesclass <- snakemake@params[["barcodesclass"]]

print(filtered_matrix)



## __2. Read data__


sce <- read10xCounts(filtered_matrix)

## __3. scDBLfinder__


sce <- scDblFinder(sce)

calls <- DataFrame(scDBLfinder=sce$scDblFinder.class, row.names=sce$Barcode)



## __4. Writing result__

write.table(calls, barcodesclass, sep="\t", col.names=F, quote=F)
