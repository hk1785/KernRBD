\name{biom}
\alias{biom}
\docType{data}
\title{
Data
}
\description{
Gut microbiome data to assess the effect of periodically restricted feeding (PRF) on gut microbiome profiles (Yanai et al., 2024, Nat Commun)
}
\usage{data("biom")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : num [1:386, 1:48] 0 0 0 0 8 0 0 0 0 0 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:386] "e920f96de1114d11639e145ceaddfaa3" "087a0803b5a078c7c547580b445f0903" "8e462eb3179eaa016862fb32cc78858a" "b04c81f410c6a2b7463350dd904fcbcb" ...
  .. .. .. .. ..$ : chr [1:48] "S1" "S2" "S3" "S4" ...
  .. .. ..@ taxa_are_rows: logi TRUE
  .. .. ..$ dim     : int [1:2] 386 48
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ : chr [1:386] "e920f96de1114d11639e145ceaddfaa3" "087a0803b5a078c7c547580b445f0903" "8e462eb3179eaa016862fb32cc78858a" "b04c81f410c6a2b7463350dd904fcbcb" ...
  .. .. .. ..$ : chr [1:48] "S1" "S2" "S3" "S4" ...
  ..@ tax_table:Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
  .. .. ..@ .Data: chr [1:386, 1:7] "k_Bacteria" "k_Bacteria" "k_Bacteria" "k_Bacteria" ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:386] "e920f96de1114d11639e145ceaddfaa3" "087a0803b5a078c7c547580b445f0903" "8e462eb3179eaa016862fb32cc78858a" "b04c81f410c6a2b7463350dd904fcbcb" ...
  .. .. .. .. ..$ : chr [1:7] "Kingdom" "Phylum" "Class" "Order" ...
  .. .. ..$ dim     : int [1:2] 386 7
  .. .. ..$ dimnames:List of 2
  .. .. .. ..$ : chr [1:386] "e920f96de1114d11639e145ceaddfaa3" "087a0803b5a078c7c547580b445f0903" "8e462eb3179eaa016862fb32cc78858a" "b04c81f410c6a2b7463350dd904fcbcb" ...
  .. .. .. ..$ : chr [1:7] "Kingdom" "Phylum" "Class" "Order" ...
  ..@ sam_data :'data.frame':	48 obs. of  2 variables:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 2
  .. .. .. ..$ : chr [1:48] "DCGG" "DCGG" "DCGG" "DCGG" ...
  .. .. .. ..$ : chr [1:48] "Baseline" "PRF 1" "PRF 2" "PRF 3" ...
  .. .. ..@ names    : chr [1:2] "Monkey.ID" "Treatment"
  .. .. ..@ row.names: chr [1:48] "S1" "S2" "S3" "S4" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree :List of 5
  .. ..$ edge       : int [1:770, 1:2] 387 388 388 387 389 390 391 391 392 392 ...
  .. ..$ edge.length: num [1:770] 0.0383 0.0789 0.1279 0.0559 0.1253 ...
  .. ..$ Nnode      : int 385
  .. ..$ node.label : chr [1:385] "" "0.752" "0.806" "0.956" ...
  .. ..$ tip.label  : chr [1:386] "e920f96de1114d11639e145ceaddfaa3" "087a0803b5a078c7c547580b445f0903" "8e462eb3179eaa016862fb32cc78858a" "b04c81f410c6a2b7463350dd904fcbcb" ...
  .. ..- attr(*, "class")= chr "phylo"
  .. ..- attr(*, "order")= chr "cladewise"
  ..@ refseq   : NULL
}
\references{
Yanai et al. Short-term periodic restricted feeding elicits metabolome-microbiome signatures with sex dimorphic persistence in primate intervention. Nat Commun, 15(1088), 2024.
}
\examples{
data(biom)
}
\keyword{datasets}
