
### Reproducible code for the 
## Potency of phenformin revealed by an unbiased study of metabolic inhibitors in a 
## large  set of pancreatic cancer patient-derived xenografts

## Hector Tejero, Cancer Bioinformatics Unit. National Spanish Cancer Center 2016-2017 

require(GSVA)
require(GSEABase)
require("hgu133plus2.db")
require("AnnotationDbi")
require(ggplot2)


## Download and unzip the expression file from GSE51798 project

# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51798/suppl/GSE51798_PDX_and_PT_data_plus_reprocessed_data.txt.gz

path = "" # include the path where the expression file is 

expression.file = file.path(path, "GSE51798_PDX_and_PT_data_plus_reprocessed_data.txt")  

expression.data = read.table(expression.file, 
               header = TRUE)

## Downlaod samples_phen_met.tsv from https://github.com/htejero/PhenforminPDXAnalysis in the path directory

response.data =read.csv(file.path( path, "samples_phen_met.tsv"), header = TRUE, stringsAsFactors = FALSE)

response.data = response.data[!is.na(response.data$Phenformin) & !(is.na(response.data$Metformin)),]

## Filter expression data 

X = expression.data[, c(1, which(colnames(expression.data) %in% response.data$Name))]

## Annotate to genes 

X$ID_REF = as.character(X$ID_REF)

idx = sapply(as.character(X$ID_REF), function(x) exists(x, hgu133plus2SYMBOL))

genes.dbi = unlist(mget(X$ID_REF[idx], hgu133plus2SYMBOL))

genes.dbi = genes.dbi[!is.na(genes.dbi)]

es = merge(X, genes.dbi, by.x = "ID_REF", by.y = 0, all.x = TRUE)

colnames(es)[ncol(es)] = "Gene.Symbol"

es =  aggregate(. ~ Gene.Symbol , data = data.frame(Gene.Symbol = es$Gene.Symbol,
                                                    data.matrix(es[, -c(1,ncol(es))]))
                , mean)

## Find the Expression of the signatures in the sample by ssGSEA

GSC = getGmt(file.path(path, "Phenformin_Limma_model_level3.gmt"))

genes = es[,1]

es = data.matrix(es[,-c(1)])

rownames(es) = genes


GSVA.result <- gsva(es , geneIds(GSC), verbose = TRUE, method = "ssgsea", mx.diff = TRUE, rnaseq = FALSE)  #lo de geneIds es por el fallo que aparece en el link de arriba 

ES.phen = GSVA.result[grep("_UP", rownames(GSVA.result)),]-GSVA.result[grep("_DN", rownames(GSVA.result)),]

data.ES = merge(response.data, ES.phen, by.x = "Name", by.y = 0, all.x = T)

colnames(data.ES)[6] = "ES"


## Remove JH033 due to dominant PALB2 mutation

data.ES = data.ES[data.ES$PDX!="JH033",]

## Correlation and tets

cor.test(data.ES$Phenformin, data.ES$ES)

# Plot the data 


ggplot(data.ES, aes(x = ES, y = Phenformin)) + geom_point(size = 3, shape = 21) + 
  xlab("Phenformin Signature (ssGSEA)") + ylab("Tumour Growth Inhibition (%)") +
  theme(legend.position="bottom") + geom_text(aes(label=PDX), hjust=-.4, vjust=-.1, size = 3)  + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylim(0,50) +
  xlim(-0.9, -0.2)

