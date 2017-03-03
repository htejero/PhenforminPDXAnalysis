#Script to obtain the Phenformin Signature from LINCS API and Raw Data 
## Hector Tejero 2016


# Drug 
N_genes = 250

limit = 200

# Get drugs distil_id from database

source("~/lincs/code/lincsApi/getGenesFromProbesDbi.R")

api_key = ""

#target =   "" "Phenformin%20hydrochloride"

 pert_desc = "Phenformin%20hydrochloride" # "PIM1" # "BRD-A91699651" #CHLOROQUINE
# 
# pert_id = ""

require('rjson')

url = paste("http://api.lincscloud.org/a2/instinfo?q={%22pert_id%22:%22", pert_id, "%22}&l=", limit,"&user_key=", api_key, sep = "")

#
url = paste("http://api.lincscloud.org/a2/instinfo?q={%22pert_desc%22:%22", pert_desc, "%22}&l=", limit,"&user_key=", api_key, sep = "")


id.list = fromJSON(file=url)

distil_id = sapply(id.list, function(x) x$distil_id)

det_plate = sapply(id.list, function(x) x$det_plate)

target_info = data.frame(distil_id = sapply(id.list, function(x) x$distil_id),
                  pert_desc =   as.character(sapply(id.list, function(x) x$pert_desc) ),
                  pert_id =  sapply(id.list, function(x) x$pert_id),
                  cell_id = sapply(id.list, function(x) x$cell_id),
                  det_plate = sapply(id.list, function(x) x$det_plate),
                  pert_time = sapply(id.list, function(x) x$pert_time),
                  pert_dose = sapply(id.list, function(x) x$pert_dose))


#  Get DMSO distil_id 

dmso.info = data.frame( distil_id = character(),                        
                        pert_desc = character(),
                        pert_id = character(),
                        cell_id = character(), 
                        det_plate = character(),
                        pert_time = character(), 
                        pert_dose = character() )

for (det_id in det_plate) {
  
  url = paste("http://api.lincscloud.org/a2/instinfo?q={%22pert_id%22:%22", "DMSO", "%22,%22det_plate%22:%22", det_id,"%22}&l=", limit,"&user_key=", api_key, sep = "")
  
  det_id.list = fromJSON(file=url)
  
  line = data.frame(distil_id = sapply(det_id.list, function(x) x$distil_id),
  pert_desc =   as.character(sapply(det_id.list, function(x) x$pert_desc)),
  pert_id =  sapply(det_id.list, function(x) x$pert_id),
  cell_id = sapply(det_id.list, function(x) x$cell_id),
  det_plate = sapply(det_id.list, function(x) x$det_plate),
  pert_time = sapply(det_id.list, function(x) x$pert_time),
  pert_dose = sapply(det_id.list, function(x) x$pert_dose))

  dmso.info = rbind(dmso.info, line)
  
}


# Make Table 


data = rbind(dmso.info, target_info)

file.data = paste("~/lincs/data/", pert_desc, "_", pert_id, "_DMSO_table.csv", sep = "")

 write.csv(data, file.data)

# Get signatures 

source("~/lincs/code/lincsC3cloud/io.R")

sig_ids = as.character(data$distil_id)
file = "~/LINCS_RawData/q2norm_n1328098x22268.gctx"  #Level 3 Rawdata


# 2. initialize a GCT object with the signatures

ds <- parse.gctx(file, cid= sig_ids) 

x = ds@mat  #get the matrix with the expression values 

H5close()

colnames(x)==data$distil_id


# Make Model 

library(limma)

#design <- model.matrix( ~ pert_desc + cell_id + det_plate + pert_desc*cell_id, data = data ) #Interaction


design <- model.matrix( ~ pert_id + cell_id + det_plate , data = data )


colnames(design)[2] = "Phen"

fit <- lmFit(x, design)

fit2 <- eBayes(fit)

TT = topTable(fit2, coef = "Phen", adjust="BH", number=dim(x)[1])  

#topTableF ranks genes on the basis of moderated F-statistics for one or more coefficients. If topTable is called and coef has two or more elements, then the specified columns will be extracted from fit and topTableF called on the result. topTable with coef=NULL is the same as topTableF, unless the fitted model fit has only one column.


probes_DN = rownames(TT)[head(order(TT$logFC), N_genes)]
probes_UP = rownames(TT)[tail(order(TT$logFC), N_genes)]


genes_dn = as.character(getGenesFromProbesDbi(probes_DN, genes2probesfile= "~/lincs/code/lincsApi/genes2probes.txt", api_key = api_key))

genes_up = as.character(getGenesFromProbesDbi(probes_UP, genes2probesfile= "~/lincs/code/lincsApi/genes2probes.txt", api_key = api_key))


#pert_desc = "Phenformin_hydrochloride"

name_UP = paste(pert_desc, "_UP", sep = "")
name_DN = paste(pert_desc, "_DN", sep = "")

desc = "level3_limma_model_generated"


library(GSEABase)


GScollect = GeneSetCollection(GeneSet(setName = name_UP, shortDescription = desc , genes_up), 
                              GeneSet(setName = name_DN, shortDescription = desc, genes_dn))


file = paste("~/LINCS_RawData/Signatures/", pert_desc, "_", pert_id, "_Limma_model_level3_N_", N_genes, ".gmt", sep = "")

toGmt(GScollect, file )


