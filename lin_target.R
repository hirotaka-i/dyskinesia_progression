args <- commandArgs(trailingOnly = TRUE)
t = strsplit(args, ";")[[1]]
MODEL=t[1];OUTCOME=t[2];COVPC=t[3];SLICE=t[4];PHENO=t[5];OUTPUT=t[6]
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
COVs = strsplit(COVPC, "\\+")[[1]]
data = fread(PHENO) 
SNPset = fread(SLICE)
SNPs = setdiff(names(SNPset), c('IID', 'GBA_N_A'))
cohort = data %>% mutate(Sex = if_else(Sex=='Male', 0, 1)) %>% arrange(IID) %>% data.frame()
cohort[COVs] = as.data.frame(scale(cohort[COVs]))
cohort_snp = inner_join(cohort, SNPset, by = "IID")
# ANALYSIS
test.listfunc = function(x){
  # Models
  MODEL = paste0(OUTCOME, "~" , "`", SNPs[x], "` + ", COVPC)
  testModel = try(glm(eval(parse(text = MODEL)), data = cohort_snp),silent = T)
  if(class(testModel)[1]=="try-error"){
    sumstat=c(SNPs[x], "NoConverge", rep(NA,5))
  }else{
    temp= summary(testModel)$coefficients
    v_interest = substr(SNPs[x],1,3)
    if(grep(v_interest, rownames(temp)) %>% length == 0){ # In this case, SNP is dropeed from the model
      sumstat=c(SNPs[x], "NoVforSNP", rep(NA, 5))
    }else{
      RES = temp[2,]
      s = cohort_snp[,c("IID", SNPs[x])] %>% distinct(IID, .keep_all = T) %>% filter(!is.na(SNPs[x]))
      sumstat <- c(SNPs[x], RES[3], RES[1], RES[2], RES[4], length(testModel$y), mean(s[,SNPs[x]])/2)
    }
  }
  return(sumstat)
}

temp = lapply(1:length(SNPs), test.listfunc)
temp2 = do.call(rbind, temp) %>% data.frame # %>%filter(complete.cases(.))
names(temp2)=c("POS_A2_A1", "Tvalue", "BETA", "SE", "P", "N", "ALT_Frq") 
temp3 = temp2 %>% separate(POS_A2_A1, c("SNP", "A2", "A1"),sep="_")
dir.create(OUTPUT, recursive = T, showWarnings = F)
FILENAME2 = sub('_', '.', basename(PHENO)) %>% sub('lin.txt', paste(MODEL, OUTCOME, 'lin.txt', sep='.'), .)
write.table(temp3, paste(OUTPUT, FILENAME2, sep="/"), row.names = F, quote = F, sep = "\t")
