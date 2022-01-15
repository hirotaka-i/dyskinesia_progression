args <- commandArgs(trailingOnly = TRUE)
t = strsplit(args, ";")[[1]]
MODEL=t[1];OUTCOME=t[2];COVPC=t[3];SLICE=t[4];PHENO=t[5];OUTPUT=t[6]
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))
COVs = strsplit(COVPC, "\\+")[[1]]
data = fread(PHENO) 
SNPset = fread(SLICE)
SNPs = setdiff(names(SNPset), c('IID', 'GBA_N_A'))
cohort = data %>% mutate(Sex = if_else(Sex=='Male', 0, 1)) %>% arrange(IID, TSTART) %>% data.frame()
cohort[COVs] = as.data.frame(scale(cohort[COVs]))
cohort_snp = inner_join(cohort, SNPset, by = "IID")
cohort_snp$SurvObj1 = with(cohort_snp, Surv(Disease_duration, Dyskinesia))
# ANALYSIS
test.listfunc = function(x){
  # Models
  MODEL = paste("SurvObj1~" , "`", SNPs[x], "` + ", COVPC, sep = "")
  testCox = try(coxph(eval(parse(text = MODEL)), data = cohort_snp),silent = T)
  if(class(testCox)[1]=="try-error"){
    sumstat=c(SNPs[x], "NoConverge", rep(NA,6))
  }else{
    temp= summary(testCox)$coefficients
    if(grep(substr(SNPs[x],1,3), rownames(temp)) %>% length == 0){ # In this case, SNP is dropeed from the model
      sumstat=c(SNPs[x], "NoVforSNP", rep(NA, 6))
    }else{
      RES = temp[1,]
      EVENT_OBS = paste(testCox$nevent, testCox$n, sep="_")
      s = cohort_snp[,c("IID", SNPs[x])] %>% distinct(IID, .keep_all = T) %>% filter(!is.na(SNPs[x]))
      sumstat <- c(SNPs[x], EVENT_OBS, as.numeric(RES[4]), RES[1], RES[3], RES[5], nrow(s), mean(s[,SNPs[x]])/2)
    }
  }
  return(sumstat)
}

temp = lapply(1:length(SNPs), test.listfunc)
temp2 = do.call(rbind, temp) %>% data.frame # %>%filter(complete.cases(.))
names(temp2)=c("POS_A2_A1", "EVENT_OBS", "Tvalue", "BETA", "SE", "P", "N", "ALT_Frq") 
temp3 = temp2 %>% separate(POS_A2_A1, c("SNP", "A2", "A1"),sep="_")
dir.create(OUTPUT, recursive = T, showWarnings = F)
FILENAME2 = sub('_', '.', basename(PHENO)) %>% sub('surv.txt', paste(MODEL, OUTCOME, 'cox.txt', sep='.'), .)
write.table(temp3, paste(OUTPUT, FILENAME2, sep="/"), row.names = F, quote = F, sep = "\t")
