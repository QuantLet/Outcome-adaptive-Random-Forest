# observational study
Hmisc::getHdata(rhc,where="https://hbiostat.org/data/repo/")
dim(rhc) # number of obs
summary(rhc$swang1) # number of treated vs. untreated
rhc <- rhc[,names(rhc)[order(names(rhc))]] # order column names alphabetically
for (cc in 1:ncol(rhc)) {
  cat(names(rhc)[cc], "----------------------------------------------------\n")
  print(summary(rhc[,cc]))
}
rhc$cat2[is.na(rhc$cat2)] <- "Missing" # remove NAs from cat2
## relabelled based on Hirano and Imbens (2002) Table 1
rhc.dt <- data.frame(
  "A"=rhc$swang1=="RHC",
  "Y"=rhc$dth30=="No",
  "age"=rhc$age,
  "sex"=rhc$sex=="Female",
  # dummy variables for race
  "raceblack"=rhc$race=="black",
  "raceother"=rhc$race=="other",
  "edu"=rhc$edu,
  # dummy variables for income
  "income1"=rhc$income=="$11-$25k",
  "income2"=rhc$income=="$25-$50k",
  "income3"=rhc$income=="> $50k",
  # dummy variables for Medical insurance
  "ins_care"=rhc$ninsclas=="Medicare",
  "ins_pcare"=rhc$ninsclas=="Private & Medicare",
  "ins_caid"=rhc$ninsclas=="Medicaid",
  "ins_no"=rhc$ninsclas=="No insurance",
  "ins_carecaid"=rhc$ninsclas=="Medicare & Medicaid",
  # dummy variables for Primary disease category
  "cat1_copd"=rhc$cat1=="COPD",
  "cat1_mosfsep"=rhc$cat1=="MOSF w/Sepsis",
  "cat1_mosfmal"=rhc$cat1=="MOSF w/Malignancy",
  "cat1_chf"=rhc$cat1=="CHF",
  "cat1_coma"=rhc$cat1=="Coma",
  "cat1_cirr"=rhc$cat1=="Cirrhosis",
  "cat1_lung"=rhc$cat1=="Lung Cancer",
  #"cat1_colon"=rhc$cat1=="Colon Cancer",
  # dummy variables for Secondary disease category
  "cat2_mosfsep"=rhc$cat2=="MOSF w/Sepsis",
  "cat2_coma"=rhc$cat2=="Coma",
  "cat2_mosfmal"=rhc$cat2=="MOSF w/Malignancy",
  "cat2_lung"=rhc$cat2=="Lung Cancer",
  "cat2_cirr"=rhc$cat2=="Cirrhosis",
  #"cat2_colon"=rhc$cat2=="Colon Cancer", 
  # Categories of admission diagnosis
  "resp"=rhc$resp=="Yes",
  "card"=rhc$card=="Yes",
  "neuro"=rhc$neuro=="Yes",
  "gastr"=rhc$gastr=="Yes",
  "renal"=rhc$renal=="Yes",
  "meta"=rhc$meta=="Yes",
  "hema"=rhc$hema=="Yes",
  "seps"=rhc$seps=="Yes",
  "trauma"=rhc$trauma=="Yes",
  #"ortho"=rhc$ortho=="Yes",
  rhc$das2d3pc,
  "dnr1"=rhc$dnr1=="Yes",
  "ca_yes"=rhc$ca=="Yes",
  "ca_meta"=rhc$ca=="Metastatic",
  rhc$surv2md1,
  rhc$aps1,
  rhc$scoma1,
  rhc$wtkilo1,
  rhc$temp1,
  rhc$meanbp1,
  rhc$resp1,
  rhc$hrt1,
  rhc$pafi1,
  rhc$paco21,
  rhc$ph1,
  rhc$wblc1,
  rhc$hema1,
  rhc$sod1,
  rhc$pot1,
  rhc$crea1,
  rhc$bili1,
  rhc$alb1,
  #rhc$urin1, # not included in Hirano and Imbens
  # Categories of comorbidities illness:
  rhc$cardiohx,
  rhc$chfhx,
  rhc$dementhx,
  rhc$psychhx,
  rhc$chrpulhx,
  rhc$renalhx,
  rhc$liverhx,
  rhc$gibledhx,
  rhc$malighx,
  rhc$immunhx,
  rhc$transhx,
  rhc$amihx
  #"wt0"=rhc$wtkilo1==0
)
names(rhc.dt) <- sapply(names(rhc.dt), function(x) {
  ifelse(grepl("rhc.",x),strsplit(x,"rhc.")[[1]][2], x)
})




