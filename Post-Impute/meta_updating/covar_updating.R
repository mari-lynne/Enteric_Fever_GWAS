#Update Covar files

#Needs to be recoded as categorical I think 
covar <- fread("all_covar.txt")

#Recode Sex
#Recode Challene
#Recode Vaccine
#Recode rechallene
#Recode (and mutate) rechallenge type (hetero or homo)

#NONE = 0
#Ty21a = 1
#M01ZH09= 2
#Vi-PS = 3
#Vi-TT= 4

covar$Vaccine <- ifelse(covar$Vaccine == "None",0,
                            ifelse(covar$Vaccine == "Ty21a",2,
                                   ifelse(covar$Vaccine == "M01ZH09", 2,
                                          ifelse(covar$Vaccine == "Vi-PS", 3,
                                                 ifelse(covar$Vaccine == "Vi-TT", 4,NA)))))

#Recode strain
covar$Challenge <- as.factor(covar$Challenge)
levels(covar$Challenge)
covar <- filter(covar, Challenge != "C")
#removed 3 rows (bicarb people probablly)

#C, P = 3, T = 1, TN = 2
levels(covar$Challenge) <- c("0", "3", "1", "2")

#sex M=1, F=2
covar$Sex <- ifelse(covar$Sex == "M", 1, ifelse(covar$Sex == "F", 2, NA))

#Re-challenge
covar$Re_Challenge <- as.factor(covar$Re_Challenge)
levels(covar$Re_Challenge)
#N=1, PP=2, PT=3, TP=4, TT=5

levels(covar$Re_Challenge) <- c("1", "2", "3", "4", "5")

#Key: ####

#Strain: C, P = 3, T = 1, TN = 2
#Sex M=1, F=2
#Re-challenge N=1, PP=2, PT=3, TP=4, TT=5
#Vaccine: NONE = 0, Ty21a = 1, M01ZH09= 2, Vi-PS = 3, Vi-TT= 4
