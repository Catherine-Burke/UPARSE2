#Plot scatterplot of OTU3 abundance in wound swabs against wound size
woundswab_sd<-rownames(sample_data(rare_mydata_sub)[sample_data(rare_mydata_sub)[,"sampletype2"]=="wound_swab"])
woundswab_OTU_3<-data.frame(otu_table(rare_mydata_sub)["OTU_3",woundswab_sd])
woundswab_OTU_3
woundswab_OTU_3<-woundswab_OTU_3[,order(names(woundswab_OTU_3))]
woundswab_OTU_3<-(woundswab_OTU_3/30000)*100
woundswab_OTU_3<-as.matrix(woundswab_OTU_3)
woundarea<-c(110,64,24,24,14,18,35,32,33,36,16,4,2,2,6,6,6,48,4,15,25,4,5,4,3,0,72,66,18,0,16,48,35,40,40,200,220,220,171)
woundswab_OTU_3<-rbind(woundswab_OTU_3,woundarea)
woundswab_OTU_3<-t(woundswab_OTU_3)
row.names(woundswab_OTU_3)<-NULL
OTU3<-woundswab_OTU_3[,"OTU_3"]
area<-woundswab_OTU_3[,"woundarea"]
plot(OTU3, area, xlab="OTU3 relative abundance", ylab="wound area (mm2)", main="OTU3 abundance in wounds vs wound size\n for all Patients")
abline(lm(area~OTU3), col="red")
#OTU_8652 abundance in wound swabs vs wound size, all patients
woundswab_OTU_8652<-data.frame(otu_table(rare_mydata_sub)["OTU_8652",woundswab_sd])
woundswab_OTU_8652
woundswab_OTU_8652<-woundswab_OTU_8652[,order(names(woundswab_OTU_8652))]
woundswab_OTU_8652<-(woundswab_OTU_8652/30000)*100
woundswab_OTU_8652<-as.matrix(woundswab_OTU_8652)
woundarea<-c(110,64,24,24,14,18,35,32,33,36,16,4,2,2,6,6,6,48,4,15,25,4,5,4,3,0,72,66,18,0,16,48,35,40,40,200,220,220,171)
woundswab_OTU_8652<-rbind(woundswab_OTU_8652,woundarea)
woundswab_OTU_8652<-t(woundswab_OTU_8652)
row.names(woundswab_OTU_8652)<-NULL
OTU8652<-woundswab_OTU_8652[,"OTU_8652"]
area<-woundswab_OTU_8652[,"woundarea"]
plot(OTU8652, area, xlab="OTU8652 relative abundance", ylab="wound area (mm2)", main="OTU8652 abundance in wounds vs wound size\n for all Patients")
abline(lm(area~OTU8652), col="red")
#Plot for OTU3 abundance in skin against wound size
diab_skin_names<-rownames(sample_data(rare_mydata_sub)[sample_data(rare_mydata_sub)[,"sampletype2"]=="diabetic_skin_contra"])
diabskin_OTU_3<-data.frame(otu_table(rare_mydata_sub)["OTU_3",diab_skin_names])
diabskin_OTU_3<-diabskin_OTU_3[,order(names(diabskin_OTU_3))]
diabskin_OTU_3<-(diabskin_OTU_3/30000)*100
diabskin_OTU_3<-as.matrix(diabskin_OTU_3)
woundarea<-c(110,64,24,24,14,18,35,32,33,36,16,4,2,6,6,48,4,15,25,4,5,4,3,0,0,72,66,18,0,0,16,48,35,40,40,200,200,220,220,220,171)
diabskin_OTU_3<-rbind(diabskin_OTU_3,woundarea)
diabskin_OTU_3<-t(diabskin_OTU_3)
row.names(diabskin_OTU_3)<-NULL
OTU3<-diabskin_OTU_3[,"OTU_3"]
area<-diabskin_OTU_3[,"woundarea"]
plot(OTU3, area, xlab="OTU3 relative abundance", ylab="wound area (mm2)", main="OTU3 abundance in skin\n vs wound size for all Patients")
abline(lm(area~OTU3), col="red")
#OTU1444
diabskin_OTU_1444<-data.frame(otu_table(rare_mydata_sub)["OTU_1444",diab_skin_names])
diabskin_OTU_1444<-diabskin_OTU_1444[,order(names(diabskin_OTU_1444))]
diabskin_OTU_1444<-(diabskin_OTU_1444/30000)*100
diabskin_OTU_1444<-as.matrix(diabskin_OTU_1444)
woundarea<-c(110,64,24,24,14,18,35,32,33,36,16,4,2,6,6,48,4,15,25,4,5,4,3,0,0,72,66,18,0,0,16,48,35,40,40,200,200,220,220,220,171)
diabskin_OTU_1444<-rbind(diabskin_OTU_1444,woundarea)
diabskin_OTU_1444<-t(diabskin_OTU_1444)
row.names(diabskin_OTU_1444)<-NULL
OTU1444<-diabskin_OTU_1444[,"OTU_1444"]
area<-diabskin_OTU_1444[,"woundarea"]
plot(OTU1444, area, xlab="OTU1444 relative abundance", ylab="wound area (mm2)", main="OTU1444 abundance in skin\n vs wound size for all Patients")
abline(lm(area~OTU1444), col="red")
#OTU15 in skin vs wound size for all patients
diabskin_OTU_15<-data.frame(otu_table(rare_mydata_sub)["OTU_15",diab_skin_names])
diabskin_OTU_15<-diabskin_OTU_15[,order(names(diabskin_OTU_15))]
diabskin_OTU_15<-(diabskin_OTU_15/30000)*100
diabskin_OTU_15<-as.matrix(diabskin_OTU_15)
woundarea<-c(110,64,24,24,14,18,35,32,33,36,16,4,2,6,6,48,4,15,25,4,5,4,3,0,0,72,66,18,0,0,16,48,35,40,40,200,200,220,220,220,171)
diabskin_OTU_15<-rbind(diabskin_OTU_15,woundarea)
diabskin_OTU_15<-t(diabskin_OTU_15)
row.names(diabskin_OTU_15)<-NULL
OTU15<-diabskin_OTU_15[,"OTU_15"]
area<-diabskin_OTU_15[,"woundarea"]
plot(OTU15, area, xlab="OTU15 relative abundance", ylab="wound area (mm2)", main="OTU15 abundance in skin\n vs wound size for all Patients")
abline(lm(area~OTU15), col="red")
#OTU_507 in skin vs wounds size for all patients
diabskin_OTU_507<-data.frame(otu_table(rare_mydata_sub)["OTU_507",diab_skin_names])
diabskin_OTU_507<-diabskin_OTU_507[,order(names(diabskin_OTU_507))]
diabskin_OTU_507<-(diabskin_OTU_507/30000)*100
diabskin_OTU_507<-as.matrix(diabskin_OTU_507)
woundarea<-c(110,64,24,24,14,18,35,32,33,36,16,4,2,6,6,48,4,15,25,4,5,4,3,0,0,72,66,18,0,0,16,48,35,40,40,200,200,220,220,220,171)
diabskin_OTU_507<-rbind(diabskin_OTU_507,woundarea)
diabskin_OTU_507<-t(diabskin_OTU_507)
row.names(diabskin_OTU_507)<-NULL
OTU507<-diabskin_OTU_507[,"OTU_507"]
area<-diabskin_OTU_507[,"woundarea"]
plot(OTU507, area, xlab="OTU507 relative abundance", ylab="wound area (mm2)", main="OTU507 abundance in skin\n vs wound size for all Patients")
abline(lm(area~OTU507), col="red")
#OTU_406
diabskin_OTU_406<-data.frame(otu_table(rare_mydata_sub)["OTU_406",diab_skin_names])
diabskin_OTU_406<-diabskin_OTU_406[,order(names(diabskin_OTU_406))]
diabskin_OTU_406<-(diabskin_OTU_406/30000)*100
diabskin_OTU_406<-as.matrix(diabskin_OTU_406)
woundarea<-c(110,64,24,24,14,18,35,32,33,36,16,4,2,6,6,48,4,15,25,4,5,4,3,0,0,72,66,18,0,0,16,48,35,40,40,200,200,220,220,220,171)
diabskin_OTU_406<-rbind(diabskin_OTU_406,woundarea)
diabskin_OTU_406<-t(diabskin_OTU_406)
row.names(diabskin_OTU_406)<-NULL
OTU406<-diabskin_OTU_406[,"OTU_406"]
area<-diabskin_OTU_406[,"woundarea"]
plot(OTU406, area, xlab="OTU406 relative abundance", ylab="wound area (mm2)", main="OTU406 abundance in skin\n vs wound size for all Patients")
abline(lm(area~OTU406), col="red")



#P1 wound samples vs OTU62 Pseudomonas
P1_woundnames<-sample_data(rare_mydata_sub)[sample_data(rare_mydata_sub)[,"sampletype1"]=="wound"]
P1_woundnames<-row.names(P1_woundnames)[P1_woundnames[,"subject"]=="P1"]
P1wounds_OTU62<-data.frame(otu_table(rare_mydata_sub)["OTU_62",P1_woundnames])
P1wounds_OTU62<-P1wounds_OTU62[,order(names(P1wounds_OTU62))]
P1wounds_OTU62<-(P1wounds_OTU62/30000)*100
P1wounds_OTU62<-as.matrix(P1wounds_OTU62)
P1_woundarea<-c(110,64,24,24,14,18)
P1wounds_OTU62<-rbind(P1wounds_OTU62,P1_woundarea)
P1wounds_OTU62<-t(P1wounds_OTU62)
row.names(P1wounds_OTU62)<-NULL
OTU62<-P1wounds_OTU62[,"OTU_62"]
P1area<-P1wounds_OTU62[,"P1_woundarea"]
plot(OTU62, P1area, xlab="OTU62 relative abundance", ylab="wound area (mm2)", main="P1 wound area vs OTU62 abundance")
abline(lm(P1area~OTU62), col="red")
#now for OTU 233
P1wounds_OTU233<-data.frame(otu_table(rare_mydata_sub)["OTU_233",P1_woundnames])
P1wounds_OTU233<-P1wounds_OTU233[,order(names(P1wounds_OTU233))]
P1wounds_OTU233<-(P1wounds_OTU233/30000)*100
P1wounds_OTU233<-as.matrix(P1wounds_OTU233)
P1_woundarea<-c(110,64,24,24,14,18)
P1wounds_OTU233<-rbind(P1wounds_OTU233,P1_woundarea)
P1wounds_OTU233<-t(P1wounds_OTU233)
row.names(P1wounds_OTU233)<-NULL
OTU233<-P1wounds_OTU233[,"OTU_233"]
P1area<-P1wounds_OTU233[,"P1_woundarea"]
plot(OTU233, P1area, xlab="OTU233 relative abundance", ylab="wound area (mm2)", main="P1 wound area vs OTU233 abundance")
abline(lm(P1area~OTU233), col="red")
model <- lm(P1area~OTU233)
par(mfrow = c(2,2))
plot(model)
