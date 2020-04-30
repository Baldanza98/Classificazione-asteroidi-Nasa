#---------------------------#
##########LIBRARY############
#---------------------------#
library(ggplot2)
theme_set(theme_classic()+
            theme(axis.line = element_line(colour = "white"))) #tema per i ggplot
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(ggcorrplot)
library(factoextra)
library(mclust)
library(GGally)
library(Rmixmod)
library(stats)
library(labstatR)
library(caret)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)} #funzione che serve per estrapolare la leggenda da un ggplot
#-----------------------------------#
########CARICAMENTO DATASET##########
#-----------------------------------#

#lettura del file
nasa <- read.csv("C:/Users/matte/Desktop/nasa.csv")
dim(nasa)
nasa=as_tibble(nasa)
colnames(nasa)
str(nasa)

#----------------------------------#
###########PRIME ANALISI############
#----------------------------------#


#Avendo troppe variabili per prima cosa devo togliere quelle che non mi servonono
#Una prima analisi semplice e veloce ? direttamente leggendole da str(nasa)
str(nasa)
#alcuni variabili ininfluenti all'analisi
colnames(nasa)
nasa=nasa %>%
  select(everything(),-c("Neo.Reference.ID","Name","Close.Approach.Date","Orbit.Determination.Date","Orbiting.Body",
                           "Equinox"))
#nessun missing value
sum(is.na(nasa))

#suddivisione labels per capire come sono suddivise le classi
table(nasa$Hazardous) #classi sbilanciat!!

#tolgo le etichette
nasa_data=nasa[-length(nasa)]

correlazioni=cor(nasa_data)#correlazioni
ggcorrplot(correlazioni,colors = c("#6D9EC1", "white","#E46726"))+
  theme(axis.text.y = element_text(size = 9))+
  labs(title = "Correlazioni")+
  scale_x_discrete(labels=NULL) #visualizzazion grafica correlazion
  
ggcorrplot(correlazioni,colors = c("#6D9EC1", "white", "#E46726"), hc.order = T, type="lower", lab=F)+
  theme(axis.text.y = element_text(size = 9))+
  labs(title = "Correlazioni")+
  scale_x_discrete(labels=NULL) #visualizzazione grafica correlazione in altro mod
#bisognera sicuramente eliminare alcune variabili,quelle troppo correlate fra di loro

#-----------------------------------------------#
##############FEATURING SELECTION################
#-----------------------------------------------#

#dalla correlazione mi accorgo che alcune variabili sono rindonanti
#spiegano infatti le stesse cose
colnames(nasa_data[,3:9])
nasa_data=nasa_data[,-c(4:9)]
#lo stesso per le altre
colnames(nasa_data[,8:11])
#via 8,9 e 11
nasa_data=nasa_data[,-c(8,9,11)]
colnames(nasa_data)
colnames(nasa_data[,c(5:7)])
#via 6,7
nasa_data=nasa_data[,-c(6,7)]
colnames(nasa_data)

#analizziamo le nuove correlazioni
correlazioni <- cor(nasa_data)

ggcorrplot(correlazioni,colors = c("#6D9EC1", "white", "#E46726"), hc.order = T, type="lower", lab=F)+
  theme(axis.text.y = element_text(size = 9))+
  labs(title = "Correlazioni")+
  scale_x_discrete(labels=NULL)

#le alte correlazioni sono state eliminate,preferisco non eliminare quelle rimaste
colnames(nasa_data)
#noto che posso togliere ancora qualche variabile 
colnames(nasa_data[,2:3])
cor(nasa$Est.Dia.in.KM.max.,nasa$Est.Dia.in.KM.min.)#max correlazione
cor(nasa$Est.Dia.in.KM.max.,nasa_data)#correlata con magnitudine(ha senso non considerarla)
nasa_data=nasa_data[,-c(2,3)]#le tolgo sono spiegate gi? da un altra variabile
length(nasa_data)
colnames(nasa_data)



#--------------------------------------------------#
########SELEZIONI DELLE VARIABILI TRAMITE PCA######
#-------------------------------------------------#
#pca
pca=prcomp(nasa_data,scale. = T)
summary(pca)

c=mean(pca$sdev^2)
pca$sdev^2>c
#scelgo le prime sei componenti

#costruisco la matrice dei punteggi
nasa_final=pca$x[,1:6]

##cerco una figura dove si vedano cluster
res.var <- get_pca_var(pca)  

#contributo alla pca 
for(i in 1:6){
  print(which.max(res.var$contrib[,i]))
}
col=colnames(nasa_data)[c(11,1,18,15,5,16)]

scelta_pca=(nasa_data[,col])
pairs(scelta_pca,col=nasa$Hazardous,lower.panel = NULL)
#nessun cluster evidente!

#---------------------------------#
#########SELEZIONE VARIABILI TRAMITE GLM############
#--------------------------------------------------#
val=ifelse(nasa$Hazardous == "False", 0,1)
#modello di regressione binomiale per la scelta delle variabili
mod=glm(val~.,nasa_data,family=binomial)
summary(mod)
final=step(mod)# stepwise backward per la scelta delle variabili piu significative
summary(final)
#colonne scelte
colonne=c("Absolute.Magnitude" ,"Orbit.ID","Orbit.Uncertainity",
          "Minimum.Orbit.Intersection" ,"Jupiter.Tisserand.Invariant" ,
          "Semi.Major.Axis" ,"Inclination" , "Perihelion.Distance","Perihelion.Time" ,"Mean.Motion")

#riduzione del dataset a colonne scelte
nasa_reduction=nasa_data%>%
  select(all_of(colonne))

#cerco dei cluster evidenti tra le variabili selezionate
pairs(nasa_reduction,col=nasa$Hazardous)

#noto possibile cluster visivo tra le prime osservazioni
Asteroide_Pericoloso=nasa$Hazardous
#ggpairs per vedere come vanno le comparation fra le variabili,e si notato cluster
ggpairs(legend=1,nasa_reduction[1:4],aes(col=Asteroide_Pericoloso),
        lower = list(continuous = wrap("points",alpha = 0.4, size=2)))+
  theme(legend.position = "bottom")

#pairs per vedere tutto
pairs(nasa_reduction,col=Asteroide_Pericoloso)

#creazione dei grafici migliori
gg1=ggplot(nasa_reduction,aes(y=Absolute.Magnitude,x=Minimum.Orbit.Intersection))+
  geom_point(aes(fill=Asteroide_Pericoloso),shape=22,alpha=1.00,colour="black")+
  scale_fill_manual(values = c("#FC4E07","chartreuse3"))

mylegend<-g_legend(gg1)

gg2=ggplot(nasa_reduction,aes(y=Inclination,x=Minimum.Orbit.Intersection))+
  geom_point(aes(fill=Asteroide_Pericoloso),shape=21,alpha=1.00,colour="black")+
  scale_fill_manual(values = c("#FC4E07","chartreuse3"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

gg3=ggplot(nasa_reduction,aes(y=Orbit.ID,x=Minimum.Orbit.Intersection))+
  geom_point(aes(fill=Asteroide_Pericoloso),shape=24,alpha=1.00,colour="grey11")+
  scale_fill_manual(values = c("#FC4E07","chartreuse3"))

#da notare che si vedono i cluster. Ho unito il tutto mantendendo la stessa variabile sulle x
#per avere una miglior visuale grafica e ho mantenuto la stessa legenda per tutti e tre
#i grafici
ggarrange(gg1, gg3,gg2,ncol=1, nrow=3, common.legend = TRUE,legend="top")

#-------------------------------------------------#
########VALUTAZIONE BOXPLOT#########
#----------------------------------#
#valutazione grafica dispersione variabili
colnames(nasa_reduction)
x=nasa_reduction %>%
  select(-Perihelion.Time)
ggplot(stack(x), aes(x = ind, y = values))+
  geom_boxplot(outlier.alpha = 0.0)+
  geom_jitter(shape=21, position=position_jitter(0.2),alpha=0.1,aes(fill=ind))+
  theme(legend.position = "NULL")+
  ylim(0,50)

#------------------------------#
########CLASSIFICATION EDDA##########
#------------------------------------#
#nasa_reduction
(n<-nrow(nasa_reduction))
set.seed(211)
test.set.labels<-sample(1:n,0.20*n)#test set
length(test.set.labels)
set.seed(211)#modello di classificazione provato sul training
pr<-mixmodLearn(scelta_pca[-test.set.labels,], nasa$Hazardous[-test.set.labels] , 
                models=mixmodGaussianModel(family="all",equal.proportions = FALSE),
                criterion=c('CV'))

pr
colnames(nasa_reduction)
pr@bestResult
pr@results[[1]]@model #piu performante
pr@results[[2]]@model #secondo piu performante
pr@results[[1]]@criterionValue [1] 
pr@results[[2]]@criterionValue [1]  

#creazione di un ciclo che mi permetta di capire quale sia il metodo migliore
#plotto graficamente tramite ggplot Cv piu basso raggiunto per modello
CV=NULL
nome=NULL
for(j in 1:100){
  pr<-mixmodLearn(scelta_pca[-test.set.labels,], nasa$Hazardous[-test.set.labels] , 
                  models=mixmodGaussianModel(family="all",equal.proportions=FALSE),
                  criterion=c('CV'))
  for (i in 1: length(pr@models@listModels)){
    ind = which(pr@results [[i]] @model == pr@models@listModels)
    CV =c(CV,pr@results [[i]] @criterionValue [1])
    nome=c(nome,(pr@models@listModels)[ind])
  }}
frame=data.frame(CV,nomi=nome)
frame$nomi=as.character(frame$nomi)
frame=as_tibble(frame)
frame %>%
  group_by(nomi) %>%
  summarize(CV=min(CV)) %>%
  ggplot()+
  geom_point(aes(x=CV,y=nomi,col=nomi),shape=21)+
  theme(legend.position = "none") +
  geom_segment(aes(x=0.150,xend=CV,y=reorder(nomi, CV),yend=reorder(nomi, CV),col=nomi))+
  xlim(0.150,0.190)+
  labs(y="Nome modello")
#miglior modello Gaussian_pk_L_I

train= mixmodLearn(scelta_pca[-test.set.labels,], nasa$Hazardous[-test.set.labels] , 
                   models=mixmodGaussianModel(listModels="Gaussian_pk_L_I",equal.proportions=FALSE),
                   criterion=c('CV'))
?mixmodGaussianModel
summary(train)

pred<- mixmodPredict(scelta_pca[test.set.labels,], classificationRule=train["bestResult"])
pred

mean(as.integer(nasa$Hazardous[test.set.labels]) == pred["partition"])

table(nasa$Hazardous[test.set.labels])
770/(770+167)
#ATTENZIONE! Il modello migliore ? quello che mi mette tutte le classi all'interno delle stessa
#SICURAMENTE PROBLEMA DI SBILANCIAMENTO NEL TRAINING SET!

valori=ifelse(pred["partition"]==1,"False","True")
valori=factor(valori)#errore nessuna classe classificata correttamente
confusionMatrix(nasa$Hazardous[test.set.labels], valori)
#-----------------------------------#
########CLASSIFICATION EDDA-GLM##########
#------------------------------------#

#nasa_reduction,stessi passaggi precedenti
(n<-nrow(nasa_reduction))
set.seed(211)
test.set.labels<-sample(1:n,0.20*n)
length(test.set.labels)
set.seed(211)
pr<-mixmodLearn(nasa_reduction[-test.set.labels,], nasa$Hazardous[-test.set.labels] , 
                models=mixmodGaussianModel(family="all",equal.proportions=FALSE),
                criterion=c('CV'))

colnames(nasa_reduction)
pr@bestResult

pr@results[[1]]@model #piu performante
pr@results[[2]]@model #secondo piu performante
pr@results[[1]]@criterionValue [1] 
pr@results[[2]]@criterionValue [1]  

#creazione di un ciclo che mi permetta di capire quale sia il metodo migliore
CV=NULL
nome=NULL
for(j in 1:100){
  pr<-mixmodLearn(nasa_reduction[-test.set.labels,], nasa$Hazardous[-test.set.labels] , 
                  models=mixmodGaussianModel(family="all",equal.proportions=FALSE),
                  criterion=c('CV'))
  for (i in 1: length(pr@models@listModels)){
    ind = which(pr@results [[i]] @model == pr@models@listModels)
    CV =c(CV,pr@results [[i]] @criterionValue [1])
    nome=c(nome,(pr@models@listModels)[ind])
  }}
frame=data.frame(CV,nomi=nome)
frame$nomi=as.character(frame$nomi)
frame=as_tibble(frame)
frame %>%
  group_by(nomi) %>%
  summarize(CV=min(CV)) %>%
  ggplot()+
  geom_point(aes(x=CV,y=nomi,col=nomi),shape=21)+
  theme(legend.position = "none") +
  geom_segment(aes(x=0.035,xend=CV,y=reorder(nomi, CV),yend=reorder(nomi, CV),col=nomi))+
  labs(y="Nome modello")

#miglior modello Gaussian_pk_Lk_Bk
train= mixmodLearn(nasa_reduction[-test.set.labels,], nasa$Hazardous[-test.set.labels] , 
                   models=mixmodGaussianModel(listModels="Gaussian_pk_Lk_Ck",equal.proportions=FALSE),
                   criterion=c('CV'))
summary(train)


pred<- mixmodPredict(nasa_reduction[test.set.labels,], classificationRule=train["bestResult"])
pred

mean(as.integer(nasa$Hazardous[test.set.labels]) == pred["partition"])

valori=ifelse(pred["partition"]==1,"False","True")
valori=factor(valori)#stavolta riesco a stampare la classificazione
confusionMatrix(nasa$Hazardous[test.set.labels], valori)

#----------------------------------#
########CLASSIFICATION MDA-GLM##########
#----------------------------------#
set.seed(123)
#modello performance prova
mod= MclustDA(nasa_reduction[-test.set.labels ,], nasa$Hazardous[-test.set.labels])
summary(mod)
# library(gtools)
# valori=permutations(5, 2, set=TRUE, repeats.allowed=T)
# bic=NULL
# for(i in 1:25){
#   mod= MclustDA(nasa_reduction[-test.set.labels ,], nasa$Hazardous[-test.set.labels],G=c(valori[i,]))
#   bic=c(bic,mod$bic)
# }
# comb=NULL
# for(i in 1:25){
#   comb=c(comb,paste(as.character(valori[i,1]),as.character(valori[i,2])))
# }
# 
# data=data.frame(bic=bic,modello=comb)
# data %>%
#   filter(modello != c("1 1")) %>%
#   ggplot(aes(y=modello,x=bic))+
#   geom_point(shape=6,aes(col=modello))+
#   geom_segment(aes(x=-30000,xend=bic,y=reorder(modello,bic),yend=reorder(modello,bic),col=modello))+
#   labs(y="K per ogni componente",x="BIC")+
#   xlim(-30000,-10000)+
#   theme(legend.position = "none")
#migliore 1 5 componenti
#verifica modelli

#con due gruppi ? migliore
summary(mod)

sum(predict(mod, nasa_reduction[test.set.labels ,])$class != nasa$Hazardous[test.set.labels])
class=predict(mod, nasa_reduction[test.set.labels ,])$class
confusionMatrix(factor(class),nasa$Hazardous[test.set.labels])


##CLASSIFICATION EDDA-SOLO3 VAR######
colnames(nasa_reduction)
#si prova con solo 3 delle variabili
num=c(1,4,7)
train= mixmodLearn(nasa_reduction[-test.set.labels,num], nasa$Hazardous[-test.set.labels] , 
                   models=mixmodGaussianModel(listModels="Gaussian_pk_Lk_Ck",equal.proportions=FALSE),
                   criterion=c('CV'))

CV=NULL
nome=NULL
for(j in 1:20){
  pr<-mixmodLearn(nasa_reduction[-test.set.labels,num], nasa$Hazardous[-test.set.labels] , 
                  models=mixmodGaussianModel(family="all",equal.proportions=FALSE),
                  criterion=c('CV'))
  for (i in 1: length(pr@models@listModels)){
    ind = which(pr@results [[i]] @model == pr@models@listModels)
    CV =c(CV,pr@results [[i]] @criterionValue [1])
    nome=c(nome,(pr@models@listModels)[ind])
  }}
frame=data.frame(CV,nomi=nome)
frame$nomi=as.character(frame$nomi)
frame=as_tibble(frame)
frame %>%
  group_by(nomi) %>%
  summarize(CV=min(CV)) %>%
  ggplot()+
  geom_point(aes(x=CV,y=nomi,col=nomi),shape=21)+
  theme(legend.position = "none") +
  geom_segment(aes(x=0.150,xend=CV,y=reorder(nomi, CV),yend=reorder(nomi, CV),col=nomi))+
  xlim(0.150,0.190)+
  labs(y="Nome modello")


summary(train)

pred<- mixmodPredict(nasa_reduction[test.set.labels,num], classificationRule=train["bestResult"])
pred

mean(as.integer(nasa$Hazardous[test.set.labels]) == pred["partition"])

valori=ifelse(pred["partition"]==1,"False","True")
valori=factor(valori)
confusionMatrix(nasa$Hazardous[test.set.labels], valori)


#prova diversa con colonne migliori-MDA###############
mod= MclustDA(nasa_reduction[-test.set.labels ,num], nasa$Hazardous[-test.set.labels])
summary(mod)

sum(predict(mod, nasa_reduction[test.set.labels ,num])$class != nasa$Hazardous[test.set.labels])
class=predict(mod, nasa_reduction[test.set.labels ,num])$class
confusionMatrix(factor(class),nasa$Hazardous[test.set.labels])
