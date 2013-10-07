fctcd <-
function(data1,gcname){
                               # gcname : est le genotype de l enfant
                               datnmv<-data1[is.na(data1[gcname])!=TRUE,] # donnees non manquantes
                               datmv<-data1[is.na(data1[gcname])==TRUE,] # donnees manquante
                               vec<-c(0,1,2)
                               if(dim(datmv)[1]!=0){
                                                   # construction du genotype de l enfant sur les donnees manquantes 
                                                     gc<-rep(vec,dim(datmv)[1])
                                                     datdm<-datmv[rep(1:nrow(datmv),rep(3,nrow(datmv))),]
                                                     datdm[gcname]<-gc; # table des donnes complete
                                                   }else{datdm<-NULL}
                               # contruction du nouveau genotype de l enfant sur les donnees complet
                               datnmv0<-datnmv
                               gc<-rep(vec,dim(datnmv0)[1])
                               datdmcp<-datnmv0[rep(1:nrow(datnmv0),rep(3,nrow(datnmv0))),]
                               datdmcp$gnpo<-gc
                               datdmcp$vdcop<-as.numeric(datdmcp[gcname]==datdmcp["gnpo"])
                               datdmcp[gcname]<-gc;datdmcp$gnpo<-NULL
                               return(list(datnmv=datnmv,datdm=datdm,datmv=datmv,datdmcp=datdmcp))
                               }
