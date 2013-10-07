fpol1 <-
function(dets,vec,vars,nom){
                                       # dets : la table vec : vecteur de variable a fixe,vars : 
                                       # variable a sommer, nom : nom de la nouvelle var
                                       # vagp : variable a garder
                                       # vars : variable a sommer 
                                       # nom : nom de la variable 
                                       n<-length(vec)
                                       ux=0
                                       for(i in 1:n){xs=max(dets[vec[i]])
                                                     if(xs!=0){ux<-ux+t(t(dets[,vec[i]]/xs))*10^{n+1-i}
                                                              }else{ux<-ux}
                                                     }
                                       ux<-round(ux,digits=4)
                                       uxf<-duplicated(ux)
                                       dats1<-data.frame(dets[,vec],ux,uxf)
                                       names(dats1)<-c(vec,"ux","uxf")
                                       ndats1<-dats1[dats1["uxf"]==FALSE,]
                                       ndets1<-ndats1[,c(vec,"ux")]                          
                                       #nvt<-dets1[duplicated(dets1)==FALSE,] 
                                       names(ndets1)<-c(vec,"ux")              
                                       vs<-tapply(dets[,vars],ux,sum)
                                       id<-t(t(as.numeric(names(vs))))
                                       vs1<-t(t(vs));rownames(vs1)<-NULL
                                       tab<-data.frame(id,vs1)
                                       names(tab)<-c("ux",nom)
                                       tab2<-merge(ndets1,tab,by=c("ux","ux"))
                                       tab2$ux<-NULL
                                       return(tab2)
                                      }
