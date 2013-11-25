Spmlficmcm <-
function(fl,N,gmname,gcname,DatfE,typ,start){ 
                  # ============================================================ 
                  # Etape 1 estimation des valeurs initiales 
                  # Valeurs initiales des parametre du modele 
                  # baterie des tests
                    genom<-DatfE[gmname];genom1<-genom[is.na(genom)!=TRUE]
                    genoen<-DatfE[gcname];genoen1<-genoen[is.na(genoen)!=TRUE]
                    nagm<-genom[is.na(genom)==TRUE]
                    dat<-DatfE[is.na(DatfE[gcname])!=TRUE,]
                    gt1<-dat[gmname];
                    gt2<-dat[gcname];
                    teg1<-ifelse(gt1==0 & gt2==2,1,0)
                    teg2<-ifelse(gt1==2 & gt2==0,1,0)
                    tegk<-teg1+teg2
                    # conditions
                    if(length(nagm)>0){print(gmname)
                                      stop("missing mother genotype value")} 
                    if(min(genom1)<0){print(gmname)
                                 stop("mother's genotype is negative")}
                    if(max(genom1)>2){print(gmname)
                                 stop("mother's genotype is greater than 2")}
                    if(min(genoen1)<0){print(gmname)
                                 stop("child's genotype is negative")}
                    if(max(genoen1)>2){print(gmname)
                                 stop("child's genotype is greater than 2")} 
                    if(max(tegk)>0){print(gmname)
                               stop("mother and child genotypes are not compatible")
                               }else{
                               if(typ==1){
                               vIn<-Est.Inpar(fl,N,gmname,gcname,DatfE,1)
                               }else{
                               vIn<-Est.Inpar(fl,N,gmname,gcname,DatfE,2)
                                }   
                    if(missing(start)){parms<-vIn$parms
                                             }else{parms<-start}
                    beta.start=parms[1:length(parms)-1];
                    theta.start=parms[length(parms)]
                  # Valeurs initiales du systeme d quation non linaire
                  ma.u<-vIn$ma.u
                  vecma.u=c(ma.u[1,],ma.u[2,]) 
                  #=============================================================
                  # Etape 2 resolution du systeme d equation  
                  RSeq<-Nlsysteq(fl,DatfE,N,gmname,gcname,beta.start,theta.start)
                  SS<-nleqslv(vecma.u,RSeq) 
                  vecma.u<-SS$x
                  #=============================================================
                  # Etape 3 ecriture de la vraisemblance profile
                  if(typ==1){
                    ftlh<-ft_likhoodCas1(fl,DatfE,N,gmname,gcname,vecma.u)
                            }else{
                    ftlh<-ft_likhoodCasM(fl,DatfE,N,gmname,gcname,vecma.u)
                                  }        
                  #=============================================================
                  # Etape 4 calcul des estimateurs  
                  # calcul du gradian
                  Grad1<-grad(ftlh,parms)
                  #calcul de la hessian
                  Hes<-hessian(ftlh,parms);
                  # estimateur des parametre
                  Parms.est=parms-solve(Hes)%*%Grad1
                  # calcul de la variance 
                  Hes1<-hessian(ftlh,Parms.est);
                  matv<-(-1)*solve(Hes1)
                  var.par<-sqrt(diag(matv))
                  #=============================================================
                  # Etape 5 preparation des resultats
                  mats<-cbind(Parms.est,var.par)
                  nma<-c("Intercept",attr(terms.formula(fl),"term.labels"),"theta")
                  nac<-c("Estimate","Std.Error")
                  colnames(mats)<-nac
                  rownames(mats)<-nma
                  loglik<-ftlh(mats[,"Estimate"])
                  rr<-list(Uim=vecma.u,MatR=mats,Matv=matv,Lhft=ftlh,Value_loglikh=loglik)
                  return(rr)
                  }
                  }
