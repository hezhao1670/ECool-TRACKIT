! everything in MKS
      implicit complex (z)
      implicit real(8) (a-h,o-y)
      implicit integer (i-n)
      parameter(nb=600000,pi=3.1415926536,clight = 299792458.0,E00=931.494)
	  parameter(epsilon=8.8541878e-12,e0=1.6021766e-19)
	  parameter(a_me=9.1093837e-31,a_mp=1.6726219e-27,a_mion=1.6605833e-27)
      parameter(nnebeam=150)
      common/inputi/np,nturns,ndim,nwrite,nresamp,nslice,nsliceIBS, nharm,nharm2,mbunch,mode,np0,nresamp2,nextra,ntdelay,idoecool
!23456
      common/rfstuff/gammat,circ,gamma0,vrf,trf,omegarf,tcoeff,pcoeff,frf,phisynch,dqmin
! x2 is horizontal
      common/inject/ampx,ampy,tauhat,dgamma_hat,coastbdpp,vrf2
      common/dynamics/beta,vrev,trev,cpermacro
      common/tran_in/tunex,tuney,chrom
      common/randomstuff/iseed
      common/species/aatom,qatom,power,pnumber,plifetime
      common/particles/y(nb),py(nb),t(nb),pt(nb),x(nb),px(nb)
      common/keepit/xkeep(7,nb)
      common/tunebeta/xtune(nb),dpwake(nb),dpwake2(nb)
      common/resden/tlob,thib,thixx,thixx2,denlon2(nb)
      common/diag/avgchrom,avgtune,csfull,dptrf(nb)
      common/tstuff/kturns,nturnon,nperturn,nwrite0
      common/diag/xpavg(nb),avgline(nb),rmsx,rmsy,dponp
      common/ibstuff/fracibstot,fracibsx,fracecool
      common/lumistuff/peakilumi,fwhmlumi
      common/burnit/cross_section,betastar,dthour,collision(nb)
      common/burnit2/nips,nebeam,idosc,ntloss
      common/ebeamstuff/ebeams,ebeamx,ebeamy,sigurestx,siguresty,sigurestt,qbeamcool,slencool
      common/detune/detunediag,detunediag2,detunecross,ddux,dduz
      common/ebeam2/betaxcool,betaycool,betavgx,betavgy
      common/speedum/speedfact,dispxe,dispxh
      common/cooler/caray1(-200:200),caray2(-200:200),caray3(-200:200)
      common/offsets/xoff,xpoff,uoff
	  
!##################
! new variables for B-M IBS model 	  
      common/fnewhe/dentran2(nb),avglinet(nb),betamax,aperture,DA
      common/ibslatt/bxa(nb),bxap(nb),byap(nb),bya(nb),dxa(nb),dxpa(nb),fraca(nb),dya(nb),dypa(nb)
      common/paribs/pnum,epsx,epsy,gamma,sigs,dponpibs,coulomb,r0,sumx,sumy,sump
      common/iparms/kloop,nsim,iichoice,iechoice
      
      
! new variables for pulsed dc cooling (parkhomchuk)     
      common/pulsedce/tie,aeb,Bg,Berror,pulseW,pulseN,riseW,gapW,vkick,temet,temel
	  common/hollowe2/nedis
	  common/hollowe/hollowsig,hollowaeb
      common/pulsedce2/dcoolx,dcooly,dcoolp,feshift,fetune,ekhi,eklo,dkicky
      common/pulsedce3/fptune,ephi,eplo,epsoff
	  common/tuning/ntchoice
      common/edistri/eldis(300),sldis(300),ssdis(300)
       
! new variables for pulsed ion beam and multi-injection  	  
	   common/ionpulse/pulseiW,gapit,spareW
	   common/ionpulse1/injN,ninjloss
   
