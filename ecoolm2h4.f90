!
!  ifort -o ecoolm2h2 ecoolm2h2.f90 getIBS.f90

		program ecoolm2h3
!-------------- ecoolm2h1 -------
! included detuning with amplitude 
! changed write in ibslong to be at nwrite
! made ibslong into real, not scaled rates
! modifying ibslong to include viscous forces. 1st pass looks reasonable
! added percentiles
! added tune measurements
! added electron cooling
! added two temperatures for electron cooling
! STRIPPING DOWN TRACKIT7
! 1st remove coherent electron force in ebeamkick
! add speedup factor to ibs and ecool
! reinsert burnoff
! put in brat
! uses x(k), px(k), y(k), py(k), t(k), dgamma(k)

!------ 2020-02-20--------ecoolm2h2
! added the B-M IBS model (Nagaitsev's mothod)
! add the file 'getIBS.f90'
        
!------ 2021-12-02--------ecoolm2h3
! rfupdate(fmix) to rfupdate(fmix,ibeam): add coasting beam
! addparticles to addparticles(iichoice): add coasting beam
! use random_number() instead of ran3(seed): better random number generator
! added coasting beam in ibs model: if(iichoice.eq.1) sigs = circ/2./sqrt(pi)
! added pulsed e-beam cooling with Parhomchuk formula (uniform distribution in long. and trans.)
!		including frequency shift and edge voltage
! added ebeam energy tuning (frequency: square-saw-triangle-sin)

!------ 2022-05-04--------ecoolm2h4
! added multi-injection of pulsed ion beam in subroutine reinjection(ibeam) for BB injection simulation.
! this multi-injection only turn on for iichoice=3
! add xoff setting in the magnetic electron cooling 
		
!------ 2022-08-05--------ecoolm2h4
! added various e-beam distribution uniform/elliptical/gaussian/hollw etc.  nedis/hollowsig/hollowaeb

!------ 2022-09-29--------ecoolm2h4
! added xoff,xpoff and uoff in non magnetized and magnetized  cooling 	
! added ion beam lifetime (exponential)
		
		
		
		include 'ecoolm2h4_c.f90'
! getinput and add particles are OK
        call random_seed()	!initialize the random seed
		call getinput
		call addparticles(iichoice)
! zero out tune sums
		call zerokeep
! check the order of magnitude
		call coolcheck
		
		nextra = 10003 x^2 + Log[(Pi - x)^2]/Pi^4 + 1
		do kturns = 0,nturns+nextra
! logic for writes
			kcheck = kturns/nwrite
			kcheck = kcheck*nwrite
			if(kcheck.eq.kturns)then
				iwrite = 1
			else
				iwrite = 0
			endif
! always write the last turn
			if(kturns.eq.nturns+nextra)then
				write(6,*)' writing last turn'
				nwrite=-1
				iwrite=1
			if(nperturn.lt.10)nperturn=10
			endif
			
!
!! write last several turns
!			if(kturns.gt.nturns) then
!				write(*,689) x(333),x(1333),x(2333),pt(333)/gamma0/beta**2,pt(1333)/gamma0/beta**2,pt(2333)/gamma0/beta**2
!			endif				
!689 	format(6e15.6)
			
! main loop
! get densities
! need to get new density subroutines
! two possibilities
! rhic an lhc here
			
			if(iichoice.eq.3) call reinjectbeam(iichoice)
            
			!write(6,*)'1--',np
			
			call zerokeep
			one=1.
			call rfupdate(one,iichoice)
            
			do k=1,nperturn
				fmix=1./nperturn
				call sptranssc(fmix)
				call keepturn
			enddo
			call sptrans2	!skew quad.
!  ibs		
			call ibslong2
			
			
! electron bunches and drifts around cooling section
				
			call etuning(ntchoice)
			
			if((idoecool.eq.1).and.(iechoice.eq.2))call ebeamkick
            
            if((idoecool.eq.1).and.(iechoice.eq.1))call ebeamkick2
            
! write if desired
			call writit(iwrite)
            
			

		enddo
	
		call getrate
		stop
    end

    
    
!234567890123456789012345678901234567890123456789012345678901234567890
    subroutine getinput
		include 'ecoolm2h4_c.f90'
! 
		open(unit=10,file='ecoolm2h4.in',status='unknown')
! number of turns, dimensionality for distribution, number of turns between writes,
! random number seed
		read(10,*)nturns,ndim,nwrite,nperturn, iseed
!		write(6,*)'wr1',nturns,ndim,nwrite,nperturn,iseed

! transition gamma, circumference, fundamental rf voltage, primary
! harmonic number, secondary rf voltage,
! secondary harmonic number
		read(10,*)gammat,circ,gamma0,vrf,nharm,vrf2,nharm2
!		write(6,*)gammat,circ,gamma0,vrf,nharm,vrf2,nharm2
! betatron tunes, rms sizes (2 is horizontal),chrom same for both planes,
! nresamp,nresamp2: bins for longitudinal density and smoothing
		read(10,*)tuney,ampy,tunex,ampx,chrom,nresamp,nresamp2
		rmsx=ampx
		rmsy=ampy
		thixx = 10.*ampx	!physical aperture
!		write(6,*)tuney,ampy,tunex,ampx,chrom,nresamp,nresamp2
! number of particles per bunch, atomic mass and number
! power=bunch shape parameter
! tauhat = half bunch length
! phisynch = synchronous phase in rf radians
! nturnon=number of turns to turn on nonlinear RF and longitudinal wakes
		read(10,*)pnumber,aatom,qatom,power,tauhat,coastbdpp,nturnon
!		write(6,*)pnumber,aatom,qatom,power,tauhat,coastbdpp,nturnon
! iichoice = 1-coasting beam/2-bunched beam 
! iechoice = cooling force mode (1-magnetic cooling;  2-non-magnetic cooling)
		read(10,*)iichoice,iechoice
        if(iichoice.eq.1) write(*,*) 'Coasting ion beam!'
        if(iichoice.eq.2) write(*,*) 'Bunched ion beam!'
		if(iichoice.eq.3) write(*,*) 'Pulsed ion beam!'
        
! 
!23456
! fraction of z and s ibs strength, other tune,other amplitude, tune split, edge of bunch
		read(10,*)fracibstot,fracibsx,fracecool,dqmin,thib
         
! electron bunch length,rms widths at center of cooling section in meters
		read(10,*)ebeams,ebeamx,ebeamy
! flag for space charge, 1=t flag for ecooling,rms vx/c,vy/c,vt/c of beta for electrons in rest frame
! magnitude of charge per bunch for cooling and length of cooling section, speedup factor
		read(10,*)idosc,idoecool,sigurestx,siguresty,sigurestt,qbeamcool,slencool,speedfact
! horizontal dispersion in cooling section for electrons and hadrons, hadron beta functions in IP
! number of cooling kicks
		read(10,*)dispxe,dispxh,betaxcool,betaycool,nslice
! physical aperture and dynamic aperture for particle loss	  
		read(10,*)betamax,DA,plifetime	
		aperture = DA * ampx
		ntloss = 0
! offsets of electron dp/p and horizontal position 
		read(10,*) xoff,xpoff,uoff
    
!============== pulsed i-injection (iichoice=3) ==================
        read(10,*)
		read(10,*) pulseiW,gapit
		read(10,*) spareW
		ninjloss = 0
		injN = 1
		
!=============== pulsed e-beam (iechoice=1) ==================	
! read pulsed e-beam parameters + Parkhomchuk force+uniform e-beam
        read(10,*)
!e-beam current(A), radius(cm), temperature_tran(eV)., temperature_long(eV)
        read(10,*) tie,aeb,temet,temel
!e-beam distribution choice/hollow beam
!nedis=1-6 (uniform/elliptical/parabolic/cos-square/guassian/hollow)
		read(10,*) nedis,hollowsig,hollowaeb
		if(nedis.eq.1) write(*,*) 'Uniform e-beam!'
        if(nedis.eq.2) write(*,*) 'Elliptical e-beam!'
		if(nedis.eq.3) write(*,*) 'Parabolic e-beam!'
		if(nedis.eq.4) write(*,*) 'Cos-square e-beam!'
        if(nedis.eq.5) write(*,*) 'Gaussian e-beam!'
		if(nedis.eq.6) write(*,*) 'Hollow e-beam!'
		write(*,*) '================================================================'
		
		
		
!Bfield(T), Bfield_error, e-beam frequency shift
        read(10,*) Bg,Berror,feshift
!pulse e-beam setting (unit=s)
        read(10,*) pulseW,pulseN,riseW,gapW

		
!=============== ebeam tuning ==================	
!ebeam energy tuning 
		read(10,*)
		read(10,*) ntchoice
!tuning frequency,high energy, low energy, tuning choice 
		read(10,*)fetune,eklo,ekhi
		read(10,*)fptune,eplo,ephi

		close(10)

		
!read the beta function of proton ring	
!		open(unit=10,file='disp.in',status='unknown')
!		read(10,*) dispxe,dispxh
!		write(*,*) "disp_input",dispxe,dispxh
!		close(10)	
		
	  
! read the lattice for the IBS calculation
		open(unit=10,file='optics.dat',status='unknown')
! number of lines of optics file, number of lines to skip after first (15 nominally)
		read(10,*)ntotal,nskip
		!write(*,*) "lattice: ", ntotal,nskip
		do k=1,nskip
			read(10,*)
		enddo
		nsliceIBS = ntotal-1-nskip
		sold=0

		do k=1,nsliceIBS
			read(10,*)s,betx,alfx,bety,alfy,dx,dpx,dy,dpy
			bxa(k)=betx
			! bxap is d beta/ds
			bxap(k)=-2*alfx
			byap(k)=-2*alfy*1.0
			bya(k)=bety
			dxa(k)=dx
			dxpa(k)=dpx
			dya(k)=dy*1.0
			dypa(k)=dpy*1.0
			fraca(k)=s-sold
			sold=s
		end do
		circ = sold		
		close(10)
        write(*,*) 'lattice loaded!'
		
		
! calculate trev, beta etc
		beta= sqrt(1.-1./gamma0**2)
		vrev = 2.998e8*beta
		trev = circ/vrev
		trf = trev/nharm
		frf = 1/trf
		omegarf = 2*pi/trf
		eta = 1/gammat**2 - 1/gamma0**2
		!write(6,*)'eta = ',eta
		if(iichoice.eq.1.OR.iichoice.eq.3) thib = trev	!thib must be trev for coasting beam for densitu calcualtion
        !write(6,*) 'thib = ',thib
		alllength = (pulseW+2.0*riseW)*pulseN + +gapW*(pulseN-1.0)! + pulseiW
		if(iichoice.eq.3.AND.alllength.gt.thib)then
			write(6,*)'wrong time structure of i/e beam (length>thib)!'
			stop
		endif
! longitudinal eom, t(k) = arrival time with synchronous point at pi rf radians
!                  pt(k) = gamma-gamma0
! d pt(k)/dn = e*qatom*vrf*omegarf*(t(k)-trf/2)/(aatom*pmass*clight**2)
! d  t(k)/dn = trev*eta*pt(k)/(beta0**2*gamma0)
! coefficients  d pt(k)/dn = pcoeff*(t(k)-trf/2) , d  t(k)/dn = tcoeff*pt(k) 
		pcoeff = qatom*vrf*omegarf/(E00*1.0e6*aatom)
		tcoeff = trev*eta/(beta**2*gamma0)
		rcoeff = pcoeff/tcoeff
		if(rcoeff.gt.0)then
			write(6,*)'wrong sign of voltage for sign of eta'
			stop
		endif
! given gap volts get (gamma-gamma0) corresponding to amplitude in arrival time variation
		dgamma_hat = tauhat*sqrt(-rcoeff)
		write(6,*)'amp of dp/p (RF) = ',dgamma_hat/gamma0/beta**2
		write(6,*)'max change in tau per turn (RF) = ',dgamma_hat*tcoeff
! synchrotron tune
		tunes = sqrt(-pcoeff*tcoeff)/(2*pi)
		tunesexact = acos(1+pcoeff*tcoeff/2)/(2*pi)
		write(6,*)'synchrotron frequency = ',tunes/trev
		write(6,*)'revolution period = ',trev
		betavgx =  circ/(2*pi*tunex)
		betavgy =  circ/(2*pi*tuney)
		write(6,*)'reference beta functions (x,y) ', circ/(2*pi*tunex),  circ/(2*pi*tuney)
		write(6,*)'rms emittances = ',ampx*ampx*2*pi*tunex/circ, ampy*ampy*2*pi*tuney/circ
		write(6,*)'rms normalized emittances = ',ampx*ampx*beta*gamma0*2*pi*tunex/circ, ampy*ampy*beta*gamma0*2*pi*tuney/circ 
!23456
		open(unit=33,file='csmon.out',status='unknown')
! normalization for longitudinal profiles
		nwrite0 = nwrite
! equivalent time
		write(6,*)'equivalent time (s/hours) = ', speedfact*nturns*trev,speedfact*nturns*trev/3600.
		dthour = speedfact*nwrite*trev/3600
! zero out memory arrays
		do k=1,nresamp
			avgline(k)=0
			avglinet(k) = 0
        enddo
        
        
        glonge=1.+temel/0.511E6
		vlonge=clight*sqrt(1.-1./glonge/glonge)
		Veff=sqrt((beta*gamma0*clight*Berror)**2+vlonge**2)
		write(*,*) 'Veff(-):',Veff,Veff/beta/clight
       
        write(*,*) '================================================================'
        !CALL EBEAM_DISTRIBUTION()
		CALL EBEAM_DISTRIBUTION_NEW(nedis)
        pulseW = pulseW/trev	!now pulseW is the ratio of W_e/circ
		riseW = riseW/trev
		gapW =gapW/trev
		
		aa = 10.0
		bb = aeb
        vkcon= tie/beta/clight/4./pi/epsilon/gamma0**2
		if(riseW.ne.0.0) vkick = vkcon/(riseW*circ)*slencool*(1.0+2.0*log(aa/bb))
		write(6,*) 'vkick_edge (V)=', vkick
		
		write(*,*) '================================================================'

	end
    
    
    
	
!--------------------------------------------------------------------	
!234567890123456789012345678901234567890123456789012345678901234567890
!energy tunning=uoff
!position tunning = epsoff
!
	subroutine etuning(ituning)
	include 'ecoolm2h4_c.f90'
	
		real(8) timek,tunet,ptest
		
	! 1 - square wave
	! 2 - saw wave
	! 3 - triangle wave
	! 4 - sine wave
		
!################################################################
		if(fetune.eq.0.0) goto 121
						
		timek = kturns*speedfact*trev
		tunet = 1.0/fetune
		tunet2 = tunet/2.0
		ptest = mod(timek,tunet)
		
		if(ituning.eq.0.) then
			goto 121
		
		elseif(ituning.eq.1) then
			if(ptest.lt.tunet2) uoff = ekhi
			if(ptest.ge.tunet2) uoff = eklo

		elseif(ituning.eq.2) then 
			uoff = (ekhi-eklo)/tunet*ptest+ eklo

		elseif(ituning.eq.3) then 
			if(ptest.lt.tunet2) uoff = (ekhi-eklo)/tunet2*ptest + eklo
			if(ptest.ge.tunet2) uoff = (eklo-ekhi)/tunet2*(ptest-tunet2) + ekhi

		elseif(ituning.eq.4) then 
			uoff = (ekhi-eklo)/2.0*sin(2.*pi*ptest/tunet) +(ekhi+eklo)/2.0

		end if
		
121		continue
		
		
!################################################################
		if(fptune.eq.0.0) goto 122
						
		timek = kturns*speedfact*trev
		tunet = 1.0/fptune
		tunet2 = tunet/2.0
		ptest = mod(timek,tunet)
		
		if(ituning.eq.0.) then
			goto 122
		
		elseif(ituning.eq.1) then
			if(ptest.lt.tunet2) epsoff = ephi
			if(ptest.ge.tunet2) epsoff = eplo

		elseif(ituning.eq.2) then 
			epsoff = (ephi-eplo)/tunet*ptest+ eplo

		elseif(ituning.eq.3) then 
			if(ptest.lt.tunet2) epsoff = (ekhi-eklo)/tunet2*ptest + eplo
			if(ptest.ge.tunet2) epsoff = (eklo-ekhi)/tunet2*(ptest-tunet2) + ephi

		elseif(ituning.eq.4) then 
			epsoff = (ephi-eplo)/2.0*sin(2.*pi*ptest/tunet) +(ephi+eplo)/2.0

		end if
		
122		continue
		!write(*,*) timek,epsoff
	
	end
	
	
	
	
	
	
!-----------------------------------------
!----   beam initialization --------------
	subroutine addparticles(ibeam)
		include 'ecoolm2h4_c.f90'
		external ran3
   		
		call beamgeneration(ibeam)

		
! initial beam distribution		
		open(unit=21,file='ini.full',status='unknown')
		do k=1,np
			write(21,100)x(k),px(k)/betavgx,y(k),py(k)/betavgy,t(k),pt(k)/gamma0/beta**2
		enddo
		close(21)
100     format(6e15.5)
        
        
		return
	end

	
	
	
!------------------------------------------------------------
!-------  6D partcle generation in phase space --------------
	
	subroutine beamgeneration(ibeam)
		include 'ecoolm2h4_c.f90'
		external ran3
        
 ! coasting beam initial       
        if(ibeam.eq.1) then 
        
			np=0
			do np=1,ndim
                
40				continue
                ! r1 and r2 are uniform on (-1,1)
				!r1 = 2*ran3(iseed)-1
				!r2 = 2*ran3(iseed)-1
                call random_number(xlc)
                r1 = 2.0*xlc-1.0
                call random_number(xlc)
                r2 = 2.0*xlc-1.0
                
				amp = r1*r1+r2*r2
! just keep the circle and cut at 5 sigma
				if((amp.ge.1).or.(amp.lt. 3.e-6))go to 40
				facc = sqrt(-2.*log(amp)/amp)

! do vertical
				y(np) = ampy*r1*facc
! px has same amplitude as x iegnuplot px = beta_L*x'
				py(np) = ampy*r2*facc
! other transverse variable, same initial emittance
 41				continue
! r1 and r2 are uniform on (-1,1)
				!r1 = 2*ran3(iseed)-1
				!r2 = 2*ran3(iseed)-1
                call random_number(xlc)
                r1 = 2.0*xlc-1.0
                call random_number(xlc)
                r2 = 2.0*xlc-1.0
				amp = r1*r1+r2*r2
! just keep the circle and cut at 5 sigma
				if((amp.ge.1).or.(amp.lt. 3.e-6))go to 41
				facc = sqrt(-2.*log(amp)/amp)
! use second amplitude for horizontal
				x(np)= ampx*r1*facc
				px(np)= ampx*r2*facc
                
! r1 and r2 are uniform on (-1,1)
                
				!r1 = 2*ran3(iseed)-1
                call random_number(xlc)
                r1 = 2.0*xlc-1.0
! use second amplitude for horizontal
                sigpt = coastbdpp*beta**2*gamma0	 !sigpt=dgamma=dp/p*gamma*beta^2
				pt(np)= gauss_random(0.0,sigpt)
				t(np)= trev*r1/2.0          
            enddo

			write(6,*)' sqrt(<(gamma-gamma0)^2>) = ',sigpt
			
			
!========================================================================        
!!!!!!!! bunched beam initial             
		else if(ibeam.eq.2) then
			    
! above transition
			phizero=pi

			np=0
			pcoeff2=pcoeff/omegarf
			pcoeff3 = pcoeff2/omegarf
			radharm2 = vrf2/vrf
! voltage for unit change in gamma
			v00 = aatom*E00*1.e6/qatom
			sigv=0

			omega0 = 2*pi/trev
			p1 = nharm*tauhat*omega0
			p2 = nharm2*tauhat*omega0

			hammax = (cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0) 
			write(6,*)' hammax,tcoeff = ',hammax,tcoeff
			ptmax = sqrt(abs(2*hammax/tcoeff))
! ndim is now the number of particles
			do np = 1,ndim
43              continue
! get longitudinal coordinates
				!tk = tauhat*(2*ran3(iseed)-1)
				!ptk = ptmax*(2*ran3(iseed)-1)
                call random_number(xlc)
				tk = tauhat*(2*xlc-1) 
                call random_number(xlc)
				ptk = ptmax*(2*xlc-1)

				p1 = nharm*tk*omega0
				p2 = nharm2*tk*omega0
				ham=0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
!  mmb may 2010      if(ham.gt.hammax)go to 43
				if(ham/hammax .gt.1)go to 43
				prob = 1-(ham/hammax)
				prob = prob**power

				test = ran3(iseed)
				if(prob.lt.test)go to 43
				pt(np) = ptk
				t(np)=tk

				sigv = sigv + pt(np)**2
 44				continue
! r1 and r2 are uniform on (-1,1)
				!r1 = 2*ran3(iseed)-1
				!r2 = 2*ran3(iseed)-1
                call random_number(xlc)
                r1 = 2.0*xlc-1.0
                call random_number(xlc)
                r2 = 2.0*xlc-1.0
				amp = r1*r1+r2*r2
! just keep the circle and cut at 5 sigma
				if((amp.ge.1).or.(amp.lt. 3.e-6))go to 44
				facc = sqrt(-2.*log(amp)/amp)

! do vertical
				y(np) = ampy*r1*facc
! px has same amplitude as x iegnuplot px = beta_L*x'
				py(np) = ampy*r2*facc
! other transverse variable, same initial emittance
 45				continue
! r1 and r2 are uniform on (-1,1)
				r1 = 2*ran3(iseed)-1
				r2 = 2*ran3(iseed)-1
				amp = r1*r1+r2*r2
! just keep the circle and cut at 5 sigma
				if((amp.ge.1).or.(amp.lt. 3.e-6))go to 45
				facc = sqrt(-2.*log(amp)/amp)
! use second amplitude for horizontal
				x(np)= ampx*r1*facc
				px(np)= ampx*r2*facc
				!write(22,*)t(np),pt(np)
			enddo
			sigv = sqrt(sigv/np)
			write(6,*)' sqrt(<(gamma-gamma0)^2>) = ',sigv
			
			
			
!========================================================================      			
! coasting pulsed beam initial       
        else if(ibeam.eq.3) then 
        
			np=0
			do np=1,ndim
                
46				continue
                ! r1 and r2 are uniform on (-1,1)
				!r1 = 2*ran3(iseed)-1
				!r2 = 2*ran3(iseed)-1
                call random_number(xlc)
                r1 = 2.0*xlc-1.0
                call random_number(xlc)
                r2 = 2.0*xlc-1.0
                
				amp = r1*r1+r2*r2
! just keep the circle and cut at 5 sigma
				if((amp.ge.1).or.(amp.lt. 3.e-6))go to 46
				facc = sqrt(-2.*log(amp)/amp)

! do vertical
				y(np) = ampy*r1*facc
! px has same amplitude as x iegnuplot px = beta_L*x'
				py(np) = ampy*r2*facc
! other transverse variable, same initial emittance
 47				continue
! r1 and r2 are uniform on (-1,1)
				!r1 = 2*ran3(iseed)-1
				!r2 = 2*ran3(iseed)-1
                call random_number(xlc)
                r1 = 2.0*xlc-1.0
                call random_number(xlc)
                r2 = 2.0*xlc-1.0
				amp = r1*r1+r2*r2
! just keep the circle and cut at 5 sigma
				if((amp.ge.1).or.(amp.lt. 3.e-6))go to 47
				facc = sqrt(-2.*log(amp)/amp)
! use second amplitude for horizontal
				x(np)= ampx*r1*facc
				px(np)= ampx*r2*facc
                
! r1 and r2 are uniform on (-1,1)
                
				!r1 = 2*ran3(iseed)-1
                call random_number(xlc)
                r1 = 2.0*xlc-1.0
! use second amplitude for horizontal
                sigpt = coastbdpp*beta**2*gamma0	 !sigpt=dgamma=dp/p*gamma*beta^2
				pt(np)= gauss_random(0.0,sigpt)
				t(np)= pulseiW*r1/2.0 + trev/2.0
				
				if(t(np).gt.trev/2.0) t(np)=-trev+t(np)    
				!write(*,*) 'pulseiW=',pulseiW
        
            enddo

!			write(6,*)' sqrt(<(gamma-gamma0)^2>) = ',sigpt			
!=========================================================================			
		endif
				
! initial number of macro-particles for normalization of current etc
		np=np-1
		np0 = np
	
	end

	
!--------------------------------------------------------
!----   beam reinjection for coasting beam --------------
! barrier bucket injection using pulse e-beam 
! half new injection and half stored beam 
	
	subroutine reinjectbeam(ibeam)
		include 'ecoolm2h4_c.f90'
		real(8) ytemp(nb),pytemp(nb),ttemp(nb),pttemp(nb),xtemp(nb),pxtemp(nb)
		real(8) yfresh(nb),pyfresh(nb),tfresh(nb),ptfresh(nb),xfresh(nb),pxfresh(nb)
		
		! time ready for new injection
		timek = kturns*speedfact*trev
		freshtest = mod(timek,gapit)
		if((timek.lt.gapit).OR.(freshtest.gt.speedfact*trev)) goto 51
		
		injN = injN + 1
		
		npnow = np
		!write(*,*) 'nptemp=',nptemp
 
		! take half of stored particles and reject particles for fresh injection 	
		kooo = 0
		ninjlosstemp = 0
		do k=1,npnow
			sparelimit = (trev-spareW)/2.0 
			! beam loss due to new injection		
			if(abs(t(k)).lt.sparelimit) then
				kooo = kooo + 1
				ninjlosstemp = ninjlosstemp
				ttemp(kooo)=t(k)
				pttemp(kooo)=pt(k)
				xtemp(kooo)=x(k)
				pxtemp(kooo)=px(k)
				ytemp(kooo)=y(k)
				pytemp(kooo)=py(k)
			else 
				ninjlosstemp = ninjlosstemp + 1 
			endif

			npstore = kooo		
		end do 
		
		ninjloss = ninjloss + ninjlosstemp
		
		call beamgeneration(ibeam)
		npnow = np0
		npfresh = npnow
		
		write(6,*) 'inject i = ',injN
		write(6,*) 'npstore,npfresh = ',npstore,npfresh
		! take half of fresh particles
		yfresh = y(1:npfresh)
		pyfresh = py(1:npfresh)
		xfresh = x(1:npfresh)
		pxfresh = px(1:npfresh)
		tfresh = t(1:npfresh)
		ptfresh = pt(1:npfresh)
		
		
		! conbined the stored and fresh beam
		npnew = npfresh + npstore
		y(1:npnew)  = [ytemp(1:npstore),  yfresh(1:npfresh)]
		py(1:npnew) = [pytemp(1:npstore), pyfresh(1:npfresh)]
		x(1:npnew)  = [xtemp(1:npstore),  xfresh(1:npfresh)]
		px(1:npnew) = [pxtemp(1:npstore), pxfresh(1:npfresh)]
		t(1:npnew)  = [ttemp(1:npstore),  tfresh(1:npfresh)]
		pt(1:npnew) = [pttemp(1:npstore), ptfresh(1:npfresh)]	
		np = npnew
		
		!np: stored particle number untill now
		!np0: particle numebr in one injection		
		!write(6,*) 'ddd',np,ninjloss
		
51		continue
	
	end
	

!234567890123456789012345678901234567890123456789012345678901234567890
!----------------------------------------------------!----------------------------------------------------
!234567890123456789012345678901234567890123456789012345678901234567890
	subroutine sptranssc(fracturn)
! this one does space charge too
		include 'ecoolm2h4_c.f90'
! bare phase advance with tune in (1,2)
		psi0x = 2*pi*(tunex -int(tunex)+1)*fracturn
		psi0y = 2*pi*(tuney -int(tuney)+1)*fracturn
! same chrom for both planes
		coeffchrom = fracturn*chrom*2*pi/(gamma0-1/gamma0)
! get space charge coefficients

!	clight = 2.998e8
! calculates the space charge kick assuming no coherent motions
		dt2 = thib/nresamp2
		thib2 = 0.5*thib
! put into units of current
		currcoeff = 1.602e-19*pnumber*qatom/(np0*dt2)
! currcoeff*avgline(k) is current
! rmsx and rmsy are rms beam sizes
! assume round
		sigrad2 = rmsx*rmsy
		if(sigrad2.eq.0)sigrad2=ampx*ampy
		etotdq = gamma0*E00*1.0e6*aatom/qatom
		wave = beta*2*pi/circ
! tune shift for full turn
		scoeff0 = 377*currcoeff/(8*pi*beta*sigrad2*(gamma0*wave)**2*etotdq*tunex)
		scoeff0 = scoeff0*fracturn*2*pi
		if(idosc.eq.0)scoeff0=0
		actsc = 0.25/sigrad2

		kooo = 0
		do k = 1,np
			tk = t(k)
			tk = (tk+thib2)/dt2
! use linear interpolation
			nlo=tk
			nhi=nlo+1
			fhi=tk-nlo
			flo=1-fhi
			scoeff = flo*avgline(nlo)+fhi*avgline(nhi)
			dqsc0 = scoeff*scoeff0


			ptk=pt(k)
			dpsichrom =  ptk*coeffchrom
			xk = y(k) 
			pk = py(k)

			x2k=x(k)
			p2k=px(k)

			act = xk*xk+pk*pk
			act2 = x2k*x2k+p2k*p2k
			dqsc = dqsc0/(1+act*actsc)
! phase for x
			psi = psi0y + dpsichrom -dqsc
			a11 = cos(psi)
			a12 = sin(psi)
			y(k)  = xk*a11 + pk*a12  
			py(k) = pk*a11 - xk*a12

			dqsc = dqsc0/(1+act2*actsc)

! phase for x2
			psi = psi0x  + dpsichrom   -dqsc  
			a11 = cos(psi)
			a12 = sin(psi)
			x(k)  = x2k*a11 + p2k*a12  
			px(k) = p2k*a11 - x2k*a12
			
			xkm = x(k)*sqrt(betamax/betavgx)
			ykm = y(k)*sqrt(betamax/betavgy)
! beam loss in transverse based on the aperture (betamax)		
			if(abs(xkm).lt.aperture .AND. abs(ykm).lt.aperture) then
				kooo = kooo + 1
				ntloss = ntloss
				t(kooo)=t(k)
				pt(kooo)=pt(k)
! other dimensions
				x(kooo)=x(k)
				px(kooo)=px(k)
				y(kooo)=y(k)
				py(kooo)=py(k)
			else 
				ntloss = ntloss + 1 
			endif
			np = kooo

		enddo
!      write(44,*)psi1,ptk*coeffchrom,ptk
		return
	end
!----------------------------------------------------

!----------------------------------------------------
	subroutine rfupdate(fmix,ibeam)
!ref: coefficients  d pt(k)/dn = pcoeff*(t(k)-trf/2) , d  t(k)/dn = tcoeff*pt(k)
!fmix: speed up factor
		include 'ecoolm2h4_c.f90'
        
        if(ibeam.eq.1.OR.ibeam.eq.3) then
			do k=1,np
				tk = t(k)
				ptk = pt(k)
! update the time
				tk = tk + fmix*tcoeff*ptk		
                if(tk.gt.thib/2.) tk=-thib/2.	!phase jump
                if(tk.lt.(-thib/2.)) tk=thib/2.
                t(k) = tk
                pt(k) = ptk 
			end do
            
        else if(ibeam.eq.2) then
            
			omega0 = 2*pi/trev
			v00 = E00*1.0e6*aatom/qatom
			koo = 0
			thib2 = thib/2
			do k=1,np
				tk = t(k)
				ptk = pt(k)
! update the time
				tk = tk + fmix*tcoeff*ptk
				p1 = tk*nharm*omega0
				p2 = tk*nharm2*omega0
				volt = vrf*sin(p1)+vrf2*sin(p2)
				dptrfk=volt
				ptk = ptk + fmix*volt/v00
! other dimensions
				xk =  y(k)
				pxk = py(k)
				x2k = x(k)
				px2k = px(k)

! particles with abs(tk).gt.thib/2 get lost
				if(abs(tk).le.thib2)then
					koo = koo + 1
					t(koo)=tk
					pt(koo)=ptk
					ham=0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
					! longitudinal hamiltonian
					xpavg(koo)=ham
! other dimensions
					y(koo)=xk
					py(koo)=pxk
					x(koo)=x2k
					px(koo)=px2k
					dptrf(koo)=dptrfk
				endif
			enddo
!	write(6,*)' np, koo',np,koo
			np = koo
 
! put the normaization in np+1
			p1 = 0.5*thib*nharm*omega0
			p2 = 0.5*thib*nharm2*omega0
			xpavg(np+1)= (cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
        end if
        
		return
  	end


!--------------------------------------------------------------------
!234567
	subroutine writit(iwrite)
	include 'ecoolm2h4_c.f90'          
	if(iwrite.ne.1)return
! line density mountain range
	call mountainr
	call writemoments 

	!fmix=0
	!call rfupdate(fmix)

	return
	end
!---------------------------------------------------
	subroutine writemoments
!234567
	include 'ecoolm2h4_c.f90'
	
	tavg = 0.
	xavg = 0.
	yavg = 0.
	ptavg = 0.
	pxavg = 0.
	pyavg = 0.
	do kk=1,np
		ptavg = ptavg + pt(kk)
		pyavg = pyavg + py(kk)
		pxavg = pxavg + px(kk)
		tavg = tavg + t(kk)
		yavg = yavg + y(kk)
		xavg = xavg + x(kk)
	enddo
	
	ptavg = ptavg/np
	pyavg = pyavg/np
	pxavg = pxavg/np
	tavg = tavg/np
	yavg = yavg/np
	xavg = xavg/np		
		
	csfull=0
	csfull2=0
	siggama=0
	siggama2=0
	sigtime=0
	sigtime2=0
!	write(6,*)' np = ',np
	
	
	do k0 = 1,np
		xk = y(k0)
		pxk = py(k0)
		xk2 = x(k0)
		pxk2 = px(k0)

		csfull = csfull + (xk-yavg)*(xk-yavg) + (pxk-pyavg)*(pxk-pyavg)	!y,py
		csfull2 = csfull2 + (xk2-xavg)*(xk2-xavg) + (pxk2-pxavg)*(pxk2-pxavg)	!x,px

		dpwake(k0)= (xk-yavg)*(xk-yavg) + (pxk-pyavg)*(pxk-pyavg)
		dpwake2(k0)= (xk2-xavg)*(xk2-xavg) + (pxk2-pxavg)*(pxk2-pxavg)

! longitudinal
		siggama = siggama + (pt(k0)-ptavg)**2
		tk = t(k0)
		sigtime = sigtime + (tk-tavg)*(tk-tavg)
		sigtime2 = sigtime2 + tk

	enddo
!	write(6,*)' sums',csfull,csfull2,siggama,sigtime,sigtime2
	!sigtime = sqrt( sigtime/np - (sigtime2/np)**2)
	sigtime = sqrt(sigtime/np)
	siggama = sqrt(siggama/np)
	csfull = csfull/np
	csfull2 = csfull2/np
    
    sigdponp = siggama/beta/beta/gamma0
    sigemitty = csfull/2./betavgy
    sigemittx = csfull2/2./betavgx
!234567
! stop the code if the beam is large enough
	if(csfull.gt.20)then
		write(6,*)'sumcsavg too big',csfull
		stop
	endif
103	format(i10,1p10e13.5)
100	format(1p14e15.5)
	write(6,104)kturns,np,ntloss,log10(sigemittx),log10(sigemitty)
104	format('kturns,np,ntloss,log10(emittx),log10(emitty)',3i10,1p10e11.3)
! peakilumi is peak current and fwhmlumi is fwhm pulse length in time
!	write(33,103)kturns,csfull,csfull2,sigtime,siggama,peakilumi,fwhmlumi,float(np)
    write(33,103)kturns,sigemittx,sigemitty,sigtime,sigdponp,xavg,yavg,ptavg/beta/beta/gamma0,pnum,float(ntloss),float(ninjloss)
!23456
! write the percentiles
	call writepercent
! write everything to tran.full if nwrite.lt.0
	if(nwrite.gt.0)return
	write(6,*)' writing csfull'
	call tunestuff
	open(unit=20,file='tran.full',status='unknown')
	do k=1,np
! 
		write(20,100)x(k),px(k)/betavgx,y(k),py(k)/betavgy,t(k),pt(k)/gamma0/beta**2,nperturn*xtune(k),dptrf(k),xpavg(k)
	enddo
	close(20)

	return
	end
!---------------------------------------------------
	SUBROUTINE SORT2(N,RA,RB)
	implicit real(8) (a-h,o-z)
	DIMENSION RA(N)
! NEXT DECLARATIONS FOR INTEGER TAG-ALONG ARRAY
! BASED ON QUICKSORT
	INTEGER RB(N),RRB
	L=N/2+1
	IR=N
10	CONTINUE
	IF(L.GT.1)THEN
		L=L-1
		RRA=RA(L)
		RRB=RB(L)
	ELSE
		RRA=RA(IR)
		RRB=RB(IR)
		RA(IR)=RA(1)
		RB(IR)=RB(1)
		IR=IR-1
		IF(IR.EQ.1)THEN
			RA(1)=RRA
			RB(1)=RRB
			RETURN
	  	ENDIF
	ENDIF
	I=L
	J=L+L
20	IF(J.LE.IR)THEN
		IF(J.LT.IR)THEN
			IF(RA(J).LT.RA(J+1))J=J+1
		ENDIF
		IF(RRA.LT.RA(J))THEN
			RA(I)=RA(J)
			RB(I)=RB(J)
			I=J
			J=J+J
      	ELSE
        	J=IR+1
      	ENDIF
    	GO TO 20
    ENDIF
	RA(I)=RRA
	RB(I)=RRB
	GO TO 10
	END
!--------------------------------------------------------------------------------
	subroutine keepturn
	include 'ecoolm2h4_c.f90'


! accumulates the running sums
	common/oldpart/xold(nb),pold(nb)
	if(nwrite.gt.0)return
	do k=1,np
		xk = y(k)
		pk = py(k)
		xok = xold(k)
		pok = pold(k)
		xkeep(1,k) = xkeep(1,k) + xk*xok
		xkeep(2,k) = xkeep(2,k) + xk*pok
		xkeep(3,k) = xkeep(3,k) + xok*xok
		xkeep(4,k) = xkeep(4,k) + xok*pok
		xkeep(5,k) = xkeep(5,k) + pok*pok
		xkeep(6,k) = xkeep(6,k) + pk*xok
		xkeep(7,k) = xkeep(7,k) + pk*pok
		pold(k)=pk
		xold(k)=xk
	enddo
	return
	end
!-----------------------------------------------------
	subroutine zerokeep
	include 'ecoolm2h4_c.f90'
	common/oldpart/xold(nb),pold(nb)

	!beam accumulation and lifetime decay 
	timek = kturns*speedfact*trev
	pnum = float(np)/float(np0)*pnumber
!	write(*,*) pnum,exp(-timek/plifetime)
	pnum = pnum*exp(-timek/plifetime)
	
	if(iwrite.ne.1)return
	if(nwrite.gt.0)return
	do k=1,np
		do m=1,7
			xkeep(m,k)=0
		enddo
		xold(k)=y(k)
		pold(k)=py(k)
	enddo
	
	return
    end

    
    
    
!-------------------------------------------------------------
	subroutine tunestuff
	include 'ecoolm2h4_c.f90'
! use the simplest way for now
!234567
	do n=1,np
		a11 = xkeep(3,n)
		a12 = xkeep(4,n)
		a21 = a12
		a22 = xkeep(5,n)
		det = a11*a22-a12*a21
		b11 = a22/det
		b22 = a11/det
		b12 = -a12/det
		b21 = -a21/det
		alf = xkeep(1,n)
		bet = xkeep(2,n)
		a = alf*b11 + bet*b12
		b = alf*b21 + bet*b22
		alf = xkeep(6,n)
		bet = xkeep(7,n)
		! = alf*b11 + bet*b12
		d = alf*b21 + bet*b22
! now modify to unit determinant
		det = abs(a*d - b*c)
		arg=0
		if(det.gt.0)arg = 0.5*(a+d)/sqrt(det)
		if(abs(arg).lt.1)then
		   tunen = acos(arg)
		else
		   tunen=0
		endif
		if(b.lt.0)tunen=-tunen
		xtune(n) = tunen/(2*pi)
	enddo
! get average tune and chrom
	sum1 = 0
	sum2 = 0
	sum3 =0
	sum4 =0
	do k=1,np
		xtk= xtune(k)
		dgk = pt(k)
		sum1 = sum1 + xtk
		sum2 = sum2 + xtk*dgk
		sum3 = sum3 + dgk
		sum4 = sum4 + dgk*dgk
	enddo
	sum1 = sum1/np
	sum2 = sum2/np
	sum3 = sum3/np
	sum4 = sum4/np
	det = sum4 - sum3*sum3
	a11 = sum4/det
	a12 = -sum3/det
	a21 = a12
	a22 = 1/det
	avgtune = a11*sum1 + a12*sum2
	avgchrom = a21*sum1 + a22*sum2
	write(6,*)'avgchrom, chrom, avgtune = ',avgchrom*gamma0,chrom,avgtune
!23456
	return
	end

!------------------------------------------------------------------------
!------------------------------------------------------------------------
REAL(8) FUNCTION RAN3(IDUM)
!         IMPLICIT REAL*4(M)
!         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
	PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
	common/ran3stuff/ MA(55),inext,inextp
	DATA IFF /0/
	IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
		IFF=1
		MJ=MSEED-IABS(IDUM)
		MJ=MOD(MJ,MBIG)
		MA(55)=MJ
		MK=1
		DO 11 I=1,54
			II=MOD(21*I,55)
			MA(II)=MK
			MK=MJ-MK
			IF(MK.LT.MZ)MK=MK+MBIG
			MJ=MA(II)
11		CONTINUE
		DO 13 K=1,4
			DO 12 I=1,55
				MA(I)=MA(I)-MA(1+MOD(I+30,55))
				IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12			CONTINUE
13      CONTINUE
		INEXT=0
		INEXTP=31
		IDUM=1
	ENDIF
	INEXT=INEXT+1
	IF(INEXT.EQ.56)INEXT=1
	INEXTP=INEXTP+1
	IF(INEXTP.EQ.56)INEXTP=1
	MJ=MA(INEXT)-MA(INEXTP)
	IF(MJ.LT.MZ)MJ=MJ+MBIG
	MA(INEXT)=MJ
	RAN3=MJ*FAC
	RETURN
	END
!234567
!------------------------------------------------------------------
	SUBROUTINE FOUR1(DATA,NN,ISIGN)
! COMMENTS BY MIKE BLASKIEWICZ
! DOES AN FFT 
!  DATA(2*K-1), K=1,..NN, ARE REAL PARTS OF INPUT ARRAY
!  DATA(2*K)  , K=1,..NN, ARE IMAGINARY PARTS OF INPUT ARRAY
! DOES FFT IN PLACE
! DPHI = ISIGN*(2*PI)/NN
! FT(J)  = SUM_{K=1}^NN  (D(2*K-1) + (0,1)*D(2*K))*EXP((0,1)*DPHI*(K-1)*(J-1))
! ON OUTPUT D(2*K-1) = REAL(FT(K)), D(2*K) = IMAG(FT(K))
	implicit real(8) (a-h,o-z)
	REAL(8) WR,WI,WPR,WPI,WTEMP,THETA
	DIMENSION DATA(*)
	N=2*NN
	J=1
	DO 11 I=1,N,2
		IF(J.GT.I)THEN
			TEMPR=DATA(J)
			TEMPI=DATA(J+1)
			DATA(J)=DATA(I)
			DATA(J+1)=DATA(I+1)
			DATA(I)=TEMPR
			DATA(I+1)=TEMPI
        ENDIF
		M=N/2
1		IF ((M.GE.2).AND.(J.GT.M)) THEN
			J=J-M
			M=M/2
			GO TO 1
		ENDIF
		J=J+M
11	CONTINUE
	MMAX=2
2	IF (N.GT.MMAX) THEN
		ISTEP=2*MMAX
		THETA=6.28318530717959D0/(ISIGN*MMAX)
		WPR=-2.D0*DSIN(0.5D0*THETA)**2
		WPI=DSIN(THETA)
		WR=1.D0
		WI=0.D0
		DO 13 M=1,MMAX,2
			DO 12 I=M,N,ISTEP
				J=I+MMAX
				TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
				TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
				DATA(J)=DATA(I)-TEMPR
				DATA(J+1)=DATA(I+1)-TEMPI
				DATA(I)=DATA(I)+TEMPR
				DATA(I+1)=DATA(I+1)+TEMPI
12			CONTINUE
			WTEMP=WR
			WR=WR*WPR-WI*WPI+WR
			WI=WI*WPR+WTEMP*WPI+WI
13 		CONTINUE
        MMAX=ISTEP
		GO TO 2
	ENDIF
	RETURN
	END
!---------------------------------------------------------

 
!234567
	subroutine ibslong
	include 'ecoolm2h4_c.f90'
	external ran3
! uses Piwinski's formulae for ibs, pg 126 in handbook

	rmsx=0
	rmsy=0
	rmsg=0
	rmst=0
! zero out the array for line density
	do k=1,nresamp2
		denlon2(k)=0
		dentran2(k)=0
	enddo
	dtsamp2 = thib/nresamp2
	thib2 = thib/2

	dxsamp2 = thixx/nresamp2
	thixx2 = thixx/2
	pnum = pnum
		
!average of the beam velocities
	ptavg = 0.
	pxavg = 0.
	pyavg = 0.
	do k=1,np
		ptavg = ptavg + pt(k)
		pyavg = pyavg + py(k)
		pxavg = pxavg + px(k)
	enddo
		
	ptavg = ptavg/np
	pyavg = pyavg/np
	pxavg = pxavg/np

	do k=1,np
! t=0 is stable fixed point
		tk = t(k)
		gk = pt(k)
		yk = y(k)
		pyk = py(k)
		rmst = rmst + tk*tk
		rmsg = rmsg + (gk-ptavg)*(gk-ptavg)
		rmsy = rmsy + yk*yk + (pyk-pyavg)*(pyk-pyavg)
		rmsx = rmsx + (px(k)-pxavg)**2 + x(k)**2
		tk = (tk+thib2)/dtsamp2
! just use nearest integer
		nk = nint(tk)
		if(nk.le.0)nk=1
		if(nk.gt.nresamp2)nk=nresamp2
		denlon2(nk)=denlon2(nk)+1

		xk = (x(k)+thixx2)/dxsamp2
		nkx = nint(xk)
		if(nxk.le.0)nxk=1
		if(nxk.gt.nresamp2)nxk=nresamp2
		dentran2(nk)=dentran2(nk)+1		
	enddo
! add into avgline(k)
! writeavg determines how many turns for counting
	writeavg = 1./nwrite0
       
	do k=1,nresamp2
		avgline(k) = (1-writeavg)*avgline(k)+writeavg*denlon2(k)
		avglinet(k) = (1.-writeavg)*avglinet(k) + writeavg*dentran2(k)
	enddo



	rmst = sqrt(rmst/np)
	rmsg = sqrt(rmsg/np)
	rmsx = sqrt(0.5*rmsx/np)
	rmsy = sqrt(0.5*rmsy/np)

! x2 is horizontal, x is vertical
	dponp = rmsg/(beta*beta*gamma0)

! nothing for no ibs
	if(fracibstot.le.0)return
      
!	clight = 2.998e8
! dispersion in smoooth approximation
	betay = circ/(2*pi*tuney)
	betax = circ/(2*pi*tunex)

	dispx = circ/(2*pi*gammat**2)
	sigh2inv =  1/dponp**2  + (dispx/rmsx)**2 
	sigh = 1/sqrt(sigh2inv)
! get capitol A
!	write(6,*)rmsx,rmsy,dponp,beta,sigh
! classical radius
	r0 = (qatom**2/aatom)*1.54e-18
	atop = r0*r0*clight*pnum

	epsx = rmsx*rmsx/betax
	epsy = rmsy*rmsy/betay
!	epsy = epsx
	sigs = rmst*clight*beta
    
    if(iichoice.eq.1.OR.iichoice.eq.3) sigs = pulseW*pulseN*circ/2./sqrt(pi)	!coasting&pulsed beam£¬for pulsed e-beam cooling, L_ion = pulseW
    
	abot = 64*pi*pi*beta**3*gamma0**4*epsx*epsy*sigs*dponp
! ca is Piwinski's A
	ca = atop/abot
! mohl's a,b,and q
! a is horizontal
	a = sigh*betax/(gamma0*rmsx)
! b is vertical
	b = sigh*betay/(gamma0*rmsy)
! log is good enough
	q = sigh*beta*sqrt(2*rmsy/r0)
! calculate fpiwin(a,b,q) with 1000 points
	npp=1000
	ceuler = 0.577215
!      pi = 3.14159265
!      write(6,*)' a, b, q ', a,b,q
	fpiwinp = fpiwin(a,b,q,npp)
	fpiwinx = fpiwin(1/a,b/a,q/a,npp)
	fpiwiny = fpiwin(1/b,a/b,q/b,npp)
	alfap0 =ca*fpiwinp*(sigh/dponp)**2
	alfax0 =ca*(fpiwinx+fpiwinp*(dispx*sigh/rmsx)**2)
	alfay0 =ca*fpiwiny
! alphas are amplitude growth rates for pnumber particles
! and no losses
! correct for small number of macro particles
! mult by 2 so it is a kick for emittance 
! use actual values times speedfact
	alfap=2*alfap0*fracibstot*speedfact
	alfax=2*alfax0*fracibsx*speedfact
	alfay=2*alfay0*fracibsx*speedfact
	  
! calculate the viscous forces ala Zenkevich
! anbunch is eq 14 with bunched beams
! coulomb log
	anbunch = log(2*rmsy/(r0*( gamma0*gamma0*( (rmsx/betax)**2 + (rmsy/betay)**2) + dponp**2)))
!      write(6,*)' anbunch = ',anbunch
! coulomb log * clight * r_i^2 * lambda
	anbunch = anbunch*clight*r0*r0*pnum/(2*sqrt(pi)*sigs)
!      write(6,*)' anbunch = ',anbunch
! divide by all but actions
	anbunch = anbunch/(2*(2*pi)**1.5 * beta**3 * gamma0**5)      
!	write(6,*)' anbunch = ',anbunch
! longitudinal      
		if(iichoice.eq.1.OR.iichoice.eq.3) then
			anbunch = anbunch/( dponp *sqrt(1.)/gamma0)   !coasting beam
		elseif(iichoice.eq.2) then
			anbunch = anbunch/( dponp *sqrt(2.)/gamma0)   !bunched beam
		endif    
!	write(6,*)' anbunch = ',anbunch
! both transverse
	anbunch = anbunch*betax*betay/(2*rmsx*rmsy)**2
!      write(6,*)' anbunch = ',anbunch
! a parameters
	a1z = gamma0**2 * sigh2inv
	a2z = betax*betay/(rmsx*rmsy)
! get the fint factors assuming a round beam
	call getfints(a1z,a2z,fint1,fint2)

    
	hatkappap = fracibstot*4*pi*fint1*anbunch*gamma0*gamma0/(dponp*dponp)
! use geometric mean 
	hatkappax = fracibsx*4*pi*fint2*anbunch*betax*betay/(rmsx*rmsy)
	hatkappay=hatkappax
! apply speedfact to hatkappas
	hatkappap = hatkappap*speedfact
	hatkappax = hatkappax*speedfact
	hatkappay = hatkappay*speedfact


! now need to add hatkappas to emittance growth rates
	alfap = alfap + hatkappap
	alfax = alfax + hatkappax
	alfay = alfay + hatkappay

! with rms growth rate for 1 dim of shm need to kick
! sqrt(2.) times harder
	if(alfap.ge.0)then
		coeffs = sqrt(2*alfap*trev)*rmsg
	else
		coeffs=0
	endif

	if(alfax.ge.0)then
		coeffx = sqrt(2*alfax*trev)*rmsx
	else
		coeffx = 0
	endif
	if(alfay.ge.0)then
		coeffy = sqrt(2*alfay*trev)*rmsy
	else
		coeffy = 0
	endif


	coeffmult = sigs*2*sqrt(pi)/(np*dtsamp2*clight)
	denlonn=0



	do k=1,np
		tk = (t(k)+thib2)/dtsamp2
! just use nearest integer
		nk = nint(tk)
		if(nk.le.0)nk=1
		if(nk.gt.nresamp2)nk=nresamp2
!          
!		denlon2k = denlon2(nk)*coeffmult
		denlon2k = avgline(nk)*coeffmult
		denlonn = denlonn+denlon2k
		denlon2k = sqrt(denlon2k)
!
		call get2gaussrv(iseed,grv1,grv2)
		dx = denlon2k*coeffx*grv1
		dg = denlon2k*coeffs*grv2
		call get2gaussrv(iseed,grv1,grv2)
		dy = denlon2k*coeffy*grv1
! decays for viscosity
		decayzenp = 1 - trev*hatkappap*denlon2k**2
		decayzenx = 1 - trev*hatkappax*denlon2k**2
		decayzeny = 1 - trev*hatkappay*denlon2k**2


		py(k)=py(k)*decayzeny +dy
		px(k)=px(k)*decayzenx +dx
		pt(k)=pt(k)*decayzenp+dg
		if(kturns.eq.3*nwrite)write(73,100)nint(tk),denlon2k,dx,dy,dg
	enddo
! get the central luminosity density
	clum=0
	do kkk=1,nresamp2
		clum = clum + denlon2(kkk)**2
	enddo
	ci = pnum/np
	clum = clum*ci*ci/(clight*dtsamp2*rmsx**2)
! turn number, rms quantities, scaled emittance growth rates,  

	itestit = kturns/nwrite
	itestit = itestit*nwrite
	if(itestit.eq.kturns)then
		write(66,100)kturns,rmsx,rmsy,rmsg/(gamma0*beta**2),alfax0,alfay0,alfap0,hatkappax,hatkappay,hatkappap
	endif
	return
 100	format(i10,1p12e14.4)
	end


!===========================================================
! B-M IBS model (Nagaitsev's method)
	subroutine ibslong2
	include 'ecoolm2h4_c.f90'
		external ran3
		integer ibstype 
		ibstype = 1

		rmsx=0
		rmsy=0
		rmsg=0
		rmst=0
! zero out the array for line density
		do k=1,nresamp2
			denlon2(k)=0
			dentran2(k)=0
		enddo
		dtsamp2 = thib/nresamp2
		thib2 = thib/2.

		dxsamp2 = thixx/nresamp2
		thixx2 = thixx/2.
		pnum = pnum

!average of the beam velocities
		ptavg = 0.
		pxavg = 0.
		pyavg = 0.
		do k=1,np
			ptavg = ptavg + pt(k)
			pyavg = pyavg + py(k)
			pxavg = pxavg + px(k)
		enddo
	
		ptavg = ptavg/np
		pyavg = pyavg/np
		pxavg = pxavg/np
		
		do k=1,np
! t=0 is stable fixed point
			tk = t(k)
			gk = pt(k)
			yk = y(k)
			pyk = py(k)
			rmst = rmst + tk*tk
			rmsg = rmsg + (gk-ptavg)*(gk-ptavg)
			rmsy = rmsy + yk*yk + (pyk-pyavg)*(pyk-pyavg)
			rmsx = rmsx + (px(k)-pxavg)**2 + x(k)**2
			tk = (tk+thib2)/dtsamp2
! just use nearest integer
			nk = nint(tk)
			if(nk.le.0)nk=1
			if(nk.gt.nresamp2)nk=nresamp2
			denlon2(nk)=denlon2(nk)+1

			xk = (x(k)+thixx2)/dxsamp2
			nkx = nint(xk)
			if(nkx.le.0)nkx=1
			if(nkx.gt.nresamp2)nkx=nresamp2
			dentran2(nkx)=dentran2(nkx)+1	
		enddo
! add into avgline(k)
		writeavg = 1./nwrite0
 
! smooth the long. distribution 
		do k=1,nresamp2
			avgline(k) = (1.-writeavg)*avgline(k) + writeavg*denlon2(k)
			avglinet(k) = (1.-writeavg)*avglinet(k) + writeavg*dentran2(k)
		enddo
		
		
		
		
		
		rmst = sqrt(rmst/np)
		rmsg = sqrt(rmsg/np)
		rmsx = sqrt(0.5*rmsx/np)
		rmsy = sqrt(0.5*rmsy/np)

		aperture = DA * rmsx
		beta = sqrt(1.-1./gamma0**2)
		
! dispersion in smoooth approximation
		betay = circ/(2*pi*tuney)
		betax = circ/(2*pi*tunex)

		dispx = circ/(2*pi*gammat**2)
		sigh2inv =  1.0/dponp**2  + (dispx/rmsx)**2		!make no sense for the wrong result
		sigh = 1.0/sqrt(sigh2inv)
	   
! x2 is horizontal, x is vertical
! emittance and bunch shape for IBS 
		dponp = rmsg/(beta*beta*gamma0)
		epsy = rmsy**2/betax
		epsx = rmsx**2/betay
		sigs = rmst*clight*beta
        
        if(iichoice.eq.1.OR.iichoice.eq.3) sigs = pulseW*pulseN*circ/2./sqrt(pi)	!coasting&pulsed beam£¬for pulsed e-beam cooling, L_ion = pulseW
		 
! nothing for no ibs
		if(fracibstot.le.0)return

! if(mod(kturns,ibstep).ne.0.0)	goto 331
		  	
		sumx=0.0
		sumy=0.0
		sump=0.0
		den=0.0
		 
		!circ = 429.727834240
		!frev = beta*clight/circ
		do k=1,nsliceIBS
		! read(10,*)bx,bxp,by,dx,dxp,frac
			bx = bxa(k)
			! bxp is d beta/ds
			bxp=bxap(k)
			byp=byap(k)
			by = bya(k)
			dx = dxa(k)
			dxp = dxpa(k)
			frac = fraca(k)
			dy = dya(k)
			dyp = dypa(k)
		!write(*,*) 'beforeget ',k,dx,frac
! get intrabeam scattering rates for emittances
!!! getrates:Nagaitsev without vertical dispertion
!!! getrates2:Borjken-Mtingua Model with vertical dispersion (Betacool)
!!! getrates3:Nagaitsev with vertical dispertion

			call getrates(alfx,alfy,alfp,bx,bxp,by,byp,dx,dxp,dy,dyp,ibstype)
		!write(*,*) 'ibslong2:  ', k,alfx,alfp
			sumx = sumx + frac*alfx
			sumy = sumy + frac*alfy
			sump = sump + frac*alfp
			den = den+frac
		!write(*,*) 'beam:  ',sumx,sumy,sump,frac
	!save the IBS heating rate
			if(kloop.eq.1.)write(67,101) frac,den,bx,by,dx,dy,alfx,alfy,alfp,sumx,sumy,sump
		enddo
		
		! normalize
		sumx = sumx/den
		sumy = sumy/den
		sump = sump/den
    
	
!331		continue		
		!write(*,*) 'ibslong2:  ',sumx,sumy,sump
		!write(*,*) 'beam:  ',epsx,epsy,dponp
! use actual values times speedfact
		alfap0=sump
		alfax0=sumx
		alfay0=sumy
	
		alfap=alfap0*fracibstot*speedfact
		alfax=alfax0*fracibsx*speedfact
		alfay=alfay0*fracibsx*speedfact
		 	
! calculate the viscous forces ala Zenkevich
! anbunch is eq 14 with bunched beams
! coulomb log
		r0 = (qatom**2/aatom)*1.54e-18
		
!-------------------------------------

		anbunch = log(2*rmsy/(r0*( gamma0*gamma0*( (rmsx/betax)**2 + (rmsy/betay)**2) + dponp**2)))
!		write(6,*)' anbunch = ',anbunch
! coulomb log * clight * r_i^2 * lambda
		!anbunch = anbunch*clight*r0*r0*pnum/(2*sqrt(pi)*sigs)
		anbunch = anbunch*clight*r0*r0*pnum/(2*sqrt(pi)*sigs)
		!write(6,*)' anbunch = ',anbunch
! divide by all but actions
		anbunch = anbunch/(2*(2*pi)**1.5 * beta**3 * gamma0**5)      
!		write(6,*)' anbunch = ',anbunch
! longitudinal      
		if(iichoice.eq.1.OR.iichoice.eq.3) then
			anbunch = anbunch/( dponp *sqrt(1.)/gamma0)   !coasting beam
		elseif(iichoice.eq.2) then
			anbunch = anbunch/( dponp *sqrt(2.)/gamma0)   !bunched beam
		endif
		       
!		write(6,*)' anbunch = ',anbunch
! both transverse
		anbunch = anbunch*betax*betay/(2*rmsx*rmsy)**2
!		write(6,*)' anbunch = ',anbunch
! a parameters
		a1z = gamma0**2 * (1.0/dponp**2  + (dispx/rmsx)**2)
		a2z = betax*betay/(rmsx*rmsy)

! get the fint factors assuming a round beam
		!write(*,*) 'test1 ',a1z,a2z,gamma0,sigh2inv
		call getfints(a1z,a2z,fint1,fint2)
        !write(*,*) a1z,a2z,fint11
		!call getfints2(a1z,a2z,fint1,fint2)
		!write(*,*) a1z,a2z,fint1
        
	!write(*,*) 'test2 ',a1z,a2z
		hatkappap = fracibstot*4.0*pi*fint1*anbunch*gamma0*gamma0/(dponp*dponp)
! use geometric mean 
		hatkappax = fracibsx*4*pi*fint2*anbunch*betax*betay/(rmsx*rmsy)
		hatkappay=hatkappax
		
! apply speedfact to hatkappas
		hatkappap = hatkappap*speedfact
		hatkappax = hatkappax*speedfact
		hatkappay = hatkappay*speedfact

! now need to add hatkappas to emittance growth rates
		alfap = alfap + 1.*hatkappap
		alfax = alfax + 1.*hatkappax
		alfay = alfay + 1.*hatkappay
	!write(*,*) "viscous2: ", hatkappax2,hatkappax,hatkappap
		alfap = alfap
		alfax = alfax
		alfay = alfay
	
! with rms growth rate for 1 dim of shm need to kick
! sqrt(2.) times harder
        
		if(alfap.ge.0)then
			coeffs = sqrt(2.*alfap*trev)*rmsg
		else
			coeffs=0
		endif

		if(alfax.ge.0)then
			coeffx = sqrt(2.*alfax*trev)*rmsx
		else
			coeffx = 0
		endif
		if(alfay.ge.0)then
			coeffy = sqrt(2.*alfay*trev)*rmsy
		else
			coeffy = 0
		endif
	!write(*,*) "ibslong: ",coeffx2,coeffx,coeffs

		coeffmult = sigs*2*sqrt(pi)/(np*dtsamp2*clight)
		denlonn=0

		do k=1,np
			tk = (t(k)+thib2)/dtsamp2
! just use nearest integer
            nk = nint(tk)
            if(nk.le.0)nk=1
            if(nk.gt.nresamp2)nk=nresamp2
		
		!wrong position
      	! write(*,*) "ibs: ",nk,avgline(nk),coeffmult
!		denlon2k = denlon2(nk)*coeffmult
			denlon2k = avgline(nk)*coeffmult
			denlonn = denlonn+denlon2k
			denlon2k = sqrt(denlon2k)
		!denlon2k = 1.0
		!write(*,*) tk,avgline(nk),denlon2k
		 
			call get2gaussrv(iseed,grv1,grv2)
			dx = denlon2k*coeffx*grv1
			dg = denlon2k*coeffs*grv2	
			call get2gaussrv(iseed,grv1,grv2)
			dy = denlon2k*coeffy*grv1
! decays for viscosity
			decayzenp = 1 - trev*hatkappap*denlon2k**2
			decayzenx = 1 - trev*hatkappax*denlon2k**2
			decayzeny = 1 - trev*hatkappay*denlon2k**2

			px(k)=px(k)*decayzenx+dx*1.0
			py(k)=py(k)*decayzeny+dy*1.0
			pt(k)=pt(k)*decayzenp+dg*1.0
            
		
		!px(k)=px(k)
		!px2(k)=px2(k)
		!pt(k)=pt(k)
			
        end do
        
        !write(*,*) px(424),px(424)*(1-decayzenx),denlon2k*coeffx
! get the central luminosity density
		clum=0
		do kkk=1,nresamp2
			clum = clum + denlon2(kkk)**2
		enddo
		ci = pnum/np
		clum = clum*ci*ci/(clight*dtsamp2*rmsx**2)
! turn number, rms quantities, scaled emittance growth rates,  

		itestit = kturns/nwrite
		itestit = itestit*nwrite
		if(itestit.eq.kturns)then
			write(66,100)kturns,rmsx,rmsy,rmsg/(gamma0*beta**2),alfax0,alfay0,alfap0,hatkappax,hatkappay,hatkappap
		endif
		return
100		format(i10,1p12e14.4)
101 	format(13e14.6)
	end


!----------------------------------------------
!-------------------------------------------------
	subroutine mountainr
! makes a mountain range plot of emittance versus 
! dgamma
		include 'ecoolm2h4_c.f90'
		dt2 = thib/nresamp2
		dt22 = thixx/nresamp2
! put into units of current
! old      currcoeff = 1.602e-19*pnumber*qatom/(np0*nwrite0*dt2)
!      currcoeff = 1.602e-19*pnumber*qatom/(nwrite0*dt2)
		currcoeff = 1.602e-19*pnum*qatom/(dt2)
		currcoeff = currcoeff/np
		
		currcoeff2 = 1.602e-19*pnum*qatom/(dt22)
		currcoeff2 = currcoeff2/np
!	write(6,*)currcoeff,nwrite0,np0
		peaki=0
		jdex=0
		do k=1,nresamp2
!			currk = avgline(k)*currcoeff
			currk = avgline(k)
			if(currk.gt.peaki)then
				peaki=currk
				jdex=k
			endif
!		avgline(k)=currk
		enddo
! get the full width half max of the peak
		do k=jdex,nresamp2
			pk = avgline(k)
			if(pk.le.peaki/2)go to 33
		enddo
 33		continue
		jhi=k
		do k=jdex,1,-1
			pk = avgline(k)
			if(pk.le.peaki/2)go to 34
		enddo
 34		continue
		jlo=k
		fwhmlumi = dt2*(jhi-jlo+1)
		peakilumi=peaki*currcoeff
		
        write(62,*) "mountainr"
		do k=1,nresamp2
			xxx = currcoeff*avgline(k)
			xxx2 = currcoeff2*avglinet(k)
			if(xxx.lt.1.e-10)xxx=0
			write(62,100)(k-nresamp2/2.0)*dt2,xxx,(k-nresamp2/2.0)*thixx2,xxx2,float(kturns)
		enddo
		
 100	format(1p10e14.4)
	return
    end
    
    
    
!------------------------------------------------------------------------------------
	subroutine get2gaussrv(iseed,grv1,grv2)
	implicit real(8) (a-h,o-z)
! output is two zero mean gaussian rvs with rms = 1
! get two uniform deviates in the unit circle
	real(8)::grv1,grv2
144     continue
! r1 and r2 are uniform on (-1,1)
	!r1 = 2*ran3(iseed)-1
	!r2 = 2*ran3(iseed)-1
    call random_number(xlc)
	r1 = 2.0*xlc-1.0 
    call random_number(xlc)
	r2 = 2.0*xlc-1.0
	amp = r1*r1+r2*r2
! just keep the circle and cut at 6 sigma
	if((amp.ge.1).or.(amp.lt. 1.e-8))go to 144
	facc = sqrt(-2.*log(amp)/amp)
	grv1 = r1*facc        
	grv2 = r2*facc
	return
	end
!234567
!-------------------------------------------------------------------------
real(8)		function fpiwin(a,b,q,np)
! direct crib from page 126 in handbook
	real(8):: a,b,q
	integer np
	ceuler = 0.577215
	pi = 3.14159265
	sum = 0
	du = 1/float(np)
	do k=0,np
		u = k*du
		cp = sqrt(a*a+(1-a*a)*u*u)
		cq = sqrt(b*b+(1-b*b)*u*u)
		dsum = 2*log(q*(1/cp+1/cq)/2) - ceuler
		dsum = dsum*(1-3*u*u)/(cp*cq)
		if(k.eq.0)dsum=dsum/2
		if(k.eq.np)dsum=dsum/2
		sum = sum + dsum
	enddo
	fpiwin = 8*pi*du*sum
!      write(6,*)'a,b,q,fpiwin',a,b,q,fpiwin
	return
	end
!----------------------------------------
	subroutine sptrans2
! does the thin skew quad
	include 'ecoolm2h4_c.f90'
	tpdqmin = 2*pi*dqmin
	do k = 1,np
		px(k)=px(k)+tpdqmin*y(k)
		py(k)=py(k)+tpdqmin*x(k)
	enddo
!      write(44,*)psi1,ptk*coeffchrom,ptk
	return
	end
!------------------------------------------------------------------------------------------
	subroutine getfints(a1z,a2z,fint1,fint2)
	implicit real(8) (a-h,o-z)
! fint1 = integral(0,inf) dx  sqrt(x/(x+b))/((x+a)*(x+b))
! fint2 = integral(0,inf) dx  sqrt(x/(x+b))/((x+a)*(x+a))
! b=a1z, a = a2z
	real(8) a1z,a2z,fint1,fint2
	a = a2z
	b = a1z
! get int1 analyically
	ratio=(b-a)/a
	if(ratio.gt.0.01)then
		x=sqrt(ratio)
		fint1 = 2*(1-atan(x)/x)/(x*x*a)
	endif
	if(ratio.lt.-0.01)then
		x=sqrt(-ratio)
		fint1 = -2*(1-0.5*log((1+x)/(1-x))/x)/(x*x*a)
	endif
	if(abs(ratio).le.0.01)then
		x = ratio
		fint1 = 1./3. - x/5 +x*x/7-x*x*x/9
		fint1 = 2*fint1/a 
	endif
!       write(6,*)' exact sum1 = ',fint1


	if(ratio.gt.0.01)then
		x=sqrt(ratio)
		fint2 = ((b/a)*atan(x)/x -1)/(b-a)
	endif

	if(ratio.lt.-0.01)then
		x=sqrt(-ratio)
		fint2 = ((b/a)*0.5*log((1+x)/(1-x))/x -1  )   /(b-a)
	endif

	if(abs(ratio).le.0.01)then
		x = ratio
		fint2=2./3 +(0.2-1/3.)*x + (0.2-1/7.)*x*x +(1/9.-1/7.)*x*x*x
		fint2 = fint2/a
	endif
	return
    end
!-------------------------------------------------------------------------------------------
!234567
!-----------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
	subroutine getfints2(a1z,a2z,fint1,fint2)
	implicit real(8) (a-h,o-z)
! fint1 = integral(0,inf) dx  sqrt(x/(x+b))/((x+a)*(x+b))
! fint2 = integral(0,inf) dx  sqrt(x/(x+b))/((x+a)*(x+a))
! b=a1z, a = a2z
	real(8) a,b,a1z,a2z,fint1,fint2
	a = a2z
	b = a1z
! get int1 analyically
	fint1 = sqrt(a*(b-a))*acos(sqrt(a/b))
    !write(*,*)"ee", fint1,b-a
    fint1 = -2.0*(a-b+fint1)
   ! write(*,*)"ee1", fint1
    fint1 = fint1/(a-b)/(a-b)
    !write(*,*)"ee2", fint1
!       write(6,*)' exact sum1 = ',fint1
    
! get int2 analyically
    fint2 = sqrt((b-a)/a)*b*acos(sqrt(a/b))
    fint2 = a-b+fint2
    fint2 = fint2/(a-b)/(a-b)  
    !write(*,*) "ee3",fint2

	return
	end
!-------------------------------------------------------------------------------------------
!234567
!-----------------------------------------------------------------------------------
	subroutine writepercent
	include 'ecoolm2h4_c.f90'
	dimension jindex(nb)
! use heap sort on dpwake and dpwake2
	call sort(np,dpwake)
	call sort(np,dpwake2)
! use simple percentiles without funky edges
	write(71,100)kturns,(dpwake( nint(k*0.05*np)),k=1,20)
	write(72,100)kturns,(dpwake2( nint(k*0.05*np)),k=1,20)
 100  format(i10,1p20e11.3)
	return
	end

!--------------------------------------------------------
	SUBROUTINE SORT(N,RA)
	implicit real(8) (a-h,o-z)
	DIMENSION RA(N)
	L=N/2+1
	IR=N
 10    CONTINUE
	IF(L.GT.1)THEN
		L=L-1
		RRA=RA(L)
	ELSE
		RRA=RA(IR)
		RA(IR)=RA(1)
		IR=IR-1
		IF(IR.EQ.1)THEN
			RA(1)=RRA
		RETURN
		ENDIF
	ENDIF
	I=L
	J=L+L
 20	IF(J.LE.IR)THEN
		IF(J.LT.IR)THEN
			IF(RA(J).LT.RA(J+1))J=J+1
	  	ENDIF
	  	IF(RRA.LT.RA(J))THEN
			RA(I)=RA(J)
			I=J
			J=J+J
	  	ELSE
			J=IR+1
	  	ENDIF
		GO TO 20
	ENDIF
	RA(I)=RRA
	GO TO 10
	END
!------------------------------------------------------
!23456
!      FUNCTION BESSI0(X)
!      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
!     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
!      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
!     *0,
!     *    0.2659732D0,0.360768D-1,0.45813D-2/
!      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
!     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
!     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
!      IF (ABS(X).LT.3.75) THEN
!        Y=(X/3.75)**2
!        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
!      ELSE
!        AX=ABS(X)
!        Y=3.75/AX
!        BESSI0=(EXP(AX)/sqrt(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
!     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
!      ENDIF
!      RETURN
!      END
!--------------------------------------------------------------------------------------------
!234567
	subroutine ebeamkick
! does all the calculations for electron cooling
! assumes the map is a full turn
! backs up all particles half the length of the IP
! drifts with electron cooling force from -slencool/2 to slencool/2
! backs up particles to the IP
	include 'ecoolm2h4_c.f90'

! the tracking goes to the center of the IP
! transform to x, ux  coordinates and track back to entrance
! x will be position in meters, not betatron coordinates
		cx = sqrt(betaxcool/betavgx)
! px with be vx/! in rest frame
		cux = gamma0/sqrt(betavgx*betaxcool)
		cy = sqrt(betaycool/betavgy)
		cuy = gamma0/sqrt(betavgy*betaycool)
		elldg = slencool/2/gamma0
!	write(6,*)' cuy = ',cuy
		do k=1,np
! get physical coordinates at center of cooling section
			xk = cx*x(k)+  dispxh*pt(k)/gamma0/beta**2
			pxk = cux*px(k)
			yk = cy*y(k)
			pyk = cuy*py(k)
! trace back to start of cooling section
			xk = xk -elldg*pxk
			yk = yk -elldg*pyk
! OK, we have physical coordinates just before interacting with electrons
			x(k)=xk
			y(k)=yk
			px(k)=pxk
			py(k)=pyk
! include longitudinal
			pt(k)=pt(k)/gamma0/beta**2
		enddo 
! there are nslice cooling kicks there are two updates 
		elldg = slencool/2/gamma0/nslice
! s coordinate with s=0 in center of IP
		sslice=-slencool/2
		do k=1,nslice
			call driftem(elldg)
!		write(6,*)'coolslice1',px(np/4),px(np/2)
! next routine performs a cooling kick
			sslice = sslice + slencool/2/nslice
			call coolslice(sslice)
!		write(6,*)'coolslice',px(np/4),px(np/2)
!		write(6,*) k,nslice
			call driftem(elldg)
			sslice = sslice + slencool/2/nslice
		enddo
! drift back to center and get back betatron coordinates
		elldg = slencool/2/gamma0
!	write(6,*)'cx,cy,cux,cuy',cx,cy,cux,cuy
		do k=1,np
! physical coordinates at of cooling section
			pxk = px(k)
			xk = x(k) 
			yk = y(k)
			pyk = py(k)
! trace back to center of cooling section
			xk = xk - elldg*pxk
			yk = yk - elldg*pyk
! OK switch to betatron coordinates, pt has not been updated yet
			x(k)=(xk-dispxh*pt(k))/cx
			y(k)=yk/cy
			px(k)=pxk/cux
			py(k)=pyk/cuy
			pt(k)= gamma0*pt(k)*beta**2
		enddo 
! all done
		return
    end
    
    
    
!---------------------------------------------------
	subroutine driftem(elldg)
    include 'ecoolm2h4_c.f90'
		do k=1,np
			x(k) = x(k) + elldg*px(k)
			y(k) = y(k) + elldg*py(k)
		enddo
		return
	end

!---------------------------------------------------------------
	subroutine coolslice(scool)
! scool is the value of s at which to evaluate e beam parameters
	include 'ecoolm2h4_c.f90'
! length of cooled section
		dscool = slencool/nslice
! beta functions at s=0
		betaxe0 = ebeamx*gamma0/sigurestx
		betaye0 = ebeamy*gamma0/siguresty
! rms unnormalized electron emittance
		epsex = ebeamx*sigurestx/gamma0
		epsey = ebeamy*siguresty/gamma0
! CS functions where cooling is applied

		betaxe = betaxe0 + scool**2 /betaxe0
		betaye = betaye0 + scool**2 /betaye0

		alfaxe = -scool/betaxe0
		alfaye = -scool/betaye0
! constant dispersion in cooling section for electrons and hadrons
		sigx2 = betaxe*epsey + (sigurestt*dispxe)**2
		sigy2 = betaye*epsey
		sigt = ebeams/beta/2.998e8
		sigt2 = sigt**2
! cool with a single bunch of length ebeams
! get the coefficients of the v/! in rest frame distribution
		sigma12 = epsex *gamma0*gamma0*(1 + (alfaxe*dispxe*sigurestt)**2 / sigx2)/betaxe
		sigma22 = epsey*gamma0*gamma0/betaye
		sigma32 = sigurestt**2 *epsex*betaxe/sigx2
		rho = alfaxe*dispxe*sigurestt/sqrt(sigx2 + (alfaxe*dispxe*sigurestt)**2)
! uzbar= coeffbaruz*x, uxbar = x*coeffbarux, uybar = y*coeffbaruy
		coeffbaruz = dispxe*sigurestt**2 / sigx2
		coeffbaruy = -gamma0*alfaye/betaye
		coeffbarux = -gamma0*alfaxe*epsex/(betaxe*epsex+(dispxe*sigurestt)**2)

! coefficients of the distribution
		coeffa = 1/(sigma12*(1-rho*rho))
		coeffb = 1/(sigma32*(1-rho*rho))
		coeffd = 1/sigma22
		coeffc = rho/(sqrt(sigma12*sigma32)*(1-rho**2))
! following look OK
!      write(6,*)' a thru d', coeffa,coeffb,coeffc,coeffd,scool,alfaxe,rho,sigurestt
!       write(6,*)'coeffbars',coeffbaruz,coeffbaruy,coeffbarux
! get vectors (v11,v12) and (v21,v22) dots with (ux-uxbar,uz-uzbar)
! axis transfer: have beta1*(ux-alfa*uz)**2 + beta2*(uz+ux*alfa)**2 = ca*ux**2 + cb*uz**2 - 2*cc*ux*uz
		call getvs(v11,v12,v21,v22,coeffa,coeffb,coeffc)
! below look OK
!      write(6,*)'vs',v11,v12,v21,v22
! now need to find the smallest of the 3 sigmas. 

		sig1 = 1/sqrt(v11**2 + v12**2)
		sig2 = 1/sqrt(v21**2 + v22**2)
		sig3 = sqrt(sigma22)
! below look OK
!	write(6,*)' sigs',sig1,sig2,sig3
		small=sig1  
		ns=1
		if(sig2.lt.small)then
			small=sig2
			ns=2
		endif

		if(sig3.lt.small)then
			small=sig3
			ns=3
		endif
! next looks OK
!	write(6,*)' ns small',ns,small

! ns is the index of the small coordinate, fill in the 3 arrays with the coefficients
		if(ns.eq.1)then
			nsum1=0
			caray1(0)= 1
		else
			call getcaray(caray1,sig1,small,nsum1)            
		endif

		if(ns.eq.2)then
			nsum2=0
			caray2(0)= 1
		else
			call getcaray(caray2,sig2,small,nsum2)            
		endif

		if(ns.eq.3)then
			nsum3=0
			caray3(0)= 1
		else
			call getcaray(caray3,sig3,small,nsum3)            
		endif
!      write(6,*)'nsums',nsum1,nsum2,nsum3
!      write(6,*)
!      write(6,*)(caray1(m),m=-nsum1,nsum1)
!      write(6,*)
!      write(6,*)(caray2(m),m=-nsum2,nsum2)
!      write(6,*)
!      write(6,*)(caray3(m),m=-nsum3,nsum3)
!      write(6,*)

! OK we now have the 3 arrays
! get the coefficient on page 13 on may 2019 notes
		re = 2.818e-15
! coloumb log = clog. just calculate at s=0
		bmax = 2*sqrt(sigy2)
! electron number density in rest frame
		ebeamrest = qbeamcool/(1.602e-19*(2*pi)**1.5 * sqrt(sigy2*sigx2) *ebeams *gamma0)

		ebeamrest = abs(ebeamrest)
!	write(6,*) 'ebeamrest = ',ebeamrest
! transverse rms of v/! for hadrons in rest frame
		rmsux2 = (rmsx*gamma0)**2/(betaxcool*betavgx)
		rmsuy2 = (rmsy*gamma0)**2/(betaycool*betavgy)
		rmsrelbeta2 = sigurestt**2+ sigurestx**2 + siguresty**2 + rmsux2+rmsuy2 + dponp**2
		debye = rmsrelbeta2/(8*pi*2.818e-15*ebeamrest)
		debye = sqrt(debye)
!	write(6,*)' debye = ',debye
		if(debye.lt.bmax)bmax=debye
! for short collisions can have smaller bmax
		bmax1 = sqrt(rmsrelbeta2)*slencool/(beta*gamma0)

		if(bmax1.lt.bmax)bmax=bmax1
		bmin = 2*qatom*re/rmsrelbeta2
		clog = log(bmax/bmin)
! max number density in lab frame
		denlabmax = qbeamcool/(1.602e-19*sqrt(sigx2*sigy2*(2*pi)**3)*ebeams)
		ratmass=1/(1836.*aatom)
		ducoeff = dscool*4*sqrt(2*pi)*qatom**2*re**2*ratmass*denlabmax*clog/(sig1*sig2*sig3*gamma0**2)
!      write(6,*)' ducoeff = ',ducoeff,dscool,clog
! unew = uold - ducoeff*sum_k C_k ... on page 13
! unit vector corresponding to sig1 multiplied by small sigma          
		h11 = v11*sig1*small
		h12 = v12*sig1*small
! unit vector corresponding to sig2           
		h21 = v21*sig2*small
		h22 = v22*sig2*small
      !write(6,*)'h s',h11,h12,h21,h22
		sqpid2 = sqrt(pi/2.)
		sqhalf = sqrt(0.5)
!234567
! apply speedup factor
		ducoeff = ducoeff*speedfact

!      write(6,*)'nsums',nsum1,nsum2,nsum3
! don't forget correction for spatial density
		do k=1,np
			xk = x(k)
			yk = y(k)
			tk = t(k)
			pxk = px(k)

			pyk = py(k)
			pzk = pt(k)
! get sum_k C_k ... on page 13
			sumx=0
			sumy=0
			sumz=0
! average electron velocitty
			uzbar= coeffbaruz*xk
			uxbar = xk*coeffbarux
			uybar = yk*coeffbaruy
! subtract average electron velocity
			 pxk1 = pxk - uxbar  - xpoff
!xpoff only in horizontal
			 pyk1 = pyk - uybar
			 pzk1 = pzk - uzbar

			 do n1 = -nsum1,nsum1
		 		do n2 = -nsum2,nsum2
		 			do n3 = -nsum3,nsum3
            			cktot = caray1(n1)*caray2(n2)*caray3(n3)
						offx = n1*h11 + n2*h21
						offz = n1*h12 + n2*h22
						offy = n3*small
						uxarg = pxk1-offx
						uyarg = pyk1-offy
! include dp/p offset of electrons
            			uzarg = pzk1-offz -uoff
						!write(*,*) uoff
! get I(R)/R**3
						r = sqrt(uxarg**2 + uyarg**2 +uzarg**2)
						rr = r/small
						fdivr3=0
						if(rr.lt.0.01)then
							fidivr3 = 1./3. - 0.1*rr*rr
						else
							if(rr.gt.6)then
								fidivr3 = sqpid2/rr**3
							else
		           				fidivr3 = (sqpid2*erf(rr*sqhalf)-rr*exp(-0.5*rr*rr))/rr**3
		          			endif
            			endif
!            if(scool.lt.-1)write(24,102)rr,fidivr3,n1,n2,n3,k
						fidivr3 = fidivr3*cktot
						sumx = sumx + uxarg*fidivr3
						sumy = sumy + uyarg*fidivr3
						sumz = sumz + uzarg*fidivr3
         			enddo
				enddo
			enddo
         
			if((abs(scool).lt.0.1).and.(kturns.eq.5))write(23,101)float(k),sumx,sumy,sumz,pxk,pyk,pzk,ducoeff
! spatial density, include offset of electron bunch
			xk11 = xk -xoff
			space = exp(-0.5*(xk11*xk11/sigx2 + yk*yk/sigy2 + tk*tk/sigt2))
			dus = ducoeff*space
			px(k)=pxk -dus*sumx
			py(k)=pyk -dus*sumy
			pt(k)=pzk -dus*sumz
        enddo      
! cooling kicks done
 101  format(10e12.3)
 102  format(2e13.3,4i5)
      return
      end

    


!----------------------------------------------------------------
	subroutine getvs(v11,v12,v21,v22,ca,cb,cc)
	implicit real(8) (a-h,o-z)
! if a and be are within 0.1% set them equal
	dar = abs(ca-cb)/(ca+cb)
	if(dar.lt.1.e-4)then
		alfa=-1
		beta1 = 0.5*(ca+cc)
		beta2 = 0.5*(ca-cc)
	else
		if(ca.gt.cb)then
			alfa = 2*cc/(abs(ca-cb)+sqrt( 4*cc**2 + (ca-cb)**2))
			beta1 = 0.5*(ca+cb)/(1+alfa**2) + 0.5*(ca-cb)/(1-alfa**2)
			beta2 = 0.5*(ca+cb)/(1+alfa**2) - 0.5*(ca-cb)/(1-alfa**2)
		else
			alfa = -2*cc/(abs(ca-cb)+sqrt( 4*cc**2 + (ca-cb)**2))
			beta1 = 0.5*(ca+cb)/(1+alfa**2) + 0.5*(ca-cb)/(1-alfa**2)
			beta2 = 0.5*(ca+cb)/(1+alfa**2) - 0.5*(ca-cb)/(1-alfa**2)
		endif
	endif
! have beta1*(ux-alfa*uz)**2 + beta2*(uz+ux*alfa)**2 = ca*ux**2 + cb*uz**2 - 2*cc*ux*uz
	v11 = sqrt(beta1)
	v12 = -v11*alfa

	v22 = sqrt(beta2)
	v21 = v22*alfa
	return
	end
!---------------------------------------------------------------------
!234567
	subroutine getcaray(caray,sigin,small,nsum)
	implicit real(8) (a-h,o-z)
	parameter(nv=401)
	real(8) caray(-200:200),fn(nv),akn(nv,nv),cn(nv)

! fits
! exp(-0.5*(x/sigin)**2) = \sum_{m=-nsum}^{nsum} c_m exp(-0.5*((m- x/small)**2)
	sig = sigin/small
	one=1
	pi = 4*atan(one)
	if(sig -1 .lt. 1.e-3)then
		nsum=0
		caray(0)=1
	return
	endif
! OK least squares
! last term is at 1.e-2 of peak, neglected term is < 1.e-3   
	nsum=1+3*sig
	ndim = 2*nsum+1
	do n=-nsum,nsum
		ndex = 1+n+nsum
		fn(ndex) = sqrt(2.*pi/(1+1/sig**2))/exp(0.5*n**2/(1+sig**2))
		do k=-nsum,nsum
			kdex = 1+k+nsum   
			akn(kdex,ndex) = sqrt(pi)/exp(0.25*(k-n)**2)
		enddo
	enddo
	  !write(*,*) 'nsum:',sig,sigin,small
! matrix is filled, use Gauss-Jordan elimination
	k1=1
	call gaussj(akn,ndim,nv,fn,k1,k1)
! solution is now in fn
	do k=1,ndim
		k1 = k-nsum-1
		caray(k1)=fn(k)
	enddo
	return
	end



       

!------------------------------------------------------------------------------
	SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
!* linear equation solution by Gauss-Jordan elimination
!* on imput A(N,N) is a matrix stored in an array NPxNP
!*          B(N,M) is a matrix containing the M rhs vectors stored in NPxMP
!* on output A is replaced by its matrix inverse and B by the solution
	implicit real(8) (a-h,o-z)
	PARAMETER (NMAX=50)
	DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
	DO 11 J=1,N
		IPIV(J)=0
11	CONTINUE
	DO 22 I=1,N
		BIG=0
        DO 13 J=1,N
			IF(IPIV(J).NE.1)THEN
			DO 12 K=1,N
				IF (IPIV(K).EQ.0) THEN
					IF (ABS(A(J,K)).GE.BIG)THEN
						BIG=ABS(A(J,K))
						IROW=J
						ICOL=K
		            ENDIF
				ELSE IF (IPIV(K).GT.1) THEN
					PAUSE 'Singular matrix'
				ENDIF
12			CONTINUE
			ENDIF
13		CONTINUE
		IPIV(ICOL)=IPIV(ICOL)+1
		IF (IROW.NE.ICOL) THEN
		DO 14 L=1,N
			DUM=A(IROW,L)
			A(IROW,L)=A(ICOL,L)
			A(ICOL,L)=DUM
14		CONTINUE
		DO 15 L=1,M
			DUM=B(IROW,L)
			B(IROW,L)=B(ICOL,L)
			B(ICOL,L)=DUM
15		CONTINUE
		ENDIF
		INDXR(I)=IROW
		INDXC(I)=ICOL
		IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
		PIVINV=1./A(ICOL,ICOL)
		A(ICOL,ICOL)=1
		DO 16 L=1,N
			A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
			B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
		IF(LL.NE.ICOL)THEN
			DUM=A(LL,ICOL)
			A(LL,ICOL)=0.
			DO 18 L=1,N
				A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
				B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
		ENDIF
21      CONTINUE
22	CONTINUE
	DO 24 L=N,1,-1
		IF(INDXR(L).NE.INDXC(L))THEN
			DO 23 K=1,N
				DUM=A(K,INDXR(L))
				A(K,INDXR(L))=A(K,INDXC(L))
				A(K,INDXC(L))=DUM
23			CONTINUE
		ENDIF
24	CONTINUE
	RETURN
	END
!-----------------------------------------------------------------
	subroutine coolcheck
! evaluates cooling assuming one kick
! scool is the value of s at which to evaluate e beam parameters
	include 'ecoolm2h4_c.f90'
! length of cooled section
	dscool = slencool
	scool=0
! beta functions at s=0
	betaxe0 = ebeamx*gamma0/sigurestx
	betaye0 = ebeamy*gamma0/siguresty
! rms unnormalized electron emittance
	epsex = ebeamx*sigurestx/gamma0
	epsey = ebeamy*siguresty/gamma0
! CS functions where cooling is applied

	betaxe = betaxe0 + scool**2 /betaxe0
	betaye = betaye0 + scool**2 /betaye0

	alfaxe = -scool/betaxe0
	alfaye = -scool/betaye0
! constant dispersion in cooling section for electrons and hadrons
	sigx2 = betaxe*epsey + (sigurestt*dispxe)**2
	sigy2 = betaye*epsey
	sigt = ebeams/beta/2.998e8
	sigt2 = sigt**2
! cool with a single bunch of length ebeams
! get the coefficients of the v/! in rest frame distribution
	sigma12 = epsex *gamma0*gamma0*(1 + (alfaxe*dispxe*sigurestt)**2 / sigx2)/betaxe
	sigma22 = epsey*gamma0*gamma0/betaye
	sigma32 = sigurestt**2 *epsex*betaxe/sigx2
	rho = alfaxe*dispxe*sigurestt/sqrt(sigx2 + (alfaxe*dispxe*sigurestt)**2)
! uzbar= coeffbaruz*x, uxbar = x*coeffbarux, uybar = y*coeffbaruy
	coeffbaruz = dispxe*sigurestt**2 / sigx2
	coeffbaruy = -gamma0*alfaye/betaye
	coeffbarux = -gamma0*alfaxe*epsex/(betaxe*epsex+(dispxe*sigurestt)**2)
	write(6,*)' kD_i =',coeffbaruz*dispxh
	write(6,*)'D_i uoff lamda = ',dispxh*uoff*xoff/sigx2
! coefficients of the distribution
	coeffa = 1/(sigma12*(1-rho*rho))
	coeffb = 1/(sigma32*(1-rho*rho))
	coeffd = 1/sigma22
	coeffc = rho/(sqrt(sigma12*sigma32)*(1-rho**2))
! following look OK
!      write(6,*)' a thru d', coeffa,coeffb,coeffc,coeffd,scool,alfaxe,rho,sigurestt
!       write(6,*)'coeffbars',coeffbaruz,coeffbaruy,coeffbarux
! get vectors (v11,v12) and (v21,v22) dots with (ux-uxbar,uz-uzbar)
	call getvs(v11,v12,v21,v22,coeffa,coeffb,coeffc)
! below look OK
!      write(6,*)'vs',v11,v12,v21,v22
! now need to find the smallest of the 3 sigmas. 
! sig1 is x, sig2 is z sig 3 is y
	sig1 = 1/sqrt(v11**2 + v12**2)
	sig2 = 1/sqrt(v21**2 + v22**2)
	sig3 = sqrt(sigma22)
! below look OK
!	write(6,*)' sigs',sig1,sig2,sig3
	small=sig1  
	ns=1
	if(sig2.lt.small)then
		small=sig2
		ns=2
	endif

	if(sig3.lt.small)then
		small=sig3
		ns=3
	endif
! next looks OK
!      write(6,*)' ns small',ns,small

! ns is the index of the small coordinate, fill in the 3 arrays with the coefficients
	if(ns.eq.1)then
		nsum1=0
		caray1(0)= 1
	else
		call getcaray(caray1,sig1,small,nsum1)            
	endif

	if(ns.eq.2)then
		nsum2=0
		caray2(0)= 1
	else
		call getcaray(caray2,sig2,small,nsum2)            
	endif

	if(ns.eq.3)then
		nsum3=0
		caray3(0)= 1
	else
		call getcaray(caray3,sig3,small,nsum3)            
	endif
	!write(6,*)'nsums',nsum1,nsum2,nsum3
	!write(6,*)
	!write(6,*)(caray1(m),m=-nsum1,nsum1)
	!write(6,*)
	!write(6,*)(caray2(m),m=-nsum2,nsum2)
	!write(6,*)
	!write(6,*)(caray3(m),m=-nsum3,nsum3)
	!write(6,*)

! OK we now have the 3 arrays
! get the coefficient on page 13 on may 2019 notes
	re = 2.818e-15
! coloumb log = clog. just calculate at s=0
	bmax = 2*sqrt(sigy2)
! electron number density in rest frame
	ebeamrest = qbeamcool/(1.602e-19*(2*pi)**1.5 * sqrt(sigy2*sigx2) *ebeams *gamma0)

	ebeamrest = abs(ebeamrest)
!      write(6,*) 'ebeamrest = ',ebeamrest
! transverse rms of v/! for hadrons in rest frame
	rmsux2 = (rmsx*gamma0)**2/(betaxcool*betavgx)
	rmsuy2 = (rmsy*gamma0)**2/(betaycool*betavgy)
	rmsrelbeta2 = sigurestt**2+ sigurestx**2 + siguresty**2 + rmsux2+rmsuy2 + dponp**2
	debye = rmsrelbeta2/(8*pi*2.818e-15*ebeamrest)
	debye = sqrt(debye)
!      write(6,*)' debye = ',debye
	if(debye.lt.bmax)bmax=debye
! for short collisions can have smaller bmax
	bmax1 = sqrt(rmsrelbeta2)*slencool/(beta*gamma0)

	if(bmax1.lt.bmax)bmax=bmax1
	bmin = 2*qatom*re/rmsrelbeta2
	clog = log(bmax/bmin)
! max number density in lab frame
	denlabmax = qbeamcool/(1.602e-19*sqrt(sigx2*sigy2*(2*pi)**3)*ebeams)
	ratmass=1/(1836.*aatom)
	ducoeff = dscool*4*sqrt(2*pi)*qatom**2*re**2*ratmass*denlabmax*clog/(sig1*sig2*sig3*gamma0**2)
	write(6,*)'du/u = ',ducoeff*speedfact/3
	sphere = ducoeff*speedfact/3
! unew = uold - ducoeff*sum_k C_k ... on page 13
! unit vector corresponding to sig1 multiplied by small sigma          
	h11 = v11*sig1*small
	h12 = v12*sig1*small
! unit vector corresponding to sig2           
	h21 = v21*sig2*small
	h22 = v22*sig2*small
	!write(6,*)'h s',h11,h12,h21,h22
	sqpid2 = sqrt(pi/2.)
	sqhalf = sqrt(0.5)
!234567
! apply speedup factor
	ducoeff = ducoeff*speedfact

!	write(6,*)'nsums',nsum1,nsum2,nsum3
! don't forget correction for spatial density

	sigxold=0
	sigxnew=0

	sigyold=0
	sigynew=0

	sigzold=0
	sigznew=0

! need to use physical variables, 
! x will be position in meters, not betatron coordinates
! px = v/! in rest frame
	cx = sqrt(betaxcool/betavgx)
! px with be vx/! in rest frame
	cux = gamma0/sqrt(betavgx*betaxcool)
	cy = sqrt(betaycool/betavgy)
	cuy = gamma0/sqrt(betavgy*betaycool)
! get physical coordinates at center of cooling section
	do k=1,np
		xk = cx*x(k)+  dispxh*pt(k)/gamma0
		pxk = cux*px(k)
		yk = cy*y(k)
		pyk = cuy*py(k)
		tk = t(k)
		pzk = pt(k)/gamma0
! need old and new emittance
		sigxold = sigxold + betaxcool*(pxk/gamma0)**2 + (xk-dispxh*pzk)**2/betaxcool
		sigyold = sigyold + betaycool*(pyk/gamma0)**2 + (yk)**2/betaycool
		sigzold = sigzold + pzk**2
! get sum_k C_k ... on page 13
		sumx=0
		sumy=0
		sumz=0
! average electron velocitty
		uzbar= coeffbaruz*xk
		uxbar = xk*coeffbarux
		uybar = yk*coeffbaruy
! subtract average electron velocity
		pxk1 = pxk - uxbar
		pyk1 = pyk - uybar
		pzk1 = pzk - uzbar

		do n1 = -nsum1,nsum1
			do n2 = -nsum2,nsum2
				do n3 = -nsum3,nsum3
					cktot = caray1(n1)*caray2(n2)*caray3(n3)
					offx = n1*h11 + n2*h21
					offz = n1*h12 + n2*h22
					offy = n3*small
					uxarg = pxk1-offx
					uyarg = pyk1-offy
! include energy offset 
					uzarg = pzk1-offz-uoff
! get I(R)/R**3
					r = sqrt(uxarg**2 + uyarg**2 +uzarg**2)
					rr = r/small
					fdivr3=0
					if(rr.lt.0.01)then
						fidivr3 = 1./3. - 0.1*rr*rr
					else
						if(rr.gt.6)then
							fidivr3 = sqpid2/rr**3
						else
							fidivr3 = (sqpid2*erf(rr*sqhalf)-rr*exp(-0.5*rr*rr))/rr**3
						endif
					endif
					!commont on 2022/3/11, file is too large 
					!write(24,101)rr,fidivr3,1.*n1,1.*n2,1.*n3,offx,offy,offz
					fidivr3 = fidivr3*cktot
					sumx = sumx + uxarg*fidivr3
					sumy = sumy + uyarg*fidivr3
					sumz = sumz + uzarg*fidivr3
				enddo
			enddo
		enddo
! include spatial offset
		xkoff = xk-xoff
		spacex = exp(-0.5*(xkoff*xkoff/sigx2 + yk*yk/sigy2))
		spacet = exp(-0.5*(tk*tk/sigt2))
         
! spatial density
		space = spacex*spacet
		dus = ducoeff*space
		pxk=pxk -dus*sumx
		pyk=pyk -dus*sumy
		pzk=pzk -dus*sumz

		dpx1 = sphere*pxk*space/(1+ (pxk*pxk+pyk*pyk+pzk*pzk)/(sig1**2 + sig2**2 + sig3**2))**1.5
		write(25,101)float(k),sumx,sumy,sumz,pxk,pyk,pzk,dus,xk,yk,space,dpx1

		sigxnew = sigxnew + betaxcool*(pxk/gamma0)**2 + (xk-dispxh*pzk)**2/betaxcool
		sigynew = sigynew + betaycool*(pyk/gamma0)**2 + (yk)**2/betaycool
		sigznew = sigznew + pzk**2
	enddo      
! ok get the change in the rms
	sigxold = sigxold/np
	sigyold = sigyold/np
	sigzold = sigzold/np

	sigxnew = sigxnew/np
	sigynew = sigynew/np
	sigznew = sigznew/np

	!write(6,*)'depsx/epsx epsx ', (sigxnew/sigxold-1), sigxold*0.5
	!write(6,*)'depsy/epsy epsy ', (sigynew/sigyold-1), sigyold*0.5
	!write(6,*)'dsigz/sigz sigz ', 0.5*(sigznew/sigzold-1), sqrt(sigzold)
	!write(6,*)' clog = ',clog

 101  format(20e12.3)
 102  format(2e13.3,4i5)
	return
	end

	
   subroutine getrate
!least squares mothed to get the cooling rate
!unit: rate1 (mm.mrad/s), rate2 (m/s)
		include 'ecoolm2h4_c.f90'
		real(8) rate1,rate2,points
		real(8) emitt,length,dpp
		open(unit=13,name='csmon.out',status='unknown')
		points = real(nturns/nwrite0)
		sumt = 0.0
		sumemit  = 0.0
		sumlong = 0.0
		sumdpp = 0.0

		sumt2   = 0.0
		sumtemit = 0.0
		sumtlong = 0.0
		sumtdpp = 0.0

		do i=1,points
			read(13,*) turns,csx,csy,css,csp
			!write(*,*) turns,csx,csy
			time = turns*trev*speedfact
			emittx = csx
			emitty = csy
			length = css/trev*circ
			dpp	= csp
			
			sumt = sumt + time
			sumemitx = sumemitx + emittx
			sumemity = sumemity + emitty
			sumlong = sumlong + length
			sumdpp = sumdpp + dpp
			
			sumt2 = sumt2 + time**2
			sumtemitx = sumtemitx + time * emittx
			sumtemity = sumtemity + time * emitty
			sumtlong = sumtlong + time * length
			sumtdpp = sumtdpp +time * dpp
		end do
	    	rate1 = (sumtemitx-((sumt*sumemitx)/points))/(sumt2-(sumt**2/points))
			rate2 = (sumtemity-((sumt*sumemity)/points))/(sumt2-(sumt**2/points))
	    	rate3 = (sumtlong-((sumt*sumlong)/points))/(sumt2-(sumt**2/points))
			rate4 = (sumtdpp-((sumt*sumdpp)/points))/(sumt2-(sumt**2/points))
		close(13)

		open(unit=13,name='rates.out',status='unknown')
		write(13,105) rate1,rate2,rate3,rate4
		close(13)
105 format(4e15.5)		
	end