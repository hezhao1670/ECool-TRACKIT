! this section is used to calcualte magnitized cooling (Parkhomchuk)
! the code is based on the old version of alpha_cool 
! with uniform (hollow) e-beam by input file 
!
! with transverse SC kick only in pulsed e-beam
    
    subroutine ebeamkick2
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
! include xoff and xpoff
			xk = cx*x(k)+  dispxh*pt(k)/gamma0/beta/beta + xoff
			pxk = cux* (px(k) + xpoff)
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
			pt(k)=pt(k)/gamma0/beta**2      !dp/p
		enddo 
! there are nslice cooling kicks there are two updates 
		elldg = slencool/2/gamma0/nslice
! s coordinate with s=0 in center of IP
		sslice=-slencool/2
		do k=1,nslice
			call driftem(elldg)
!		write(6,*)'coolslice1',px(np/4),px(np/2)
! next routine performs a cooling kick
! now pt(k) is the dp/p in lab frame
			sslice = sslice + slencool/2/nslice
			call coolslice2(sslice)
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
			x(k)=(xk-dispxh*pt(k) - xoff)/cx
			y(k)=yk/cy
			px(k)=pxk/cux -xpoff
			py(k)=pyk/cuy
			pt(k)= beta**2*gamma0*pt(k)
		enddo 
! all done
    return
    end
    
  
    
    
subroutine coolslice2(scool)
include 'ecoolm2h4_c.f90'
! scool is the value of s at which to evaluate e beam parameters
! scool is the value of s at which to evaluate e beam parameters
    do k=1,np
! physical coordinates at of cooling section
		pxk = px(k)
		xk = x(k) 
		yk = y(k)
		pyk = py(k)
        ptk = pt(k)
        tk = t(k)
        phiparticle = 2.0*pi*tk/trev
        call pulse_cooling(xk,yk,pxk,pyk,ptk,phiparticle)
       ! write(*,*) tk,dcoolp
        
        px(k) = pxk + dcoolx
		x(k)  = xk 
		y(k)  = yk
		py(k) = pyk + dcooly
        pt(k) = ptk + dcoolp + dkicky/beta/clight
        tk = t(k)
		
!        write(*,*) ptk,dcoolp,dkicky
	end do
	!write(*,*) t(100),pt(100),dkicky/beta/clight
end

   
    
    
SUBROUTINE pulse_cooling(x0,z0,xp0,zp0,vy0,phi0)
include 'ecoolm2h4_c.f90'

    REAL(8) x0,z0,xp0,zp0,vy0,phi0
    REAL(8) beginpoint,endpoint,kple

    dcoolx=0.0
    dcoolz=0.0
    dcooly=0.0
    dkickx=0.0
	dkicky=0.0
	dkickz=0.0
    
    kple = 0.0
!!parameters for the cooling force calculation	
!############ pulsed electron beam #########################
    totalpulse=pulseN*pulseW+(pulseN-1)*gapW+2*pulseN*riseW
    If(totalpulse.gt.1.) then
        WRITE(*,*) "The pulse electron beam is too long!!!!!!!!!!"
        stop
    endif
            
    mm=1
    pkvalue1=pulseN-2.0*(mm-1)
    pkvalue2=(pulseN-2.0*(mm-1))*0.5
    pkvalue3=(pulseN-2.0*(mm-1)-1)*0.5                    
    !beginpoint=-2.0*pi*(pkvalue1*riseW+pkvalue2*pulseW+pkvalue3*gapW)+0.0
		
    !write(*,*) "in:1",vx/v0/gamma0,vz/v0/gamma0,vy/gamma0/gamma0
    !write(*,*) "in",dcoolx,dcoolz,dcooly
!pulsed electron beam cooling (square distribution in longitudinal)
	
    Do j=1,pulseN
        mm=j
        pkvalue1=pulseN-2.0*(mm-1)
        pkvalue2=(pulseN-2.0*(mm-1))*0.5
        pkvalue3=(pulseN-2.0*(mm-1)-1)*0.5
        !write(*,*) pkvalue1,pkvalue2,pkvalue3
        beginpoint=-2.0*pi*(pkvalue1*riseW+pkvalue2*pulseW+pkvalue3*gapW)+epsoff
        beginpoint=beginpoint-2.0*pi*kturns*nperturn*(1.0-(feshift+1.0/trev)*trev)
		
		!write(*,*) epsoff,beginpoint
        if(beginpoint.gt.pi) then
           beginpoint= -pi+mod(beginpoint,2.0*pi)
            !beginpoint= -beginpoint+2.0*mod(beginpoint,pi)
            !write(*,*) 'here'
        else if(beginpoint.lt.(-pi)) then
            beginpoint=pi+mod(Abs(beginpoint),-2.0*pi)
        end if
          
       ! write(*,*) beginpoint
                        
        endpoint=beginpoint+2.*pi*(2.0*riseW+pulseW)

        if(endpoint.gt.pi) then
            endpoint=-pi+mod(endpoint,2.0*pi)
        else if (beginpoint.lt.(-pi)) then
            endpoint=pi+mod(endpoint,-2.0*pi)
        end if      
        
        
         
!#################################################   
        !write(*,*) beginpoint,endpoint                 
        if (endpoint.gt.beginpoint) then
            if(phi0.gt.beginpoint) then
                if(phi0.lt.(beginpoint+2.0*pi*riseW)) then
                    kple=(phi0-beginpoint)/(2.0*pi*riseW)
                    dkicky=-nperturn*speedfact*qatom/aatom*vkick/beta**2/(gamma0*E00*1.e6)*clight*beta
					!write(*,*) phi0,dkicky/beta/clight
                else if (phi0.lt.(beginpoint+2.0*pi*(riseW+pulseW))) then
                    kple=1.0
                else if (phi0.lt.(beginpoint+2.*pi*(2.0*riseW+pulseW))) then
                    kple=(endpoint-phi0)/(2.0*pi*riseW)
                    dkicky=nperturn*speedfact*qatom/aatom*vkick/beta**2/(gamma0*E00*1.e6)*clight*beta
					
                end if
                
                CALL PARKHOMCHUK_FORCE(x0,z0,xp0,zp0,vy0,phi0,kple)
            end if                                  
            
            
        else if(endpoint.lt.beginpoint) then 
            if(phi0.gt.beginpoint) then
                if(phi0.lt.(beginpoint+2.0*pi*riseW)) then
                    kple=(phi0-beginpoint)/(2.0*pi*riseW)
                    dkicky=-nperturn*speedfact*qatom/aatom*vkick/beta**2/(gamma0*E00*10**6)*clight*beta     
                else if (phi0.lt.pi) then
                    kple=1.0
                end if
                CALL PARKHOMCHUK_FORCE(x0,z0,xp0,zp0,vy0,phi0,kple)
                
            else if(phi0.lt.endpoint) then
                if(phi0.gt.(endpoint-2.0*pi*riseW)) then
                    kple=(endpoint-phi0)/(2.0*pi*riseW)
                    dkicky=nperturn*speedfact*qatom/aatom*vkick/beta**2/(gamma0*E00*10**6)*clight*beta
                else if(phi0.gt.-pi) then
                    kple=1.0
                end if 
                CALL PARKHOMCHUK_FORCE(x0,z0,xp0,zp0,vy0,phi0,kple)
            end if
        end if
	end do
    END
    
   
    
    
    
    
SUBROUTINE PARKHOMCHUK_FORCE(x0,z0,xp0,zp0,vy0,phi0,kple)
include 'ecoolm2h4_c.f90'
REAL(8) x0,z0,xp0,zp0,vy0,phi0,kple
! xp0,zp0,vy0 are the parameters in lab frame
! without space charge kicks 
    dscool = slencool/nslice
!parameters for the cooling force calculation
    etacool=dscool/circ
    glonge=1.+temel/0.511E6
    vlonge=clight*sqrt(1.-1./glonge/glonge)
    Veff=sqrt((beta*gamma0*clight*Berror)**2+vlonge**2)
    gtranse=1.+temet/0.511E6
    tett=Sqrt(1.-1./(1.+temet/0.511e6)**2)/beta
     
	
    X_real = x0
    r_xz2=X_real*X_real + z0*z0
    r_xz = sqrt(r_xz2)           
         
! jel: locate the particle's position inside the e-beam
    acu=0.0
    jel=299.*r_xz/(aeb*0.01)+1.0
    !write(*,*) jel
    if(jel.lt.1) jel=1
    if(jel.gt.300) goto 33
    acu=eldis(jel)*kple     
    if(acu.eq.0.0) goto 33
	
	
        
    !ple=tie/e0/clight/beta/gamma0/pi/(aeb*0.01)**2*kple !density
	ple=tie/e0/clight/beta/gamma0*eldis(jel)*kple   !real density
	
	!write(*,*) phi0,ple
	
	
    
    dek=4.*(2.8179E-15*qatom)**2*ple*clight**4/1836./aatom*etacool  !coefficient in the friction force
    tauf=dscool/clight/beta/gamma0      !flight of time in PRF
    armin=2.8179E-15*qatom*clight**2    !rho_min
    we=clight*sqrt(4.*pi*ple*2.8179E-15)    !plasma frequency
    
    rlar=tett*beta*gamma0*0.511E+06/clight/Bg   !rho_larmor
    dec=nperturn*trev*dek/gamma0*speedfact  !one turn period in PRF
    
    vx=xp0*gamma0*clight*beta/gamma0   !already in PRF frame
    vz=zp0*gamma0*clight*beta/gamma0   
    !vy=vy0*gamma0*gamma0
	vy = vy0*clight*beta

	 

    
    Vss=uoff*beta*clight
    v_shift = ssdis(jel)*kple      
    Vsd = Vss + v_shift     !Vss is the set energy tuning 
    !Vxd = sldis(jel)*kple*(-z0)/r_xz
    !Vzd = sldis(jel)*kple*X_real/r_xz
	
	Vxd = sldis(jel)*kple*X_real/r_xz
    Vzd = sldis(jel)*kple*z0/r_xz
	
    v2= (vx-Vxd)**2 + (vz-Vzd)**2+(vy-Vsd)**2+Veff**2
    v=sqrt(v2)
!    write(*,*) Vsd,Vxd,Vzd
    
    rmax=v/(we+1.0/tauf)
    rmin=armin/v2
    Cln=log((rmax+rlar+rmin)/(rlar+rmin))
	
	!write(*,*) 'rmax','rlar','rmin'
	!write(*,*) rmax,',',rlar,',',rmin
    !write(*,*) Cln
	
	
	
    !dd=exp(-dec*acu*Cln/v2/v)   !cooling step:  exp(-lambda*delta_t)
    dd=exp(-dec*Cln/v2/v)   !cooling step:  exp(-lambda*delta_t)
	
	
    
    vx=(vx-Vxd)*dd+Vxd
    vz=(vz-Vzd)*dd+Vzd
    vy=(vy-Vsd)*dd+Vsd

	 
	!space charge transverse (need to be confirmed)
    re = 2.818e-15  !electron
	rp = re/1836.	!proton
	rion= e0**2/4.0/pi/epsilon/a_mion/clight**2

	!dvx=Abs(x0)*qatom/aatom*2.0*pi*ple*rion*dscool*clight/beta/gamma0
	!dvz=Abs(z0)*qatom/aatom*2.0*pi*ple*rion*dscool*clight/beta/gamma0
    !dvx=x0*qatom/aatom*2.0*pi*ple*rion*dscool/beta**2/gamma0**3 !delta in lab frame delta_x'
	!dvz=z0*qatom/aatom*2.0*pi*ple*rion*dscool/beta**2/gamma0**3 !delta in lab frame delta_y'
	
	! 考虑任意电子束横向分布
	dvx = x0/r_xz*qatom/aatom*4.0*pi*rion*dscool/beta**2/gamma0**3*sldis(jel) !delta in lab frame delta_x'
	dvz	= z0/r_xz*qatom/aatom*4.0*pi*rion*dscool/beta**2/gamma0**3*sldis(jel) !delta in lab frame delta_y'
	!write(*,*) (x0**2+z0**2)**0.5,(dvx**2+dvz**2)**0.5
	
    !cooling effect in PRF
    dcoolx=(vx/clight/beta-xp0)*fracecool - dvx*gamma0*0.0
    dcooly=(vz/clight/beta-zp0)*fracecool - dvz*gamma0*0.0
    dcoolp=(vy/clight/beta-vy0)*fracecool   
    
    !write(*,*) phi0,dvx
    !write(*,*) v,dd
	
	!write(*,1121) vx/clight/beta,dcoolx,vy/clight/beta,dcoolp
1121 format(4e15.5)
	
33	continue
end