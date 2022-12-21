 
    
SUBROUTINE EBEAM_DISTRIBUTION()
!Generate the initial e-beam distribution uniform
include 'ecoolm2h4_c.f90'
!################# E-Beam distribution in transverse ####################
    ple=tie/e0/clight/beta/gamma0/pi/(aeb*0.01)**2  !unit:1/m^2
    a_drift=ple/Bg/2./epsilon*e0
    a_sp_ch=2.818e-15*ple*3.14159/beta**2*clight*beta/gamma0**3
    OPEN(unit=1,file='eledis.dat')
    s=0.
    ss=0.
    DO i=1,300
        READ(1,*) eldis(i)
        s=s+eldis(i)*i
        sldis(i)=s              !sldis=[1,1+2,1+2+3,...,1+2+3+...+100]
        ss=ss+s/i
        ssdis(i)=ss             !ssdis=[ss(n)=ss(n-1)+(n+1)/2  OR  ss(n)=(n^2+2n+4)/4]
    END DO
    CLOSE(unit=1)
    An=45150./sldis(300)
    dscalx=1
    dscaly=1
    !write(*,*) "cool:", ple,a_drift,a_sp_ch
    !##################### Eletron beam velpulsewidthocity distribution in transverse #########
    OPEN(unit=4,file='elbeam.dat')
    DO i=1,300
        eldis(i)=eldis(i)*An
        sldis(i)=2.*An/300.*aeb*0.01*a_drift*sldis(i)/i         ![0,r]*a_drift
        ssdis(i)=0.04/3.0*An/300.*(aeb*0.01)**2*a_sp_ch*ssdis(i)    ![0,r*r]*a_sp_ch
        y1=eldis(i)*ple
        y2=sldis(i)
        y3=ssdis(i)
        WRITE(4,99) 0.01*i*aeb/300,y1,y2,y3          !实际座标系
    END DO
    CLOSE(unit=4)  
99     format(4e15.5)
END
    
    
	
    
SUBROUTINE EBEAM_DISTRIBUTION_NEW(nechoice)
!Generate the initial e-beam distribution uniform/elliptical/cos/gaussian/hollow
include 'ecoolm2h4_c.f90'
    integer nechoice
	real(8) s
	real(8) sl(300),ss(300)
!################# E-Beam distribution in transverse ####################
    ple=tie/e0/clight/beta/gamma0  !unit:1/m^2, DC-beam
    a_drift=ple*e0/Bg/epsilon
    a_sp_ch=ple*e0/epsilon/beta**2/gamma0**3/0.511E6*beta*clight
!	hollowsig = 0.1	!cm
!	hollowaeb = 1.0	!cm

	DO i=1,300
		s = i*0.01*aeb/300.0
		sl(i)=s
		if(nedis.eq.1) then    !uniform
		    ss(i) = 1.0/pi/(aeb*0.01)**2
			
		elseif(nedis.eq.2) then    !Elliptical
			ss(i) = 3.0/2.0/pi/(aeb*0.01)**2 * (1.0-(s/(aeb*0.01))**2)**0.5
		elseif(nedis.eq.3) then     !Parabolic
			ss(i) = 2.0/pi/(aeb*0.01)**2 * (1.0-(s/(aeb*0.01))**2)
		elseif(nedis.eq.4) then     !cos-square
			ss(i) = 2.0*pi/(pi**2-4)/(aeb*0.01)**2 * cos(pi*s/2/(aeb*0.01))**2
		elseif(nedis.eq.5) then     !Gaussian
			! here sigma is equal to aeb/3
			sigmaer1 = aeb*0.01/3.0
			ss(i) = 1.0/2.0/pi/sigmaer1**2*exp(-s**2/2.0/sigmaer1**2)
		elseif(nedis.eq.6) then     !hollow
			! aeb here is the positon offset 
			sigmaer2 = hollowsig*0.01
			Ann = hollowaeb*0.01*sqrt(pi/2.0)*(1+erf(hollowaeb*0.01/sqrt(2.0)/sigmaer2))
			Ann = Ann + sigmaer2*exp(-(hollowaeb*0.01)**2/2.0/sigmaer2**2)
			Ann = sigmaer2/Ann
			ss(i) = Ann/2.0/pi/sigmaer2**2*exp(-(s-hollowaeb*0.01)**2/2.0/sigmaer2**2)
		endif
		
		
	END DO
	
	!##################### Eletron beam velpulsewidthocity distribution in transverse #########
	OPEN(unit=4,file='elbeam.dat')
	
	tempdrift = 0.0
	tempdeviation = 0.0
	DO i=1,299
		
		eldis(i) = ss(i)
		!tempdrift = tempdrift + (sl(i)+sl(i+1))*(ss(i)+ss(i+1))/4.0*(aeb*0.01)/300.0
		tempdrift = tempdrift + sl(i)*ss(i)*(aeb*0.01)/300.0
		sldis(i) = tempdrift/sl(i)*a_drift
		
		tempdeviation = tempdeviation + sldis(i)/a_drift*aeb*0.01/300.0
		ssdis(i) = tempdeviation*a_sp_ch
			
!		sldis(i) = sldis(i)*0.0
		
		y1=eldis(i)*ple
        y2=sldis(i)
        y3=ssdis(i) 
        WRITE(4,99) 0.01*i*aeb/300,y1,y2,y3          !实际座标系
	END DO
    CLOSE(unit=4)  
99     format(4e15.5)
	END SUBROUTINE
	
			