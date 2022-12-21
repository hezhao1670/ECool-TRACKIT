
!------------------------------------------------------------------
	subroutine getrates(alfx,alfy,alfp,bx,bxp,by,byp,dxin,dxp,dy,dyp,ibstype)
	!parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19)
	!common/parms/pnum,qatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb,r0
	!include 'trackit8az_c.f90'
	include 'ecoolm2h4_c.f90'
! alfx,alfy,alfp are emittance growth rates for bunched beams in sec^{-1}
! first get some other parameters
	gamma = gamma0
	beta = dsqrt(1-1/gamma**2)
	dx = dxin

! vertical dispersion is zero
!   write(6,*) '  gamma*sig(thetax) =',dsqrt(epsx/bx)*gamma
	sigx = dsqrt(bx*epsx + dponp**2 * dx**2 )
	sigy = dsqrt(by*epsy)
!      write(6,*)'sigx,sigy',sigx,sigy
!      write(6,*)'95% normalized emittances ',sigx*sigx*6*gamma/bx,
!     : sigy*sigy*6*gamma/by
!23456
	sigurest2 = dponp**2 + gamma*gamma*(epsx/bx +epsy/by)
! classical radius
	
	r0 = (qatom**2/aatom)*1.54e-18
	rhomin = r0/sigurest2
	rhomax = 2*sigy

	rhomax2 = dsqrt(sigurest2*sigx*sigy*sigs/(pnum*r0))
	if(rhomax2.lt.rhomax)rhomax=rhomax2

	coulomb = log(rhomax/rhomin)
!      write(6,*)' coloumb log = ',coloumb,rhomax/rhomin
! Nagaitsev's coefficient
	coeffn = pnum*r0*r0*clight*coulomb/(12*pi*beta**3 * gamma**5 * sigs)
	
! now need the other parms

	phi = dxp - 0.5*bxp*dx/bx
! three kinds of ibs
! ibstype=1, keep full algorithm
	if(ibstype.eq.2)phi=0.
	if(ibstype.eq.3)then
	dx2 = dx*dx + (bx*dpx-0.5*bxp*dx)**2
	dx = dsqrt(dx2)
	phi=0
	endif

	ax = bx/epsx
	ay = by/epsy
	as = ax*(dx**2/bx**2 + phi**2)+1/dponp**2
	a1 = 0.5*(ax+as*gamma**2)
	a2 = 0.5*(ax-as*gamma**2)
! eigenvalues of eq 2 in Nagaitsev 2005
	flam1 = ay
	disc = dsqrt(a2*a2 + (gamma*ax*phi)**2)
	flam2 = a1+disc
	flam3 = a1-disc
! eq 25 thru 27
	arg1 = 1/flam1
	arg2 = 1/flam2
	arg3 = 1/flam3
	r1 = rd(arg2,arg3,arg1,ierr1)*arg1
	r2 = rd(arg3,arg1,arg2,ierr2)*arg2
	r3 = rd(arg1,arg2,arg3,ierr3)*arg3
!write(*,*) arg1,arg2,arg3
	if(ierr1+ierr2+ierr3.gt.0)write(6,*)'trouble  ',ierr1,ierr2,ierr3


! eq 33 thru 35
	sp = 0.5*gamma**2 *(2*r1 - r2*(1-3*a2/disc) -r3*(1+3*a2/disc))
	sx = 0.5*(2*r1 - r2*(1+3*a2/disc) -r3*(1-3*a2/disc))
	sxp = 3*(gamma*phi)**2 * ax *(r3-r2)/disc
! eq 30 thru 32 but normalized to emittance to give growth rates
	alfp = coeffn*sp/(sigx*sigy*dponp**2)
	alfy = coeffn*by*(r3+r2-2*r1)/(sigx*sigy*epsy)
	alfx = coeffn*bx*(sx + (dx**2/bx**2 + phi**2)*sp + sxp)/(sigx*sigy*epsx)
	return
	end




!============================================================================
!=============   Borjken-Mtingua Model with vertical dispersion ============
! Same model with BETACOOL
	
	subroutine getrates2(alfx,alfy,alfp,bx,bxp,by,byp,dxin,dxp,dy,dyp,ibstype)
	    
	!parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19)
	!common/parms/pnum,qatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb,r0
	!include 'trackit8az_c.f90'
	include 'ecoolm2h4_c.f90'
	!real(8)::a(2,2),b(2,2),c(2,2)
	real(8),dimension(3,3)::Lh,Lv,Lp,L
	real(8),dimension(3,3)::II,LL,LB,KK,TT,LLh,LLv,LLp
	real(8)::lamLimit,stepLam,lambda,det
	gamma = gamma0
	r0 = (qatom**2/aatom)*1.54e-18

	stepLam = 45000.
	lamLimit = 4.0
	beta = dsqrt(1-1/gamma**2)
	dx = dxin
!write(6,*) '  gamma*sig(thetax) =',dsqrt(epsx/bx)*gamma
	sigx = dsqrt(bx*epsx + dponp**2 * dx**2 )
	sigy = dsqrt(by*epsy + dponp**2 * dy**2)

!23456
	sigurest2 = dponp**2 + gamma*gamma*(epsx/bx +epsy/by)
! classical radius
	r0 = (qatom**2/aatom)*1.54e-18
	rhomin = r0/sigurest2
	rhomax = 2*sigy

	rhomax2 = dsqrt(sigurest2*sigx*sigy*sigs/(pnum*r0))
	if(rhomax2.lt.rhomax)rhomax=rhomax2

	coulomb = log(rhomax/rhomin)
!coefficient
	coeffn = clight*r0**2*pnum*coulomb/(8.0*pi*beta**3*gamma**4*epsx*epsy*sigs*dponp)
    
	phix = dxp - 0.5*bxp*dx/bx
	phiy = dyp - 0.5*byp*dy/by
	Hx = (dx*dx + bx*bx*phix*phix)/bx
	Hy = (dy*dy + by*by*phiy*phiy)/by    
    
	Lh(1,1) = bx/epsx
	Lh(1,2) = -gamma*bx*phix/epsx
	Lh(2,1) = Lh(1,2)
	Lh(2,2) = gamma*gamma*Hx/epsx
	Lh(3,3) = 0.

	Lv(1,1) = 0.
	Lv(2,2) = gamma*gamma*Hy/epsy
	Lv(2,3) = -gamma*by*phiy/epsy
	Lv(3,2) = Lv(2,3)
	Lv(3,3) = by/epsy

	Lp(1,1) = 0.
	Lp(2,2) = gamma*gamma/dponp**2
	Lp(3,3) = 0.
    
	uplimit = Lh(1,1)
	do i=1,3
		do j=1,3
			L(i,j) = Lp(i,j) + Lh(i,j) + Lv(i,j)
			if(uplimit.lt.abs(L(i,j))) uplimit=abs(L(i,j))
			!The integral limit is the biggest one in the matrix L

			KK(i,j) = 0.
			if(i.eq.j) Then
			II(i,j) = 1.0
			else
			II(i,j) = 0.0
			end if             
		end do
	end do
	uplimit = uplimit * lamLimit !LamLimit is a given number at the initial
    
	do k=0,stepLam
    	lambda = uplimit * float(k) / stepLam  !stepLam is the total integral steps
    	LL = L + II*lambda
		
!det = M33DET(LL)
    det = LL(1,1)*LL(2,2)*LL(3,3)- LL(1,1)*LL(2,3)*LL(3,2)- LL(1,2)*LL(2,1)*LL(3,3)+ LL(1,2)*LL(2,3)*LL(3,1)+ LL(1,3)*LL(2,1)*LL(3,2) - LL(1,3)*LL(2,2)*LL(3,1)
!det = det + LL(1,2)*LL(2,3)*LL(3,1)+ LL(1,3)*LL(2,1)*LL(3,2) - LL(1,3)*LL(2,2)*LL(3,1)
!write(*,*) LL(1,1)

    	call M33INV (LL, LB, OK_FLAG) 
        
		!Tr = LB.Trace()
        
		do i=1,3
			do j=1,3
    			Tr0 = LB(1,1) + LB(2,2) + LB(3,3)
				if(i.eq.j) then
				    ame = 1.0
				else
				    ame = 0.0
				end if

				KK(i,j) = KK(i,j) + dsqrt(lambda/det)*(Tr0 * ame - 3.0*LB(i,j))
			end do
		end do
        
	end do
    
	do i=1,3
		do j=1,3
		    KK(i,j) = KK(i,j) * uplimit / stepLam
		end do
	end do
    
	alfx = 0.0
	alfy = 0.0
	alfp = 0.0
	do i=1,3
		do j=1,3
			alfx = alfx + coeffn * (KK(i,j) * Lh(i,j))
			alfy = alfy + coeffn * (KK(i,j) * Lv(i,j))
			alfp = alfp + coeffn * (KK(i,j) * Lp(i,j))    
		end do
	end do

	
   !write(*,*) bx,alfx,alfy,alfp    
        
	!do k=0,stepLam
 !       lambda = uplimit * real(k) / stepLam  !stepLam is the total integral steps
	!	!write(*,*) uplimit,stepLam
 !       LL = L + II*lambda
	!	
 !       !det = M33DET(LL)
 !       det = LL(1,1)*LL(2,2)*LL(3,3)- LL(1,1)*LL(2,3)*LL(3,2)- LL(1,2)*LL(2,1)*LL(3,3)+ LL(1,2)*LL(2,3)*LL(3,1)+ LL(1,3)*LL(2,1)*LL(3,2) - LL(1,3)*LL(2,2)*LL(3,1)
	!	!write(*,*) LL(1,1)
	!	
 !       call M33INV (LL, LB, OK_FLAG) 
 !       
 !       !Tr = LB.Trace()
 !       Tr0 = LB(1,1) + LB(2,2) + LB(3,3)
 !       Trh = Lh(1,1) + Lh(2,2) + Lh(3,3)
 !       Trv = Lv(1,1) + Lv(2,2) + Lv(3,3)
 !       Trp = Lp(1,1) + Lp(2,2) + Lp(3,3)
 !       
 !       LLh = matmul(Lh,LB)
 !       LLv = matmul(Lv,LB)
 !       LLp = matmul(Lp,LB)
	!	
 !       Trh1 = LLh(1,1) + LLh(2,2) + LLh(3,3)
 !       Trv1 = LLv(1,1) + LLv(2,2) + LLv(3,3)
 !       Trp1 = LLp(1,1) + LLp(2,2) + LLp(3,3)
 !       
 !       alfx = coeffn * dsqrt(lambda/det) * (Trh * Tr0 - 3. * Trh1) * uplimit/stepLam
 !       alfy = coeffn * dsqrt(lambda/det) * (Trv * Tr0 - 3. * Trv1) * uplimit/stepLam
 !       alfp = coeffn * dsqrt(lambda/det) * (Trp * Tr0 - 3. * Trp1) * uplimit/stepLam 
	!	
 !       write(*,*) alfx,alfy,alfp
	!end do
    
	end
    

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!=============   S. Naigaitsev with vertical dispersion ============
!With simplified integral for speed

	subroutine getrates3(alfx,alfy,alfp,bx,bxp,by,byp,dxin,dxp,dy,dyp,ibstype)

	!parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19)
	!common/parms/pnum,qatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb,r0
	!include 'trackit8az_c.f90'
	include 'ecoolm2h4_c.f90'
	real(8),dimension(3,3):: BBx,BBy,BBp,GG,GGi,LL,LLi,Si,Si1
	real(8),dimension(3,3):: RR,psi,MM,MM1,MM2,MM3
	real(8),dimension(3):: DD
	real(8) rx,ry,rp
	
	gamma = gamma0
	r0 = (qatom**2/aatom)*1.54e-18

	rx = 0.0
	ry = 0.0
	rp = 0.0
	beta = dsqrt(1.-1./gamma**2)
	dx = dxin
!write(6,*) '  gamma*sig(thetax) =',dsqrt(epsx/bx)*gamma
	sigx = dsqrt(bx*epsx + dponp**2 * dx**2 )
	sigy = dsqrt(by*epsy + dponp**2 * dy**2)
    
!23456
	sigurest2 = dponp**2 + gamma*gamma*(epsx/bx +epsy/by)
    
!coefficient
	phix = dxp - 0.5*bxp*dx/bx
	phiy = dyp - 0.5*byp*dy/by
    
	Ax = dx * dx / bx + phix*phix * bx
	Ay = dy * dy / by + phiy*phiy * by
	ff1 = dx*dx*dponp*dponp / (bx * epsx)
	ff2 = dy*dy*dponp*dponp / (by * epsy)
	ff = 1. + ff1 + ff2
	
	BBx(1,1) = bx
	BBx(1,3) = -phix * bx
	BBx(3,1) = BBx(1,3)
	BBx(3,3) = Ax
	
	BBy(1,1) = 0.
	BBy(2,2) = by
	BBy(2,3) = -phiy * by
	BBy(3,2) = BBy(2,3)
	BBy(3,3) = Ay
	
	BBp(1,1) = 0.0
	BBp(2,2) = 0.0
	BBp(3,3) = sigs / dponp
	
	GG(1,1) = 1.
	GG(2,2) = 1.
	GG(3,3) = 1./gamma
	
	LL(1,1) = bx / epsx
	LL(1,2) = 0.0
	LL(1,3) = -(phix * bx / epsx)
	LL(2,1) = 0.0
	LL(2,2) = by / epsy
	LL(2,3) = -(phiy * by / epsy)
	LL(3,1) = LL(1,3)
	LL(3,2) = LL(2,3)
	LL(3,3) = 1.0/dponp/dponp + Ax/epsx + Ay/epsy
	     

	call M33INV(LL,LLi,OK_FLAG)
	Si1 = matmul(LLi,GG)
	Si = matmul(transpose(GG),Si1)
	call jacobi_eigenvalue(3,Si,50000, RR, DD, it_num, rot_num)
!call Jacobi(Si,RR,1.E-19,3)
   	 
	qq = DD(1) + DD(2) + DD(3)
	dd1 = DD(1)
	dd2 = DD(2)
	dd3 = DD(3)
	!write(*,*) 'RERE',dd1,dd2,dd3
	r1 = RSD(dd2,dd3,dd1)
	r2 = RSD(dd3,dd1,dd2)
	r3 = RSD(dd1,dd2,dd3)
    

	psi1 = -2.*dd1*r1 + dd2*r2 + dd3*r3
	psi2 = -2.*dd2*r2 + dd1*r1 + dd3*r3
	psi3 = -2.*dd3*r3 + dd1*r1 + dd2*r2
	
	psi(1,1) = psi1
	psi(2,2) = psi2
	psi(3,3) = psi3
    
	call M33INV(GG,GGi,OK_FLAG)

	MM1 = matmul(transpose(RR),GGi)
	MM2 = matmul(psi,MM1)
	MM3 = matmul(RR,MM2)
	MM = matmul(transpose(GGi),MM3)


	coe = 1.0 / dsqrt(ff * bx * by)
	do i=1,3
		do j=1,3
			rx = rx + coe * BBx(i,j) * MM(i,j) 
			ry = ry + coe * BBy(i,j) * MM(i,j)
			rp = rp + coe * BBp(i,j) * MM(i,j)
		end do
	end do
  
  	
!Naigatsev coulomb log	
	rhomin = r0/qq/(gamma**2*beta**2)
	rhomax = sigx
	if (rhomax.gt.sigy) rhomax = sigy
	if (rhomax.gt.gamma*sigs) rhomax = gamma*sigs     
	coulomb = log(rhomax/rhomin)
    !write(*,*) coulomb

	coeffn = clight*r0**2*pnum/(12.0*pi*beta**3*gamma**5*dsqrt(epsx)*dsqrt(epsy)*sigs)
	alfx = coeffn * coulomb * rx / epsx
	alfy = coeffn * coulomb * ry / epsy
	alfp = coeffn * coulomb * rp / sigs / dponp
	 !write(*,*) "rr:", alfx,alfy,alfp
	end

	
	
!-------------------------------------------------------------------------------------------------
!same as the RD()
	REAL(8) FUNCTION RSD(XX,YY,ZZ)	
	IMPLICIT NONE
	REAL(8) C1,C2,C3,C4,C5,C6
	REAL(8) XX,YY,ZZ,XT,YT,ZT,DETX,DETY,DETZ,ERRTOL
	REAL(8) SUM,FAC,EPSLON,EA,EB,EC,ED,EE
	REAL(8) SQRTX,SQRTY,SQRTZ,ALAMB,AVE,RD1,RD2

	C1 = 3.0E0/14.0E0
	C2 = 1.0E0/6.0E0
	C3 = 9.0E0/22.0E0
	C4 = 3.0E0/26.0E0
	C5 = 9.0E0/88.0E0
	C6 = 9.0E0/52.0E0
	!XX=0.2
	!YY=0.3
	!ZZ=0.4
	ERRTOL = 0.05
	XT = XX
	YT = YY
	ZT = ZZ
	DETX = 1.
	DETY = 1.
	DETZ = 1.

	SUM = 0.0
	FAC = 1.0
	
55	SQRTX = dsqrt(XT)
	SQRTY = dsqrt(YT)
	SQRTZ = dsqrt(ZT) 

		
	ALAMB = SQRTX * (SQRTY + SQRTZ) +SQRTY*SQRTZ
	SUM = SUM + FAC / (SQRTZ * (ZT + ALAMB))
	FAC = 0.25 * FAC
	XT = 0.25 * (XT + ALAMB)
	YT = 0.25 * (YT + ALAMB)
	ZT = 0.25 * (ZT + ALAMB)
	AVE = 0.2 * (XT + YT + 3.0*ZT)
		
	DETX = (AVE	- XT)/AVE
	DETY = (AVE - YT)/AVE
	DETZ = (AVE - ZT)/AVE
	
	EPSLON = MAX(ABS(DETX),ABS(DETY),ABS(DETZ))
	IF(EPSLON.LT.ERRTOL) GO TO 66
	
	GO TO 55
!	
66	EA = DETX * DETY
	EB = DETZ * DETZ
	EC = EA - EB
	ED = EA - 6.0 * EB
	EE = ED + EC + EC
	
	RD1 = 1.0 + ED * (-C1 + C5*ED - C6*DETZ*EE)
	RD2 = C2*EE + DETZ*(-C3*EC + DETZ*C4*EA)
	RSD = 3.0 * SUM + FAC * (RD1 + DETZ*RD2) / (AVE * dsqrt(AVE))
	END	FUNCTION

    
!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

	SUBROUTINE M33INV (A, AINV, OK_FLAG)

	IMPLICIT NONE

	REAL(8), DIMENSION(3,3), INTENT(IN)  :: A
	REAL(8), DIMENSION(3,3), INTENT(OUT) :: AINV
	LOGICAL, INTENT(OUT) :: OK_FLAG

	REAL(8), PARAMETER :: EPS = 1.0D-10
	REAL(8):: DET
	REAL(8), DIMENSION(3,3) :: COFACTOR

	DET =   A(1,1)*A(2,2)*A(3,3)  &
	- A(1,1)*A(2,3)*A(3,2)  &
	- A(1,2)*A(2,1)*A(3,3)  &
	+ A(1,2)*A(2,3)*A(3,1)  &
	+ A(1,3)*A(2,1)*A(3,2)  &
	- A(1,3)*A(2,2)*A(3,1)

	IF (ABS(DET) .LE. EPS) THEN
		AINV = 0.0D0
		OK_FLAG = .FALSE.
		RETURN
	END IF

	COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
	COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
	COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
	COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
	COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
	COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
	COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
	COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
	COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

	AINV = TRANSPOSE(COFACTOR) / DET

	OK_FLAG = .TRUE.

	RETURN
	END SUBROUTINE M33INV



!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----  get the eigenvactors and eigenvalue of one matrix --------------------
	subroutine jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num )
!*****************************************************************************80
!
!! JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
!								    
!  Discussion:
!
!    This function computes the eigenvalues and eigenvectors of a
!    real symmetric matrix, using Rutishauser's modfications of the classical
!    Jacobi rotation method with threshold pivoting. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 September 2013
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A(N,N), the matrix, which must be square, real,
!    and symmetric.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, real ( kind = 8 ) V(N,N), the matrix of eigenvectors.
!
!    Output, real ( kind = 8 ) D(N), the eigenvalues, in descending order.
!
!    Output, integer ( kind = 4 ) IT_NUM, the total number of iterations.
!
!    Output, integer ( kind = 4 ) ROT_NUM, the total number of rotations.
!
	implicit none

	integer ( kind = 4 ) n

	real ( kind = 8 ) a(n,n)
	real ( kind = 8 ) bw(n)
	real ( kind = 8 ) c
	real ( kind = 8 ) d(n)
	real ( kind = 8 ) g
	real ( kind = 8 ) gapq
	real ( kind = 8 ) h
	integer ( kind = 4 ) i
	integer ( kind = 4 ) it_max
	integer ( kind = 4 ) it_num
	integer ( kind = 4 ) j
	integer ( kind = 4 ) k
	integer ( kind = 4 ) l
	integer ( kind = 4 ) m
	integer ( kind = 4 ) p
	integer ( kind = 4 ) q
	integer ( kind = 4 ) rot_num
	real ( kind = 8 ) s
	real ( kind = 8 ) t
	real ( kind = 8 ) tau
	real ( kind = 8 ) term
	real ( kind = 8 ) termp
	real ( kind = 8 ) termq
	real ( kind = 8 ) theta
	real ( kind = 8 ) thresh
	real ( kind = 8 ) v(n,n)
	real ( kind = 8 ) w(n)
	real ( kind = 8 ) zw(n)

	do j = 1, n
		do i = 1, n
	  		v(i,j) = 0.0D+00
		end do
		v(j,j) = 1.0D+00
	end do

	do i = 1, n
		d(i) = a(i,i)
	end do

	bw(1:n) = d(1:n)
	zw(1:n) = 0.0D+00
	it_num = 0
	rot_num = 0

	do while ( it_num < it_max )
    	it_num = it_num + 1
!
!  The convergence threshold is based on the size of the elements in
!  the strict upper triangle of the matrix.
!
		thresh = 0.0D+00
		do j = 1, n
			do i = 1, j - 1
				thresh = thresh + a(i,j) ** 2
			end do
		end do

    	thresh = dsqrt ( thresh ) / real ( 4 * n, kind = 8 )

		if ( thresh == 0.0D+00 ) then
			exit 
		end if

		do p = 1, n
			do q = p + 1, n
				gapq = 10.0D+00 * abs ( a(p,q) )
				termp = gapq + abs ( d(p) )
				termq = gapq + abs ( d(q) )
!
!  Annihilate tiny offdiagonal elements.
!
				if ( 4 < it_num .and. &
					termp == abs ( d(p) ) .and. &
					termq == abs ( d(q) ) ) then
					a(p,q) = 0.0D+00
!
!  Otherwise, apply a rotation.
        		else if ( thresh <= abs ( a(p,q) ) ) then

					h = d(q) - d(p)
					term = abs ( h ) + gapq

					if ( term == abs ( h ) ) then
		        		t = a(p,q) / h
					else
		        		theta = 0.5D+00 * h / a(p,q)
		        		t = 1.0D+00 / ( abs ( theta ) + dsqrt ( 1.0D+00 + theta * theta ) )
						if ( theta < 0.0D+00 ) then 
							t = - t
		        	end if
				end if

				c = 1.0D+00 / dsqrt ( 1.0D+00 + t * t )
				s = t * c
				tau = s / ( 1.0D+00 + c )
				h = t * a(p,q)
!
!  Accumulate corrections to diagonal elements.
				zw(p) = zw(p) - h                  
				zw(q) = zw(q) + h
				d(p) = d(p) - h
				d(q) = d(q) + h

				a(p,q) = 0.0D+00
!
!  Rotate, using information from the upper triangle of A only.
!
				do j = 1, p - 1
					g = a(j,p)
					h = a(j,q)
					a(j,p) = g - s * ( h + g * tau )
					a(j,q) = h + s * ( g - h * tau )
				end do

				do j = p + 1, q - 1
					g = a(p,j)
					h = a(j,q)
					a(p,j) = g - s * ( h + g * tau )
					a(j,q) = h + s * ( g - h * tau )
				end do

				do j = q + 1, n
					g = a(p,j)
					h = a(q,j)
					a(p,j) = g - s * ( h + g * tau )
					a(q,j) = h + s * ( g - h * tau )
				end do
!
!  Accumulate information in the eigenvector matrix.
!
				do j = 1, n
					g = v(j,p)
					h = v(j,q)
					v(j,p) = g - s * ( h + g * tau )
					v(j,q) = h + s * ( g - h * tau )
				end do

				rot_num = rot_num + 1

        		end if
  			end do
    	end do

		bw(1:n) = bw(1:n) + zw(1:n)
		d(1:n) = bw(1:n)
		zw(1:n) = 0.0D+00

		end do
!
!  Restore upper triangle of input matrix.
!
		do j = 1, n
			do i = 1, j - 1
				a(i,j) = a(j,i)
			end do
		end do
!
!  Ascending sort the eigenvalues and eigenvectors.
!
	  	do k = 1, n - 1

		m = k

		do l = k + 1, n
			if ( d(l) < d(m) ) then
			m = l
			end if
		end do

    	if ( m /= k ) then

			t    = d(m)
			d(m) = d(k)
			d(k) = t

			w(1:n)   = v(1:n,m)
			v(1:n,m) = v(1:n,k)
			v(1:n,k) = w(1:n)

    	end if

  	end do

  	return
    end




!*DECK RD
	REAL(8) FUNCTION RD (X, Y, Z, IER)
!***BEGIN PROLOGUE  RD
!***PURPOSE  Compute the incomplete or complete elliptic integral of the
!            2nd kind.  For X and Y nonnegative, X+Y and Z positive,
!             RD(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -3/2
!                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X or Y is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      SINGLE PRECISION (RD-S, DRD-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. !.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     RD
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the second kind
!          Standard FORTRAN function routine
!          Single precision version
!          The routine calculates an approximation result to
!          RD(X,Y,Z) = Integral from zero to infinity of
!                              -1/2     -1/2     -3/2
!                    (3/2)(t+X)    (t+Y)    (t+Z)    dt,
!          where X and Y are nonnegative, X + Y is positive, and Z is
!          positive.  If X or Y is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling Sequence
!
!          RD( X, Y, Z, IER )
!
!          Parameters on Entry
!          Values assigned by the calling routine
!
!          X      - Single precision, nonnegative variable
!
!          Y      - Single precision, nonnegative variable
!
!                   X + Y is positive
!
!          Z      - Real, positive variable
!
!
!
!          On Return     (values assigned by the RD routine)
!
!          RD     - Real approximation to the integral
!
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine.  It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!
!          X, Y, Z are unaltered.
!
!   3.    Error Messages
!
!         Value of IER assigned by the RD routine
!
!                  Value Assigned         Error Message Printed
!                  IER = 1                MIN(X,Y) .LT. 0.0E0
!                      = 2                MIN(X + Y, Z ) .LT. LOLIM
!                      = 3                MAX(X,Y,Z) .GT. UPLIM
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y, and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                    Not less  than 2 / (machine maximum) ** (2/3).
!
!          UPLIM  - Upper limit of valid arguments
!
!                    Not greater than (0.1E0 * ERRTOL / machine
!                    minimum) ** (2/3), where ERRTOL is described below.
!                    In the following table it is assumed that ERRTOL
!                    will never be chosen smaller than 1.0E-5.
!
!
!                    Acceptable Values For:   LOLIM      UPLIM
!                    IBM 360/370 SERIES   :   6.0E-51     1.0E+48
!                    CDC 6000/7000 SERIES :   5.0E-215    2.0E+191
!                    UNIVAC 1100 SERIES   :   1.0E-25     2.0E+21
!                    CRAY                 :   3.0E-1644   1.69E+1640
!                    VAX 11 SERIES        :   1.0E-25     4.5E+21
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!          ERRTOL    Relative error due to truncation is less than
!                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
!
!
!
!              The accuracy of the computed approximation to the inte-
!              gral can be controlled by choosing the value of ERRTOL.
!              Truncation of a Taylor series after terms of fifth order
!              introduces an error less than the amount shown in the
!              second column of the following table for each value of
!              ERRTOL in the first column.  In addition to the trunca-
!              tion error there will be round-off error, but in prac-
!              tice the total error from both sources is usually less
!              than the amount given in the table.
!
!
!
!
!          Sample Choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0E-3    4.0E-18
!                           3.0E-3    3.0E-15
!                           1.0E-2    4.0E-12
!                           3.0E-2    3.0E-9
!                           1.0E-1    4.0E-6
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   RD Special Comments
!
!
!
!          Check: RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y)
!          = 3 /  dsqrt(X * Y * Z), where X, Y, and Z are positive.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral RD(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, and Z are unaltered.
!
!
!
!          ********************************************************
!
!           WARNING: Changes in the program may improve speed at the
!                    expense of robustness.
!
!
!
!    -------------------------------------------------------------------
!
!
!   Special Functions via RD and RF
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
!                  ----------------------------------------------
!
!
!                                            2         2   2
!                  E(PHI,K) = SIN(PHI) RF(COS (PHI),1-K SIN (PHI),1) -
!
!                     2      3            2         2   2
!                  -(K/3) SIN (PHI) RD(COS (PHI),1-K SIN (PHI),1)
!
!
!                                 2        2           2
!                  E(K) = RF(0,1-K ,1) - (K/3) RD(0,1-K ,1)
!
!
!                         PI/2     2   2      1/2
!                       = INT  (1-K SIN (PHI) )  D PHI
!                          0
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
!                  ----------------------------------------------
!
!                                              2 2    2
!                  EL2(X,KC,A,B) = AX RF(1,1+KC X ,1+X ) +
!
!                                              3         2 2    2
!                                 +(1/3)(B-A) X RD(1,1+KC X ,1+X )
!
!
!
!                  Legendre form of alternative ELLIPTIC INTEGRAL of 2nd
!                  -----------------------------------------------------
!                        kind
!                        ----
!
!                            Q     2       2   2  -1/2
!                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
!                            0
!
!
!
!                                   3          2     2   2
!                  D(Q,K) =(1/3)(SIN Q)  RD(COS Q,1-K SIN Q,1)
!
!
!
!
!
!                  Lemniscate constant B
!                  ---------------------
!
!
!
!                       1    2    4 -1/2
!                  B = INT  S (1-S )    DS
!                       0
!
!
!                  B =(1/3)RD (0,2,1)
!
!
!
!
!                  Heuman's LAMBDA function
!                  ------------------------
!
!
!
!                  (PI/2) LAMBDA0(A,B) =
!
!                                       2                2
!                     = SIN(B) (RF(0,COS (A),1)-(1/3) SIN (A) *
!
!                               2              2         2       2
!                      *RD(0,COS (A),1)) RF(COS (B),1-COS (A) SIN (B),1)
!
!                               2       3            2
!                     -(1/3) COS (A) SIN (B) RF(0,COS (A),1) *
!
!                             2         2       2
!                      *RD(COS (B),1-COS (A) SIN (B),1)
!
!
!
!                  Jacobi ZETA function
!                  --------------------
!
!
!                             2                2       2   2
!                  Z(B,K) = (K/3) SIN(B) RF(COS (B),1-K SIN (B),1)
!
!
!                                      2            2
!                             *RD(0,1-K ,1)/RF(0,1-K ,1)
!
!                               2       3          2       2   2
!                            -(K /3) SIN (B) RD(COS (B),1-K SIN (B),1)
!
!
!    -------------------------------------------------------------------
!
!***REFERENCES  B. !. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. !. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. !. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Modify calls to XERMSG to put in standard form.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RD
	CHARACTER*16 XERN3, XERN4, XERN5, XERN6
	INTEGER IER
	REAL(8) LOLIM, UPLIM, EPSLON, ERRTOL
	REAL(8) C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
	REAL(8) MU, POWER4, SIGMA, S1, S2, X, XN, XNDEV
	REAL(8) XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, ZNROOT
	LOGICAL FIRST
	SAVE ERRTOL, LOLIM, UPLIM, C1, C2, C3, C4, FIRST
	DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  RD
	IF (FIRST) THEN
!         ERRTOL = (R1MACH(3)/3.0E0)**(1.0E0/6.0E0)
!         LOLIM  = 2.0E0/(R1MACH(2))**(2.0E0/3.0E0)
!         TUPLIM = R1MACH(1)**(1.0E0/3.0E0)
!         TUPLIM = (0.10E0*ERRTOL)**(1.0E0/3.0E0)/TUPLIM
!         UPLIM  = TUPLIM**2.0E0
!
	C1 = 3.0E0/14.0E0
	C2 = 1.0E0/6.0E0
	C3 = 9.0E0/22.0E0
	C4 = 3.0E0/26.0E0
	ENDIF
	FIRST = .FALSE.
!
!         CALL ERROR HANDLER IF NECESSARY.
!
	RD = 0.0E0
	IF( MIN(X,Y).LT.0.0E0) THEN
	IER = 1
	WRITE (6, '(1PE15.6)') X
	WRITE (6, '(1PE15.6)') Y
!         CALL XERMSG ('SLATEC', 'RD',
!     *      'MIN(X,Y).LT.0 WHERE X = ' // XERN3 // ' AND Y = ' //
!     *      XERN4, 1, 1)
	RETURN
	ENDIF
!

	IER = 0
	XN = X
	YN = Y
	ZN = Z
	SIGMA = 0.0E0
	POWER4 = 1.0E0
	errtol = 1.e-5
!
30 	MU = (XN+YN+3.0E0*ZN)*0.20E0
	XNDEV = (MU-XN)/MU
	YNDEV = (MU-YN)/MU
	ZNDEV = (MU-ZN)/MU
	EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
	IF (EPSLON.LT.ERRTOL) GO TO 40
	XNROOT = dsqrt(XN)
	YNROOT = dsqrt(YN)
	ZNROOT = dsqrt(ZN)
	LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
	SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
	POWER4 = POWER4*0.250E0					
	XN = (XN+LAMDA)*0.250E0
	YN = (YN+LAMDA)*0.250E0
	ZN = (ZN+LAMDA)*0.250E0
	GO TO 30
!
40 	EA = XNDEV*YNDEV
	EB = ZNDEV*ZNDEV
	EC = EA - EB
	ED = EA - 6.0E0*EB
	EF = ED + EC + EC
	S1 = ED*(-C1+0.250E0*C3*ED-1.50E0*C4*ZNDEV*EF)
	S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
	RD = 3.0E0*SIGMA + POWER4*(1.0E0+S1+S2)/(MU* dsqrt(MU))
!
	RETURN
    END


    
    function gauss_random(average,sigma)
		real(8) gauss_random1,gauss_random2
		real(8) average,sigma,mm,nn,w
2		call random_number(xsl)
		mm = 2.0*xsl-1.0
		call random_number(xsl)
		nn = 2.0*xsl-1.0
		w = mm*mm + nn*nn
		if(w.gt.1.0) goto 2
		if(w.lt.3.0e-7) goto 2
		gauss_random = mm*sqrt((-2.0*log(w))/w)*sigma+average
	end
