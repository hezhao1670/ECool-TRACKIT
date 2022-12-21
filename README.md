Updated 2022/09-(ecoolm2h4 build)

# ecool-trackit 
electron cooling simulation (Multi-particles tracking code)


# Features: 
1. basic betatron motion, synchrotron motion, chromaticity, etc..
2. cooling: 3D non-magnetic cooling and dispersion
3. IBS: B-M IBS model with vertical dispersion
4. cooling section: several sections for hadron beam evolution (e-beam no change)  
5. add magnetized cooling force (Parkhomchuk, uniform e-beam) and energy/frequency tuning 
6. multi-injection (iichoice=3)
7 added various e-beam distribution uniform/elliptical/gaussian/hollw etc. 
		
		
# Publications:
Phys. Rev. Accel. Beams 23, 074201
Phys. Rev. Accel. Beams 24, 083502


# How to use:

1. In Linux machine:
Compile:
	1) move to the directory where the source code is located
	2) run the command
	: ifort -o ecool ecoolm2h4.f90 getIBS.f90 e-distribution.f90 pulsecoolh1.f90
	3) ecool is the output run file
Run:
	1) create a subdirectory
	 :mkdir s1
	2) move to this directory
	 :cd s1
	3) copy input and optics file 
	 :cp ../ecoolm2h4.in .
	 :cp ../optics.dat .
	4) edit the input file and add optics.dat
	 :gedit ecoolm2h4.in &
	5) run code
	 :../ ecoolm2h2
	 
2. In Window machine:
	1) Just prepare all the source files and compile it by any compilors such as ifort or IFC
	2) Prepare the optics file and the input file, and then run the exe file.


# File explanation
1. The input file looks like:
400000 6000 1000  1 19042727   nturns,ndim,nwrite,nperturn,iseed
2.0402  128.8  1.322  3.e3  1  -0.e6  7560   gammat,circ,gamma0,vrf,nharm,vrf2,nharm2
2.898  3.76e-3  1.679  4.94e-3  2. 100 100 tuney,ampy,tunex,ampx,chrom,nresamp,nresamp2
5.0e7  12.   6.  3.  5.0e-8   4.8e-4  1  pnumber,aatom,qatom,power,tauhat,coastbdpp,nturnon
1  1   iichoice,iechoice
1. 1.  1.  0.00  70.e-8 fracibstot, fracibsx ,fracecool,dqmin,thib
20.  1e-3  1e-3 ebeams,ebeamx, ebeamy
0 1 0.5e-3 0.5e-3 4.e-4 18.e-7 3.4 380 idosc,idoecool,sigurestx,siguresty,sigurestt,qbeamcool,slencool,speedfact
0 0  10. 10. 1  dispxe,dispxh,betaxcool,betaycool,nslice
8.0  40.5 100.0 	betamax,DA,lifetime
0.0e-3 0.0e-4  0.0e-4  xoff, xpoff, uoff
=============== pulsed i-injection (iichoice=3) ==================
0.25e-6  10      pulseiW, gapit  
0.3e-6	      spareW
=============== pulsed e-beam (iechoice=1) ==================
0.1  1.0   0.5  1e-5  tie(A), aeb(cm), temet(eV), temel(eV)
1     0.1   1.0  nedis(1-uni;2-elli;3-para;4-cossq;5-gau;6-holl), hollowsig(cm), hollowaeb(cm)
0.178   5e-5    0.    Bg(T), Berror, feshift
0.6568e-6   1.0   0.e-9   0e-9  0.   pulseW, pulseN, riseW, gapW
=============== ebeam energy tuning ======================
1.e0  -4.0e-4  4.e-4   0  fetune,eklo,ekhi,ntchoice(square-saw-tri-sin)
0    ntchoice (off-square-saw-tri-sin)	#tuning choice
0.e0  -4.0e-4  4.e-4    fetune,eklo,ekhi	#energy tuning
800.e3  -2.0e-8  2.e-8    fptune,eplo,ephi	#position tuning  
        


nturns = number of turns to simulate
ndim = number of macro-particles in ion bunch
nwrite = write beam parameters to csmon.out every nwrite turns
nperturn= tracking step (here nperturn always be 1, no other choice)
Iseed = random number seed
gammat = transition gamma for ring
circ   = circumference of ring [unit=m]
gamma0 = central gamma for hadrons
vrf    = amplitude of rf voltage [unit=V]
nharm  = harmonic number (120 for 9 MHz)
vrf2   = amplitude of rf1 voltage [unit=V]
nharm2 = second harmonic (360 for 28 MHz) 
tuney  = verical betatron tune
ampy   = initial rms beam radius y at beta function = circ/(2*pi*tuney) [unit=m]
tunex  = horizontal betatron tune
ampx   = initial rms beam radius x at beta function = circ/(2*pi*tunex) [unit=m]
chrom  = chromaticty
nresamp  = number of longitudinal bins for hadron beam 
nresamp2 = same as nresamp
pnumber  = number of ions per
aatom  = atomic mass of ion
qatom  = atomic number of ion
taupart = longituinal smoothing length for writing average profiles. this is used to smooth the longitudinal distribution [unit=s]
power  = bunch shape parameter to generate initial longitudinal distribution
tauhat = half bunch length (~3sigma) (dpp_ion will matched based on tauhat and RF setting) [unit=s]
coastbdpp = dp/p for coasting beam (iichoice = 1)
nturnon= 1 number of turns to turn on nonlinear RF and longitudinal wakes
iichoice = ion beam choice (1-coasting beam;  2-bunched beam; 3-pulsed ion for multi-injection)
	note: all choices using the same input for beam emittance, longitudinal is different!
		coasting beam: using coastbdpp
		bunched beam: using tauhat/ vrf/nharm/vrf2/nharm2	
iechoice = cooling force mode (1-magnetic cooling;  2-non-magnetic cooling)
	note:	magnetic cooling: using the parameters in pulsed e-beam
		non-magnetic cooling: using single bunch parameters: sigurestx,siguresty,sigurestt, etc
fracibstot= 1, fraction of lonigtudinal IBS, 2 means 2 times larger of IBS rate, 0 turn off IBS 
fracibsx  = 1, fraction of transverse IBS, 0 turn off IBS 
fracecool = 1/0, turn on/off magnetic cooling for heating test
dqmin  = minimum actual tune difference with betatron coupling
thib   = length of interval, typically the RF period (for longitudinal loss) [unit=s]
ebeams = ebeam rms bunch length [unit=m] (only for bunched e-beam)
ebeamx = ebeam rms size x at betax_star [unit=m]
ebeamy = ebeam rms size y at betay_stat [unit=m]
idosc = 2 if you want ion space charge included, 0 otherwise
idoecool = 1 for electron cooling, 0 for no cooling
sigurestx = rms angle spread = sqrt(emitt_x/betax_star)
siguresty = rms angle spread = sqrt(emitt_y/betay_star)
sigurestt = rms value of v_z/c in rest frame = rms(p)/p in lab
qbeamcool = bunch charge used for cooling
slencool = total length of cooling section
speedfact = speed up factor (total tracking turns = nturns*speedfact)
dispxe = electron dispersion
dispxh = hadron dispersion
betaxcool = hori. beta function for hadron ring (beta_star)
betaycool = vert.. beta function for hadron ring (beta_star)
nslice = split the cooling section into nslice sections (for hadron beam)
betamax = maximum beta function for particle loss checking sigma0 = sqrt(emitt*betamax)
DA = particle loss if position_x > DA*sigma0
plifetime = ion beam lifetime manually set
xoff = hori. position offset between e_beam and hadron beam
xpoff = hori. Velocity offset between e_beam and hadron beam
uoff = energy offset between e_beam and hadron beam
=============== pulsed i-injection (iichoice=3) ==================
pulseiW = Ion beam pulse width [unit=s]
gapit =   injection interval time [unit=s]
spareW= spare space for the new injection [unit=s]
=============== pulsed e-beam (iechoice=1) ==================
tie(A) = pulse ebeam current
aeb(cm) = pulse ebeam radius (half width)
temet(eV) = temperature_tran. 
temel(eV) = temperature_long.
nedis= ebeam trans. distribution (1-uniform;2-elliptical;3-parabolic;4-cos-square;5-gau;6-hollow) hollowsig = hollow beam sigma [unit=cm]
hollowaeb =hollow beam center position [unit=cm]
Bg(T) = magnetic field
Berror = field error 
feshift = ebeam frequency shift/f_0
pulseW = pulse width [unit=s]
pulseN = number of pulses
riseW = edge rising width [unit=s]
gapW = gap width between pulses [unit=s]
vkick = longitudinal kick voltage at the edge[unit=V]
=============== ebeam tuning ======================
ntchoice = ebeam tuning choice (0-off, 1-square, 2-saw, 3-tri, 4-sin)
fetune = energy tuning frequency of e-beam [unit=Hz]
eklo = energy tuning low limit (dE/E)
ekhi = energy tuning high limit (dE/E)
fptune= position tuning frequency of e-beam [unit=Hz]
eplo = position tuning low limit [unit=m] (only in horizontal offset tuning) 
ephi = position tuning high limit [unit=m]



2. As an input file, the lattice file “optics.dat” must have the format as below:
The first two parameters (for example: 3011 2) define the length of the lattice data and the first row of the data.
Data format:  S(m)    BETX(m)    ALFX    BETY(m)    ALFY    DX(m)    DPX    DY(m)    DPY

3. the output “elbeam.dat” is beam distribution of e-beam based on the beam dis choice
1-aeb = ebeam radius
2-elsdis(i)*ple  = density 
3-sldis(i) = transverse velocity distribution vs. r
4-ssdis(i) = longitudinal velocity distribution (dv/v) vs. r 
   

4. the output file csmon.out: give the evolusion of hadron beam parameters
1-kturns  = number of turns (In fact, total turns= kturns*speedfact )
2-csfull  = emitty
3-csfull2 = emittx
4-sigtime = bunch length [unit=s]
5-siggama = dp/p
6-peakilumi = estimate peak lumi.
7-fwhmlumi  = estimate lumi. 
8-float(np) = number of remaining macro-particles
9-float(ntloss) = number of particle lost in transverse





