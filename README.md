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

