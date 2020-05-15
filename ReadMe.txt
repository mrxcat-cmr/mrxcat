==========================================================================
    MRXCAT - A Matlab software for numerical simulation of cardiac MRI
==========================================================================

MRXCAT is a Matlab toolbox for realistic numerical simulation of cardiac 
MRI. MRXCAT employs the XCAT phantom to obtain realistic anatomical masks
including cardiac contraction and respiratory motion options. MRXCAT 
simulates dynamic tissue contrast, MR signal models, multiple receiver 
coils, object noise, Cartesian and other trajectories based on XCAT masks.
MRXCAT is currently available for cardiac cine and myocardial perfusion 
MRI simulation.


MRXCAT Showcase / Demo: HowTo
=============================
1. 	Download the MRXCAT .zip file from www.biomed.ee.ethz.ch/mrxcat,
	unpack it and add the folder to your Matlab path (only the MRXCAT 
	folder, not the @MRXCAT_CMR_PERF and @MRXCAT_CMR_CINE folders).
2.	Download the XCAT perfusion and/or cine example .zip files and unpack
	to any working directory. 
3. 	Download and install the GUI Layout Toolbox for your Matlab version
	For Matlab 2014a and earlier: 
	www.mathworks.com/matlabcentral/fileexchange/27758-gui-layout-toolbox
	For Matlab 2014b and later:
	www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
3.	Type MRD = MRXCAT_Showcase; in Matlab to start the Showcase GUI, 
	adapt parameters in the graphical user interface, click "Create 
	Dataset" and wait for the result.
	
	
MRXCAT normal mode: HowTo
=========================	
1. 	Download the MRXCAT .zip file from www.biomed.ee.ethz.ch/mrxcat,
	unpack it and add the folder to your Matlab path (only the MRXCAT 
	folder, not the @MRXCAT_CMR_PERF and @MRXCAT_CMR_CINE folders).
2.	Download the XCAT perfusion and/or cine example .zip files and unpack
	to any working directory. 
3. 	Adapt the MRXCAT parameters in CINEpar.m in @MRXCAT_CMR_CINE or 
	PERFpar.m in @MRXCAT_CMR_PERF to your needs. For a first try, go
	with the predefined parameters.
4.	Start cine or perfusion MRXCAT by typing
	MRXCAT_CMR_CINE; or MRXCAT_CMR_PERF; into the command line.
5. 	Select the first XCAT .bin file in the file selection dialog
	(cine_act_1.bin for cine, perfusion_act_1.bin for perfusion). 
	Once the simulation is done, you get the following files:
	*.cpx		MRXCAT phantom data
	*.msk		XCAT mask data
	*.sen		MRXCAT coil sensitivity maps
	*.noi		MRXCAT noise only
	*_par.mat	MRXCAT parameters
6.	To display the produced phantom, run DisplayMRXCAT; and select
	the *.cpx file in the file selection dialog.


Dependencies
============
This code was tested in Matlab R2011a and R2012a.

Some methods require the IRT toolbox by Jeffrey Fessler and the 
NUFFT wrapper by Miki Lustig for non-uniform FFT. 
If needed, download and extract them into your Matlab path from:
http://web.eecs.umich.edu/~fessler/
http://www.eecs.berkeley.edu/~mlustig/Software.html


Distribution of MRXCAT
======================
MRXCAT is distributed in Matlab source code format to support its 
propagation for research purposes. Anyone is allowed to use the MRXCAT
original and/or modified software for non-commercial purposes. If you 
publish research using MRXCAT, please cite the MRXCAT publication. 
Please do not share code derived from MRXCAT without permission of 
MRXCAT authors.

MRXCAT is available from www.biomed.ee.ethz.ch/mrxcat

(c) Copyright University and ETH Zurich; 2014 Institute for Biomedical Engineering

	
Contact details
===============
For anything related to MRXCAT:
Lukas Wissmann
mrxcat.cmr@gmail.com

For project updates: 
Sebastian Kozerke 
Institute for Biomedical Engineering 
University and ETH Zurich 
kozerke@biomed.ee.ethz.ch


Note about XCAT
===============
XCAT is a separate software, which offers a lot more than the
example *.bin datasets included in MRXCAT. The XCAT software is 
available here: https://olv.duke.edu/xcat


History
=======
yymmdd  au  version	description
----------------------------------------------------------
130127  SK  v0.1    CINE/PERF: INITIAL VERSIONS
130128  SK          CINE/PERF: COIL MAPS ADDED
130207  SK          PERF: DXCAT2 W/ ANGULATION
130219  LW          PERF: SAVE PAR STRUCT FOR RECON 
130305  LW          PERF: POPULATION AVG AIF; SIGNAL MODEL UPDATE; TSHIFT-AIF ADDED
130326  LW  v0.1    CINE: DXCAT2 W/ ANGULATION
130315  LW  v0.2    PERF: NOISE ADDITION (CNR)
130327  LW  v0.4    PERF: BIOT-SAVART FOR COIL MAPS IN 3D, NOISE UPDATE
130416  LW  v0.5    PERF: COIL MAPS DEBUG
130503  LW  v0.6    PERF: PROFILES BUGFIX, REALISTIC CONC AIF
130623  LW  v0.7    PERF: EXCHANGE ORDER OF NOISE AND COIL SENSITIVITIES ADD
130624  LW  v0.8    PERF: LOWPASS FILTER OPTION
130625  LW  v0.8    CINE: ADAPTATION TO MRXCAT CMR PERF
130828  LW  v0.9    CINE: RADIAL TRAJECTORIES
131121  LW  v0.9    PERF: FREE BREATHING W/ DIFFERENT I/O #FRAMES
140130  LW  v1.0    CINE/PERF: OO IMPLEMENTATION
140202  LW  v1.1    CINE: SEGMENTED ACQ, CINE/PERF: LOW-PASS FILTER
140421  LW  v1.2    CINE/PERF: VERSION W/ COMMENTS
150930  LW  v1.3    CINE/PERF: added MRXCAT_Showcase GUI
170111  LW  v1.4    BUGFIX coil calculus (credit: Javier Royuela del Val)