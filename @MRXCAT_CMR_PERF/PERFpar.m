function PERFpar( MRX, filename )
% This function is the parameter file for MRXCAT_CMR_PERF.
% Change parameters in section "MRXCAT settings" only
%
% Note: Not all combinations of any parameter values are possible.
%       Some parameter changes require changes in the XCAT files.
%       E.g. if you want to change orientation of your phantom, you
%       need to change RotationXYZ parameter, but you also need to
%       create new XCAT masks (*.bin files)


% --------------------------------------------------------------------
%   MRXCAT settings
% --------------------------------------------------------------------
RhoMuscle   = 80.0;                         % Proton density muscle [%]
RhoFat      = 70.0;                         % Proton density fat    [%]
RhoBlood    = 95.0;                         % Proton density blood  [%]
RhoLiver    = 90.0;                         % Proton density liver  [%]
RhoBone     = 12.0;                         % Proton density bone   [%]

T1Muscle    = 800.0;                        % T1 muscle             [ms]
T1Fat       = 350.0;                        % T1 fat                [ms]
T1Blood     = 1200.0;                       % T1 blood              [ms]
T1Liver     = 800.0;                        % T1 liver              [ms]
T1Bone      = 250.0;                        % T1 bone               [ms]

Relaxivity  = 5.6 / 1000;                   % Gd relaxivity         [l/(mmol*ms)]
MBFrest     = 1.0 / 60;                     % Rest MBF              [ml/g/s]
MBFstress   = 3.5 / 60;                     % Stress MBF            [ml/g/s]
Fermi_alpha = 0.25;                         % Fermi model parameter alpha
Fermi_beta  = 0.25;                         % Fermi model parameter beta
Tshift      = 3;                            % temporal LV-myo shift [s]

TR          = 2.0;                          % Repetition time       [ms]
Trrc        = 1.0;                          % Duration r-r cycle    [s]
Flip        = 15.0;                         % Flip angle            [deg]
Tsat        = 150.0;                        % Saturation delay      [ms]
Nky0        = max(floor((Tsat-70)/TR),1);   % Number excitations to ky0 (70 ms: eff. prepulse delay)
Frames      = 32;                           % Number of Dynamics (default: 32; 0=use # XCAT frames (.bin files))

BoundingBox = [0.0,0.7;0.2,0.8;0.0,1.0];    % BoundingBox in rel. units
RotationXYZ = [115.0;35.0;240.0];           % Rotations about x,y,z [deg]  (default: 115/35/240) => z=225 for comp-based
                                            % x=(RL) y=(AP) z=(FH)
CropProfs   = 0;                            % Crop profiles along x around heart (Recon time)
LowPassFilt = 0;                            % low-pass filter images
FilterStr   = 1.2;                          % low-pass filter strength (default: 1.2)

RespMotion  = 0;                            % 0=no motion;1=resp motion
RestStress  = 2;                            % 1=rest; 2=stress
Dose        = 0.075;                        % [mmol/kg b.w.]
SNR         = 30;                           % signal-to-noise ratio (CNR in PERF case!)
Coils       = 4;                            % number coils (Biot-Savart)
CoilDist    = 350;                          % body radius           [mm]    = distance of coil centres from origin
CoilsPerRow = 8;                            % number of coils on 1 "ring" or row of coil array (default: 8)

Trajectory  = 'Cartesian';                  % k-space trajectory (Cartesian), only Cartesian for perfusion so far
Undersample = 1;                            % undersampling factor (not for Cartesian right now)

% ---------------------------------------------------------------------------------------------------------------------
% Define AIF at half dose (pop avg of 6 healthy volunteers, 0.05 mmol/kg b.w. dose, first-pass by gamma-variate fit)
% ---------------------------------------------------------------------------------------------------------------------
aif005 = [0,0,0,0,0.00088825,0.017793,0.13429,0.5466,1.453,2.8276,4.3365,5.5112,6.0165,5.7934,5.0205,3.9772,2.9161, ...
    1.9988,1.2913,0.79165,0.46315,0.25983,0.14035,0.073254,0.037056,0.018216,0.0087227,0.0040768,0.0018632, ...
    0.00083405,0.00036621,0.00015793,6.6972e-05,2.7956e-05,1.1499e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

% ---------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------



% --------------------------------------------------------------------
%   Display title
% --------------------------------------------------------------------
fprintf ( '------------------------------------------\n' );
fprintf ( '      MRXCAT-CMR-PERFUSION (VER %3.1f)      \n' , MRX.Version );
fprintf ( '------------------------------------------\n' );

% --------------------------------------------------------------------
%   Open window, select file
% --------------------------------------------------------------------
if ~exist('filename','var') || ~exist(filename,'file')
    [filename,pathname] = uigetfile({'*.bin;*.par','XCAT2 files (*.par,*.bin)'});
else
    pathname = pwd;
    if ispc, pathname = [pathname,'\']; else pathname = [pathname,'/']; end
end

% --------------------------------------------------------------------
%   Generate XCAT2 *.bin files
% --------------------------------------------------------------------
if strcmp(filename(end-2:end),'par')
    %         fname = filename(1:end-4); exe = './dxcat2 ';
    fname = 'perfusion'; exe = './dxcat2 ';
    fprintf ('Generating XCAT2 bin file --> %s ... ',fname);
    if ispc  exe = 'dxcat2.exe '; end
    if ismac exe = 'dxcat2mac ';  end
    x = num2str(RotationXYZ(1));
    y = num2str(RotationXYZ(2));
    z = num2str(RotationXYZ(3));
    s = [exe,filename,' --phan_rotx ',x,' --phan_roty ',y,' --phan_rotz ',z,' ',fname];
    cd(pathname); system(s); clear x y z s;
    fprintf ('ok\n');
    filename = [fname,'_act_1.bin'];
    clear x y z s fname;
end

% --------------------------------------------------------------------
%   Read log file
% --------------------------------------------------------------------
MRX.Filename = [pathname filename];
MRX.readLogFile;

% --------------------------------------------------------------------
%   Store tissue, contrast and sequence parameters
% --------------------------------------------------------------------
MRX.Par.tissue.rhomuscle    = RhoMuscle;
MRX.Par.tissue.rhofat       = RhoFat;
MRX.Par.tissue.rhoblood     = RhoBlood;
MRX.Par.tissue.rholiver     = RhoLiver;
MRX.Par.tissue.rhobone      = RhoBone;
MRX.Par.tissue.t1muscle     = T1Muscle;
MRX.Par.tissue.t1fat        = T1Fat;
MRX.Par.tissue.t1blood      = T1Blood;
MRX.Par.tissue.t1liver      = T1Liver;
MRX.Par.tissue.t1bone       = T1Bone;

MRX.Par.contrast.ry         = Relaxivity;
MRX.Par.contrast.qr         = MBFrest;
MRX.Par.contrast.qs         = MBFstress;
MRX.Par.contrast.rs         = RestStress;
MRX.Par.contrast.falpha     = Fermi_alpha;
MRX.Par.contrast.fbeta      = Fermi_beta;
MRX.Par.contrast.tshift     = Tshift;
MRX.Par.contrast.dose       = Dose;
MRX.Par.contrast.aif        = aif005;

MRX.Par.scan.trep           = TR;
MRX.Par.scan.trrc           = Trrc;
MRX.Par.scan.flip           = pi*Flip/180;
MRX.Par.scan.tsat           = Tsat;
MRX.Par.scan.nky0           = Nky0;
MRX.Par.scan.bbox           = BoundingBox;
MRX.Par.scan.crop           = CropProfs;
MRX.Par.scan.lowpass        = LowPassFilt;
MRX.Par.scan.lowpass_str    = FilterStr;
MRX.Par.scan.resp           = RespMotion;
MRX.Par.scan.snr            = SNR;
MRX.Par.scan.coils          = Coils;
MRX.Par.scan.coildist       = CoilDist;
MRX.Par.scan.coilsperrow    = CoilsPerRow;
MRX.Par.scan.trajectory     = Trajectory;
MRX.Par.scan.undersample    = Undersample;
MRX.Par.scan.rotation       = pi*RotationXYZ/180;
if Frames>0 %&& ~RespMotion % only overwrite Par.scan.frames, if Frames ~= 0
    MRX.Par.scan.frames     = Frames;
else
    MRX.Par.scan.frames     = MRX.Par.scan.frames_xcat;
end

% --------------------------------------------------------------------
%   Error checks
% --------------------------------------------------------------------
if mod(MRX.Par.scan.frames,1)>0 % check if #frames is whole number
    error('Number of frames must be an integer value!')
end

end
