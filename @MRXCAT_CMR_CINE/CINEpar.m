function CINEpar( MRX, filename )
% This function is the parameter file for MRXCAT_CMR_CINE.
% Change parameters in section "MRXCAT settings" only
%
% Note: Not all combinations of any parameter values are possible.
%       Some parameter changes require changes in the XCAT files.
%       E.g. if you want to increase the number of segments, you need
%       more XCAT heart phases for the additional segments, 
%       i.e. additional	*.bin files.


% --------------------------------------------------------------------
%   MRXCAT settings
% --------------------------------------------------------------------
RhoMuscle   = 80.0;                         % Proton density muscle [%]
RhoFat      = 70.0;                         % Proton density fat    [%]
RhoBlood    = 95.0;                         % Proton density blood  [%]
RhoLiver    = 90.0;                         % Proton density liver  [%]
RhoBone     = 12.0;                         % Proton density bone   [%]

T1Muscle    = 900.0;                        % T1 muscle             [ms]
T1Fat       = 350.0;                        % T1 fat                [ms]
T1Blood     = 1200.0;                       % T1 blood              [ms]
T1Liver     = 800.0;                        % T1 liver              [ms]
T1Bone      = 250.0;                        % T1 bone               [ms]

T2Muscle    = 50.0;                         % T2 muscle             [ms]
T2Fat       = 30.0;                         % T2 fat                [ms]
T2Blood     = 100.0;                        % T2 blood              [ms]
T2Liver     = 50.0;                         % T2 liver              [ms]
T2Bone      = 20.0;                         % T2 bone               [ms]

TR          = 3.0;                          % Repetition time       [ms]
TE          = 1.5;                          % Echo time             [ms]
Flip        = 60.0;                         % Flip                  [deg]
Frames      = 24;                           % Number of heart phases (default: 24; 0=use # XCAT frames (.bin files))
Segments    = 1;                            % Number of segments

BoundingBox = [0.2,0.6;0.3,0.7;0.0,1.0];    % BoundingBox in rel. units
RotationXYZ = [115.0;35.0;240.0];           % Rotations about x,y,z [deg]  (default: 115/35/240)
                                            % x=(RL) y=(AP) z=(FH)
LowPassFilt = 1;                            % low-pass filter images
FilterStr   = 1.2;                          % low-pass filter strength (default: 1.2)

SNR         = 20;                           % signal-to-noise ratio
Coils       = 4;                            % number coils (Biot-Savart)
CoilDist    = 450;                          % body radius           [mm]    = distance of coil centres from origin
CoilsPerRow = 8;                            % number of coils on 1 "ring" or row of coil array (default: 8)

Trajectory  = 'Cartesian';                  % k-space trajectory (Cartesian, Radial, GoldenAngle)
Undersample = 1;                            % undersampling factor (only for Radial/GoldenAngle right now)

% ---------------------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------------------------


% --------------------------------------------------------------------
%   Display title
% --------------------------------------------------------------------
fprintf ( '------------------------------------------\n' );
fprintf ( '        MRXCAT-CMR-CINE (VER %3.1f)      \n' , MRX.Version );
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
    fprintf ('Generating XCAT2 bin files...');
    fname = 'cine'; exe = 'dxcat2 ';
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

MRX.Par.tissue.t2muscle     = T2Muscle;
MRX.Par.tissue.t2fat        = T2Fat;
MRX.Par.tissue.t2blood      = T2Blood;
MRX.Par.tissue.t2liver      = T2Liver;
MRX.Par.tissue.t2bone       = T2Bone;

MRX.Par.scan.tr             = TR;
MRX.Par.scan.te             = TE;
MRX.Par.scan.flip           = pi*Flip/180;
MRX.Par.scan.segments       = Segments;
MRX.Par.scan.bbox           = BoundingBox;
MRX.Par.scan.lowpass        = LowPassFilt;
MRX.Par.scan.lowpass_str    = FilterStr;
MRX.Par.scan.snr            = SNR;
MRX.Par.scan.coils          = Coils;
MRX.Par.scan.coildist       = CoilDist;
MRX.Par.scan.coilsperrow    = CoilsPerRow;
MRX.Par.scan.rotation       = pi*RotationXYZ/180;
if Frames>0 % only overwrite Par.scan.frames, if Frames ~= 0
    MRX.Par.scan.frames     = Frames;
    xcat_segments           = round(MRX.Par.scan.scan_dur/MRX.Par.scan.heartbeat_length);
    frames_max              = MRX.Par.scan.frames_xcat/xcat_segments;
    MRX.Par.scan.phases     = round( linspace(1,frames_max,Frames) );
else
    MRX.Par.scan.frames     = MRX.Par.scan.frames_xcat/MRX.Par.scan.segments;
    MRX.Par.scan.phases     = 1:MRX.Par.scan.frames;
end
MRX.Par.scan.trajectory     = Trajectory;
MRX.Par.scan.undersample    = Undersample;

% --------------------------------------------------------------------
%   Error checks
% --------------------------------------------------------------------
if mod(max(MRX.Par.scan.frames),1)>0 % check if #frames is whole number
    error('Number of frames must be an integer value. Check number of segments in CINEpar.m and number of XCAT .bin files!')
end
if exist('frames_max','var') && frames_max < MRX.Par.scan.frames
    error('Number of XCAT phases < desired number of phases. Set Frames <= %d in CINEpar.m',frames_max)
end

end