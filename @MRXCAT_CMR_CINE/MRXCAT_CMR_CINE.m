%=========================================================================
% 
% MRXCAT_CMR_CINE MRXCAT class for cine cardiac imaging
% 
% The MRXCAT_CMR_CINE class contains specific cine MRXCAT methods.
% Methods common with other classes are in the MRXCAT superclass.
% 
% PROPERTIES:	
%
%           All properties are listed in MRXCAT.m
%       
% METHODS:  (MRX refers to the MRXCAT instance)
% 
%           CINEpar( MRX, filename )    External parameter function (CINEpar.m)
%                                       Use CINEpar.m to modify MRXCAT parameters
%           computeNoiseStdDev( MRX, sen )                 
%                                       Compute standard deviation of object noise based on desired SNR
%                                       in MRX.Par.scan.snr. Signals are averaged over the heart only 
%                                       to get "cardiac SNR".
%           mapTissueProps( MRX, data ) Assign tissue properties (T1, T2, rho) and apply signal model.
%                                       For cine, the signal model is a balanced SSFP sequence. 
%                                       Create tissue masks with different XCAT mask values.
%           extractSegment( MRX, segm, img )                 
%                                       Extract k-space segment segm from image img (in image space).
%                                       Total no. of segments: MRX.Par.scan.segments 
%           radialResample( MRX, img )  Calculate and apply radial trajectory. The Cartesian input
%                                       is resampled to a radial trajectory (MRX.Par.scan.trajectory)
%                                       including optional undersampling (factor MRX.Par.scan.undersample).
% 
% STATIC METHODS:
% 
%           getRadWeights2D( no_samples, no_profiles, alt_prof, gafl, normfl )
%                                       Calculate sampling weights for radial trajectory. For more details,
%                                       cf. function help.
%           buildRadTraj2D( no_samples, no_profiles, alt_prof, gafl, normfl, dim_t, t_offset, dim_z, z_offset)
%                                       Create radial trajectory. For more details, cf. function help.
% 
% 
% WEBSITE: 	http://www.biomed.ee.ethz.ch/mrxcat
% 
%=========================================================================

%=========================================================================
% VERSION HISTORY:						
%		130129SK INITIAL VERSION
%		130326LW DXCAT2 W/ ANGULATION - v0.1
%       130625LW ADAPTATION TO MRXCAT CMR PERF - v0.8
%       130828LW RADIAL TRAJECTORIES - v0.9
%       140130LW OO IMPLEMENTATION - v1.0
%       140202LW SEGMENTED ACQ, LOW-PASS FILTER - v1.1
% 
% AUTHORS:  Lukas Wissmann, Sebastian Kozerke, Claudio Santelli
%           Institute for Biomedical Engineering, University and ETH Zurich
%
%=========================================================================

classdef MRXCAT_CMR_CINE < MRXCAT
    
    properties
        
    end
    
    methods
        function MRX = MRXCAT_CMR_CINE( filename, varargin )
            
            MRX = MRX@MRXCAT();
            
            % ------------------------------------------------------------
            % Load and assign parameters (file 
            % ------------------------------------------------------------
            if ~exist('filename','var'), filename = ''; end
            MRX.CINEpar( filename );
            
            % ------------------------------------------------------------
            % Check for additional inputs
            % ------------------------------------------------------------
            for i = 1:length( varargin )
                if any( strcmpi( varargin{i}, {'snr'} ))
                    MRX.Par.scan.snr = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'frames'} ))
                    MRX.Par.scan.frames = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'coils'} ))
                    MRX.Par.scan.coils = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'flip'} ))
                    MRX.Par.scan.flip = pi/180*varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'demo_gui'} )) % overwrite PARs from demo GUI
                    par = varargin{i+1};
                    names1 = fieldnames(par.contrast);
                    for n=1:length(names1)
                        if isfield( par.contrast, names1{n} )
                            MRX.Par.contrast.(names1{n}) = par.contrast.(names1{n});
                        end
                    end
                    names2 = fieldnames(par.scan);
                    for n=1:length(names2)
                        if isfield( par.scan, names2{n} )
                            MRX.Par.scan.(names2{n}) = par.scan.(names2{n});
                        end
                    end
                    % currently no tissue parameters modifiable by GUI
                end
            end
            
            % --------------------------------------------------------------------
            %   Calculate coil sensitivies (Biot-Savart)
            % --------------------------------------------------------------------
            sen = MRX.calculateCoilMaps; sen = single(sen);
            
            % --------------------------------------------------------------------
            % Compute standard deviation factor for noise addition (SNR way):
            % stddev = mean ROI signal / desired SNR
            % --------------------------------------------------------------------
            MRX.Par.scan.noisestd = MRX.computeNoiseStdDev( sen );
            
            st = 0; % initialize waitbar
            
            % --------------------------------------------------------------------
            %   Produce MRXCAT phantom loops over heart phases & k-space segments
            % --------------------------------------------------------------------
            for t=MRX.Par.scan.phases % heart phases 
                for s=1:MRX.Par.scan.segments % k-space segments
                    
                    % --------------------------------------------------------------------
                    %   Read data
                    % --------------------------------------------------------------------
                    xcat_no = t+((s-1)*MRX.Par.scan.frames);
                    data = MRX.readImgData(xcat_no); data = single(data);
                    
                    % ----------------------------------------------------------------
                    %   Map MR tissue properties
                    % ----------------------------------------------------------------
                    [img,msk] = MRX.mapTissueProps(data); img = single(img); msk = single(msk);
                    
                    % ----------------------------------------------------------------
                    %   Low Pass Filter (Blur)
                    % ----------------------------------------------------------------
                    [img,msk] = MRX.lowPassFilter(img,msk);
                    
                    % ----------------------------------------------------------------
                    %   Add coils
                    % ----------------------------------------------------------------
                    img       = MRX.multiplyCoilMaps(img,sen);
                    
                    % ----------------------------------------------------------------
                    %   Add noise
                    % ----------------------------------------------------------------
                    [img,noi] = MRX.addNoise(img);
                    
                    % ----------------------------------------------------------------
                    %   Extract needed k-space segment
                    % ----------------------------------------------------------------
                    if s==1; MRX.Ksp = zeros(size(img)); end
                    MRX.extractSegment( s, img );
                    if s==MRX.Par.scan.segments %% FFT after k-space is filled
                        img = MRXCAT.k2i( MRX.Ksp, [1 2 3]);
                    end
					
                end
                MRX.Ksp = [];
                
                % -----------------------------------------------------------------
                %   Regrid data for radial trajectory
                % -----------------------------------------------------------------
                if strcmpi(MRX.Par.scan.trajectory, 'Radial') || strcmpi(MRX.Par.scan.trajectory, 'GoldenAngle')
                    img = MRX.radialResample( img );
                end
                
                % -----------------------------------------------------------------
                %   Save data to .mat file
                % -----------------------------------------------------------------
                MRX.saveImgData(img, msk, noi, sen, t);
                st = MRX.textwaitbar( t, st, 'Writing MRXCAT output data'); % update waitbar
            end
            
            % --------------------------------------------------------------------
            %   Save Parameters for Recon
            % --------------------------------------------------------------------
            MRX.saveParsForRecon( img );
            
        end
                
        
        %=========================================================================
		% Calc noise std deviation based on desired SNR
        function stdev = computeNoiseStdDev( MRX, sen )
            % STDEV = COMPUTENOISESTDDEV( MRX, SEN )
            %   Compute noise standard deviation based on desired SNR and signal
            %   intensity in ROI
            %
            %   Input:  MRX ...     MRXCAT instance
            %           sen ...     coil sensitivity maps
            % 
            %   Output: STDEV ...   Std. deviation of noise
            % 
            % Note: noise is not added per segment, because it doesn't have to be the 
            % true noise, but just simulated noise
            
            for t=MRX.Par.scan.phases
                data        = MRX.readImgData(t); data = single(data);
                [img,msk]   = MRX.mapTissueProps(data); img = single(img); msk = single(msk);
                img         = MRX.multiplyCoilMaps(img,sen);
                roi         = msk==MRX.Par.act.myoLA_act|msk==MRX.Par.act.myoLV_act| ... 
                    msk==MRX.Par.act.myoRA_act|msk==MRX.Par.act.myoRV_act;
                % find biggest connected mask = heart (if > 1 connected region)
                conn        = bwconncomp(roi);
                if conn.NumObjects > 1
                    for k=1:conn.NumObjects
                        len(k)  = length(conn.PixelIdxList{k});
                    end
                    [~,maxk]    = max(len);
                    idxs        = conn.PixelIdxList{maxk};
                    roi         = zeros(size(roi));
                    roi(idxs)   = 1;
                end
                roiall(:,:,:,:,t) = repmat(roi,[1,1,1,size(sen,4)]).*img;
            end
            
            for k=1:MRX.Par.scan.coils, roik = roiall(:,:,:,k,:); smean(k) = mean(roik(roik~=0)); end
            sosmean     = sqrt(sum(abs(smean).^2));
            stdev       = 1/MRX.Par.scan.snr*sosmean; 
            fprintf('adding noise w/ std dev : %f\n',stdev);
            
        end
        
        
        %=========================================================================
		% Calc signal intensities using tissue and sequence pars
        function [img,msk] = mapTissueProps(MRX, data)
            
            act = cell2mat(struct2cell(MRX.Par.act));
            tis = fieldnames(MRX.Par.act);
            img = single(zeros(size(data)));
            msk = uint8(zeros(size(data)));
            
            for i=1:length(act)
                
                rho = 0;
                % ----------------------------------------------------------------
                %   Select tissue type
                % ----------------------------------------------------------------
                switch char(tis(i))
                    
                    % ------------------------------------------------------------
                    %   Myocardium
                    % ------------------------------------------------------------
                    case {'myoLV_act','myoRV_act','myoLA_act','myoRA_act'}
                        rho = MRX.Par.tissue.rhomuscle;
                        t1  = MRX.Par.tissue.t1muscle;
                        t2  = MRX.Par.tissue.t2muscle;
                        % ------------------------------------------------------------
                        %   Blood
                        % ------------------------------------------------------------
                    case {'bldplLV_act','bldplRV_act','art_activity','vein_activity','bldplLA_act','bldplRA_act'}
                        rho = MRX.Par.tissue.rhoblood;
                        t1  = MRX.Par.tissue.t1blood;
                        t2  = MRX.Par.tissue.t2blood;
                        % ------------------------------------------------------------
                        %   Body
                        % ------------------------------------------------------------
                    case {'body_activity','pericardium_activity'}
                        rho = MRX.Par.tissue.rhofat;
                        t1  = MRX.Par.tissue.t1fat;
                        t2  = MRX.Par.tissue.t2fat;
                    case 'muscle_activity'
                        rho = MRX.Par.tissue.rhomuscle;
                        t1  = MRX.Par.tissue.t1muscle;
                        t2  = MRX.Par.tissue.t2muscle;
                    case 'liver_activity'
                        rho = MRX.Par.tissue.rholiver;
                        t1  = MRX.Par.tissue.t1liver;
                        t2  = MRX.Par.tissue.t2liver;
                    case {'rib_activity','cortical_bone_activity','spine_activity','bone_marrow_activity'}
                        rho = MRX.Par.tissue.rhobone;
                        t1  = MRX.Par.tissue.t1bone;
                        t2  = MRX.Par.tissue.t2bone;
                    otherwise
                        rho = 0;
                end
                
                % ----------------------------------------------------------------
                %   Signal model (Balanced SSFP)
                % ----------------------------------------------------------------
                a   = MRX.Par.scan.flip;
                te  = MRX.Par.scan.te;
                sig = rho*sin(a)/((t1/t2+1)-cos(a)*(t1/t2-1))*exp(-te/t2);
                
                % ----------------------------------------------------------------
                %   Update tissue compartment
                % ----------------------------------------------------------------
                img(find(data(:)==act(i))) = sig;
                
                % ----------------------------------------------------------------
                %   Update tissue masks
                % ----------------------------------------------------------------
                msk(find(data(:)==act(i))) = act(i);
                
            end
        end
        
        
        %=========================================================================
		% Extract k-space segment from fully sampled k-space
        function extractSegment(MRX, segm, img)
            nsegm    = MRX.Par.scan.segments;                   % no. of segments
            kyrange  = ceil((segm-1)*size(img,2)/nsegm)+1:min(ceil(segm*size(img,2)/nsegm),size(img,2));
            temp     = MRXCAT.i2k(img,[1 2 3]);                 % FFT
            MRX.Ksp(:,kyrange,:,:,:) = temp(:,kyrange,:,:,:);   % extract segment
        end
        
        
        %=========================================================================
        % Radial Trajectory Resampling + Calc
        function [img_rad,ksp] = radialResample( MRX, img )
            
            % standard radial or golden angle trajectory
            if strcmpi(MRX.Par.scan.trajectory, 'goldenAngle')
                ga = 1;     % golden angle
            else
                ga = 0;     % stamdard radial
            end
            
            % calculate radial trajectory
            samp = size(img,1);
            % Number of profiles acc. to radial Nyquist & undersampling factor.
            % Cf. MRXCAT.computeBoundingBox for larger FOV in radial than Cartesian.
            prof = round(size(img,1)*2/pi*1/MRX.Par.scan.undersample);  % # prof = 1/1.57*FOVx
            w    = MRXCAT_CMR_CINE.getRadWeights2D(samp,prof,0,ga,1);
            k    = MRXCAT_CMR_CINE.buildRadTraj2D(samp,prof,0,ga,1);
            % apply trajectory to image
            % Check if NUFFT toolbox (J.Fessler) and wrapper (M.Lustig) are present
            errormsg = ['NUFFT toolbox by J. Fessler and/or NUFFT wrapper ' ...
                'by M. Lustig not found in Matlab path. ' ...
                'Please download and install from ' ...
                ' http://web.eecs.umich.edu/~fessler/code/index.html' ...
                ' and http://www.eecs.berkeley.edu/~mlustig/Software.html'];
            if exist('NUFFT')<2 || exist('nufft')<2
                error(errormsg);
            else % in case of Matlab version < 2011b, exist is not case-sensitive
                [~,f0]=fileparts(which('nufft'));
                [~,f1]=fileparts(which('NUFFT'));
                % case-sensitive check
                if ~strcmp(f0,'nufft') || ~strcmp(f1,'NUFFT')
                    error(errormsg);
                end
            end

%             e    = NUFFT(k,w,[0,0],[samp,samp]);        % Claudio Santelli version (no phase argin)
            e    = NUFFT(k,w,1,[0,0],[samp,samp],1);    % Miki Lustig wrapper (phase=1)
            
            % radially resample for all coil elements
            for k=1:MRX.Par.scan.coils
                ksp(:,:,:,k)     = e*double(img(:,:,:,k));
                img_rad(:,:,:,k) = e'*ksp(:,:,:,k);
            end
        end
        
    end
    
    
    
    methods ( Static )
        
        %=========================================================================
        % Calculate sampling weights for radial trajectory
        function w = getRadWeights2D(no_samples, no_profiles, alt_prof, gafl, normfl)
            %==========================================================================
            % Function which calculates analytically the sampling weights of a radial
            % trajectory.
            %
            % Inputs:
            % -------
            % no_samples:       Number of samples along each projection.
            % no_profiles:      Number of projections per frame and slice.
            % alt_prof:         Alternating profiles flag.
            % gafl:             Golden angle flag.
            % normfl:           Normalization flag.
            %
            % Outputs:
            % --------
            % w:                Radial Weights
            %
            % Function calls:    none
            % ---------------
            %
            % Claudio Santelli, 06/11/2012
            %==========================================================================
            
            % Check input
            %--------------------------------------------------------------------------
            error(nargchk(3,5, nargin))
            if nargin<5, normfl = false; end
            if nargin<4, gafl   = false; end
            
            % Initialize weight and build up kr vector (radial coordinate along spoke)
            %--------------------------------------------------------------------------
            w = []; kr = abs(-floor(no_samples/2):ceil(no_samples/2-1));
            
            if gafl
                % Calculate angle for every spoke and map it onto interval [0,pi]
                %----------------------------------------------------------------------
                phi      = (2*pi/(sqrt(5)+1))*(0:no_profiles-1)+pi/2;
                phi      = ((phi./pi)-floor(phi./pi))*pi;
                
                % Sort spokes according to their angles on the interval [0,pi],
                % calculate relative angular distances to neighboring spokes, and
                % finally, get corresponding total relative angles.
                %----------------------------------------------------------------------
                [phi, I] = sort(phi);
                dPhi1    = [phi(2:end) (phi(1)+pi)]-phi; % Left relative angular distance
                dPhi2    = circshift(dPhi1,[0 1]); % Right relative angular distance
                dPhi     = 0.5*(dPhi1+dPhi2);
                
                % Build up weighting matrix, d, where the i-th column corresponds to
                % the i-th spoke in the sorted order.
                %----------------------------------------------------------------------
                for i=1:no_profiles, w = [w, (dPhi(i)*kr).']; end
                
                % Modify zero point and bring columns of d back into the original order,
                % i.e. i-th column corresponds to the weighting of the i-th spoke.
                %----------------------------------------------------------------------
                w(kr==0,:) = pi/(4*no_profiles);
                w(:,I)     = w;
            else
                % Calculate weight for one spoke and modify zero-point (every spoke has
                % the same weighting).
                %----------------------------------------------------------------------
                w        = (pi/no_profiles)*kr;
                w(kr==0) = pi/(4*no_profiles);
                w        = repmat(w.', [1 no_profiles]);
                if alt_prof
                    w(:,2:2:end) = w(end:-1:1,2:2:end);
                end
            end
            
            if normfl, w = w./max(w(:)); end
            
        end
        
        
        %=========================================================================
        % Create radial trajectory for NUFFT
        function k = buildRadTraj2D(no_samples, no_profiles, alt_prof, gafl, normfl, dim_t, t_offset, dim_z, z_offset)
            %==========================================================================
            % Returns radial trajectory in a complex valued data array. The real part
            % corresponds to the x- and the imaginary part to y-component respectively.
            %
            % Inputs:
            % -------
            % no_samples:   Number of samples along each projection.
            % no_rofiles:   Number of projections per frame and slice.
            % alt_prof:     Alternating profiles flag.
            % gafl:         Golden angle flag.
            % normfl:       Coordinates normalization flag.
            % dim_t:        Number of time frames.
            % t_offset:     Profile offset between two adjacent time frames.
            % dim_z:        Number of slices in z-direction.
            % z_offset:     Profile offset between two adjacent slices.
            %
            % Outputs:
            % --------
            % k:            k-space coordinates.
            %
            % Function calls:    none
            % ---------------
            %
            % Claudio Santelli, 11/06/2012
            %==========================================================================
            
            % Check input
            %--------------------------------------------------------------------------
            error(nargchk(3,9, nargin))
            if nargin<9 || isempty(z_offset), z_offset = 0; end
            if nargin<8 || isempty(dim_z),    dim_z    = 1; end
            if nargin<7 || isempty(t_offset), t_offset = 0; end
            if nargin<6 || isempty(dim_t),    dim_t    = 1; end
            if nargin<5 || isempty(normfl),   normfl   = true; end
            if nargin<4 || isempty(gafl),     gafl     = false; end
            
            % Build trajectory
            %--------------------------------------------------------------------------
            k = [];
            
            % Initial spoke along ky-axis
            k0 = [zeros(1,no_samples); -floor(no_samples/2):ceil(no_samples/2-1)];
            
            % Angle increment
            if gafl
                goldenRatio = (sqrt(5)+1)/2;
                dPhi        = pi/goldenRatio;
            else
                dPhi = pi/no_profiles;
            end
            
            for z=1:dim_z
                for t=1:dim_t
                    for i=1:no_profiles
                        % Update rotation matrix
                        rot_angle = ((i-1)+(t-1)*t_offset+(z-1)*z_offset)*dPhi;
                        if alt_prof && ~mod(i,2)
                            rot_angle = rot_angle+pi;
                        end
                        R = [cos(rot_angle), -sin(rot_angle);
                            sin(rot_angle),  cos(rot_angle)];
                        % Rotate k0 vector accordingly
                        ktmp       = (R*k0).';
                        k(:,i,z,t) = ktmp(:,1)+1i*ktmp(:,2);
                    end
                end
            end
            
            if normfl, k = k./no_samples; end
            
        end
        
    end
    
end
