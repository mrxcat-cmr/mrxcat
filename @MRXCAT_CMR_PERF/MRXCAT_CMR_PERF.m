%=========================================================================
% 
% MRXCAT_CMR_CINE MRXCAT class for myocardial perfusion imaging
% 
% The MRXCAT_CMR_PERF class contains specific myocardial perfusion MRXCAT methods.
% Methods common with other classes are in the MRXCAT superclass.
% 
% PROPERTIES:	
%
%           All properties are listed in MRXCAT.m
%       
% METHODS:  (MRX refers to the MRXCAT instance)
% 
%           PERFpar( MRX, filename )    External parameter function (PERFpar.m)
%                                       Use PERFpar.m file to modify MRXCAT parameters
%           computeDynamicConc( MRX )   Compute dynamic contrast agent concentration in input (LV, AIF)
%                                       and tissue of interest (myocardium). 
%           computeNoiseStdDev( MRX, sen )                 
%                                       Compute standard deviation of object noise based on desired CNR
%                                       in MRX.Par.scan.snr. Contrast is determined by the maximum
%                                       difference in myocardial signal along time during acquisition.
%                                       To obtain the same std deviation for same CNR, a reference 
%                                       contrast agent dose of 0.075 mmol/kg b.w. is used here.
%           mapTissueProps( MRX, data ) Assign tissue properties (T1, rho) and apply signal model.
%                                       For myocardial perfusion, the signal model is a spoiled 
%                                       saturation recovery gradient echo sequence. 
%                                       Create tissue masks with different XCAT mask values.
%           updateContrastConc( MRX, t, ca, cm )                
%                                       Get contrast concentration at time t for RA, RV, LA, LV from 
%                                       LV concentration ca(t) and myocardial concentration cm(t),
%                                       update MRX. This method is used for phantom creation loop over t.
%           fermiFunction( MRX, t )     Computes the Fermi function 
%                                           f = (1+b) / ( 1+b*exp(a*t) ), 
%                                       with
%                                       a = MRX.Par.contrast.falpha
%                                       b = MRX.Par.contrast.fbeta
%                                       t = input t (usually a vector)
%           shiftAIF( MRX, t, aif )     Shift aif by dT along t.
%                                       dT is specified in MRX.Par.contrast.tshift (PERFpar).
%                                       Used to generate the time-shifted myocardial curves.
% 
% STATIC METHODS:
% 
%           convolve( t, c, h, q )      Convolves c with h along time t and scales by q. Usually,
%                                       c ... input function (AIF)
%                                       h ... impulse residue function (Fermi function)
%                                       t ... time
%                                       q ... blood flow rate
% 
% WEBSITE: 	http://www.biomed.ee.ethz.ch/mrxcat
% 
%=========================================================================

%=========================================================================
%	VERSION HISTORY:						
%		130127SK INITIAL VERSION - v0.1
%		130128SK COIL MAPS ADDED
%		130207SK DXCAT2 W/ ANGULATION
%       130208SK SINGLE PRECISION TO SAVE MEMORY
%		130219LW SAVE PAR STRUCT FOR RECON 
%       130305LW POPULATION AVG AIF; SIGNAL MODEL UPDATE; TSHIFT-AIF ADDED
%       130315LW NOISE ADDITION (CNR) - v0.2
%       130327LW BIOT-SAVART FOR COIL MAPS IN 3D, NOISE UPDATE - v0.3/v0.4
%       130416LW COIL MAPS DEBUG - v0.5
%       130503LW PROFILES BUGFIX, REALISTIC CONC AIF - v0.6
%       130623LW EXCHANGE ORDER OF NOISE AND COIL SENSITIVITIES ADD - v0.7
%       130624LW LOWPASS FILTER OPTION - v0.8
%       131121LW FREE BREATHING PERF W/ DIFFERENT I/O #FRAMES - v0.9
%       140130LW OO IMPLEMENTATION - v1.0
% 
% AUTHORS:  Lukas Wissmann, Sebastian Kozerke, Claudio Santelli
%           Institute for Biomedical Engineering, University and ETH Zurich
%
%=========================================================================

classdef MRXCAT_CMR_PERF < MRXCAT
    
    properties
        
    end
    
    methods
        function MRX = MRXCAT_CMR_PERF( filename, varargin )
            
            MRX = MRX@MRXCAT();
            
            % ------------------------------------------------------------
            % Load and assign parameters (file 
            % ------------------------------------------------------------
            if ~exist('filename','var'), filename = ''; end
            MRX.PERFpar( filename );
            
            % ------------------------------------------------------------
            % Check for additional inputs
            % ------------------------------------------------------------
            for i = 1:length( varargin )
                if any( strcmpi( varargin{i}, {'snr'} ))
                    MRX.Par.scan.snr = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'dose'} ))
                    MRX.Par.contrast.dose = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'frames'} ))
                    MRX.Par.scan.frames = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'coils'} ))
                    MRX.Par.scan.coils = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'tshift'} ))
                    MRX.Par.contrast.tshift = varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'flip'} ))
                    MRX.Par.scan.flip = pi/180*varargin{i+1};
                end
                if any( strcmpi( varargin{i}, {'crop'} ))
                    MRX.Par.scan.crop = varargin{i+1};
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
            %   Prepare arterial input function and tissue residue function
            % --------------------------------------------------------------------
            [ca, cm] = MRX.computeDynamicConc;
            
            
            % --------------------------------------------------------------------
            %   Calculate coil sensitivies (Biot-Savart)
            % --------------------------------------------------------------------
            sen = MRX.calculateCoilMaps; sen = single(sen);
            
            % --------------------------------------------------------------------
            % Compute standard deviation factor for noise addition (CNR way):
            % stddev = mean ROI signal difference between max&min myo enhancement / desired SNR
            % --------------------------------------------------------------------
            MRX.Par.scan.noisestd = MRX.computeNoiseStdDev( sen );
            
            st = 0; % initialize waitbar
            
            for t=1:MRX.Par.scan.frames % time frames (dynamics)
                
                % ----------------------------------------------------------------
                %   Read data
                % ----------------------------------------------------------------
                data = MRX.readImgData(t); data = single(data);
                
                % ----------------------------------------------------------------
                %   Update contrast concentration
                % ----------------------------------------------------------------
                MRX.updateContrastConc(t, ca, cm);
                
                % ----------------------------------------------------------------
                %   Map MR tissue properties
                % ----------------------------------------------------------------
                [img,msk]   = MRX.mapTissueProps(data); img = single(img); msk = single(msk);
                
                % ----------------------------------------------------------------
                %   Low Pass Filter (Blur)
                % ----------------------------------------------------------------
                [img,msk]   = MRX.lowPassFilter(img,msk);
                
                % ----------------------------------------------------------------
                %   Add coils
                % ----------------------------------------------------------------
                img         = MRX.multiplyCoilMaps(img, sen);
                
                % ----------------------------------------------------------------
                %   Add noise
                % ----------------------------------------------------------------
                [img,nois]  = MRX.addNoise(img);
                
                % -----------------------------------------------------------------
                %   Save data to .mat file
                % -----------------------------------------------------------------
                % Crop data in readout direction (x) around heart
                if MRX.Par.scan.crop && ~MRX.Par.scan.resp % crop
                    % find RV/LV myocardial+ventricular x indices
                    [xi,~,~] = ind2sub(size(data),find(data>=1&data<=8));
                    ranx     = [min(xi)-5:max(xi)+5];
                    img      = img(ranx,:,:,:);
                    msk      = msk(ranx,:,:);
                    nois     = nois(ranx,:,:,:);
                else % no crop
                    ranx     = 1:size(data,1);
                end
                
                MRX.saveImgData(img, msk, nois, sen(ranx,:,:,:), t);
                st = MRX.textwaitbar( t, st, 'Writing MRXCAT output data'); % update waitbar
            end
            if MRX.Par.scan.crop && ~MRX.Par.scan.resp % crop message
                fprintf('Cropped data in x to %d:%d.\n------------------------------------------\n', ... 
                    ranx(1),ranx(end));
            end
            MRX.Par.scan.crop_xprofs = ranx;
            
            % --------------------------------------------------------------------
            %   Save Parameters for Recon
            % --------------------------------------------------------------------
            MRX.saveParsForRecon( img );
            
        end
        
        
        %=========================================================================
		% Calc dynamic concentration based on conc AIF
        function [ca, cm] = computeDynamicConc( MRX )
            % Compute Gd concentration at input and ROI as function of time
            %
            % INPUT:        MRX     MRXCAT object
            %
            % OUTPUT:       ca      arterial input at dose specified in Par.contrast.dose [mmol/l]
            %               cm      myocardial tissue concentration [mmol/l]
            %
            % Note:
            % Tissue density is not used here, because the conversion from dose [mmol/kg] 
            % to c_Gd [mmol/ml] is done implicitely, i.e. we inject dose*b.w./c_Gd [ml] 
            % of contrast medium, i.e. b.w. drops out in the calc.
            

            % crop population average AIF 
            aif005 = MRX.Par.contrast.aif(1:MRX.Par.scan.frames);
            
            % --------------------------------------------------------------------
            %   Convert to absolute contrast agent concentration, using assumptions:
            %   - pure Gadovist: c_Gd = 1.0 mmol/ml = 1000 mmol/l
            %   - c = c_Gd*(dilution factor d)
            %   - estimation of dilution factor d:
            %       - 80-100 ml stroke volume
            %       - 80% ejection fraction (lost 20% not considered twice)
            %           - 64-80 ml eff. ejection volume
            %       - looking at normalized AIF shape (aif/trapz(aif)) 
            %         => peak value d ~ 0.12 = 12% (assume heart rate = 60 bpm)
            %   - absolute dose calc
            %       - e.g. patient weight 75 kg = 75*dose ml Gadovist injected 
            %         (e.g. d=0.1 => 7.5 ml Gd)
            %       - e.g. 75 ml eff. ejection volume
            %       - max. 12% of 75*dose ml Gd per heart beat 
            %         => e.g. 0.9 ml Gd / 75 ml blood => d = 12/1000
            %   - max c = c_Gd*d = 1000*12/1000 = 12 mmol/l
            % --------------------------------------------------------------------
            aif01 = aif005*12/max(aif005);
            ca = aif01*MRX.Par.contrast.dose/0.1;   % scale to desired dose
            
            % --------------------------------------------------------------------
            %   Upsample AIF to "pseudo-continuous" AIF
            % --------------------------------------------------------------------
            t       = linspace(0,(MRX.Par.scan.frames-1)*MRX.Par.scan.trrc,length(ca));
            tinf    = 0:1/100:t(end);               % upsample by factor 100
            cainf   = interp1(t,ca,tinf,'pchip','extrap');
            
            % --------------------------------------------------------------------
            %   Scale flow and calculate impulse residue function
            % --------------------------------------------------------------------
            qfl     = MRX.Par.contrast.qr;          % rest flow
            if (MRX.Par.contrast.rs>1), qfl = MRX.Par.contrast.qs; end % stress flow
            qfl     = qfl*length(ca)/round(MRX.Par.scan.frames*MRX.Par.scan.trrc);
            irf     = MRX.fermiFunction(tinf);      % impulse residue function
            
            % --------------------------------------------------------------------
            %   Calculate myocardial concentration, apply Tshift and downsample
            % --------------------------------------------------------------------
            cm = MRXCAT_CMR_PERF.convolve(tinf,cainf,irf,qfl); % myocardial conc
            cm = MRX.shiftAIF(tinf, cm); % apply Tshift
            ca = interp1(tinf,cainf,t);  % downsample AIF (discrete measurement)
            cm = interp1(tinf,cm,t);     % downsample MYO
            
        end
        
        
        %=========================================================================
		% Calc noise std deviation
        function stdev = computeNoiseStdDev( MRX, sen )
            % Calculate noise standard deviation based on desired CNR
            
            % get conc AIF/MYO for reference dose
            dosebkp     = MRX.Par.contrast.dose;
            MRX.Par.contrast.dose = 0.075; % reference dose for SNR
            [ca, cm]    = MRX.computeDynamicConc;
            MRX.Par.contrast.dose = dosebkp;
            
            % Maximum myocardial enhancement mean signal
            [~,tmax]    = max(cm);
            data        = MRX.readImgData(tmax); data = single(data);
            MRX.updateContrastConc(tmax,ca,cm);
            [img,msk]   = MRX.mapTissueProps(data); img = single(img); msk = single(msk);
            img         = MRX.multiplyCoilMaps(img,sen);
            roi         = msk==MRX.Par.act.myoLA_act|msk==MRX.Par.act.myoLV_act;
            roi         = repmat(roi,[1,1,1,size(sen,4)]).*img;
            for k=1:MRX.Par.scan.coils, roik = roi(:,:,:,k); smax(k) = mean(roik(roik~=0)); end
            sumsqmax    = sqrt(sum(abs(smax).^2));
            
            % Minimum myocardial enhancement mean signal
            [~,tmin]    = min(cm);
            data        = MRX.readImgData(tmin); data = single(data);
            MRX.updateContrastConc(tmin,ca,cm);
            [img,msk]   = MRX.mapTissueProps(data); img = single(img); msk = single(msk);
            img         = MRX.multiplyCoilMaps(img,sen);
            roi         = msk==MRX.Par.act.myoLA_act|msk==MRX.Par.act.myoLV_act;
            roi         = repmat(roi,[1,1,1,size(sen,4)]).*img;
            for k=1:MRX.Par.scan.coils, roik = roi(:,:,:,k); smin(k) = mean(roik(roik~=0)); end
            sumsqmin    = sqrt(sum(abs(smin).^2));
            
            % Calculate std dev based on CNR and contrast = mean(max)-mean(min)/CNR
            stdev       = 1/MRX.Par.scan.snr*(sumsqmax-sumsqmin); 
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
                
                % ----------------------------------------------------------------
                %   Select tissue type
                % ----------------------------------------------------------------
                switch char(tis(i))
                    
                    % ------------------------------------------------------------
                    %   Myocard
                    % ------------------------------------------------------------
                    case 'myoLV_act'
                        rho = MRX.Par.tissue.rhomuscle;
                        r1  = 1/MRX.Par.tissue.t1muscle+MRX.Par.contrast.cm.lv*MRX.Par.contrast.ry;
                    case 'myoRV_act'
                        rho = MRX.Par.tissue.rhomuscle;
                        r1  = 1/MRX.Par.tissue.t1muscle+MRX.Par.contrast.cm.rv*MRX.Par.contrast.ry;
                    case 'myoLA_act'
                        rho = MRX.Par.tissue.rhomuscle;
                        r1  = 1/MRX.Par.tissue.t1muscle+MRX.Par.contrast.cm.la*MRX.Par.contrast.ry;
                    case 'myoRA_act'
                        rho = MRX.Par.tissue.rhomuscle;
                        r1  = 1/MRX.Par.tissue.t1muscle+MRX.Par.contrast.cm.ra*MRX.Par.contrast.ry;
                    % ------------------------------------------------------------
                    %   Blood
                    % ------------------------------------------------------------
                    case 'bldplLV_act'
                        rho = MRX.Par.tissue.rhoblood;
                        r1  = 1/MRX.Par.tissue.t1blood+MRX.Par.contrast.ca.lv*MRX.Par.contrast.ry;
                    case {'bldplRV_act','art_activity','vein_activity'}
                        rho = MRX.Par.tissue.rhoblood;
                        r1  = 1/MRX.Par.tissue.t1blood+MRX.Par.contrast.ca.rv*MRX.Par.contrast.ry;
                    case 'bldplLA_act'
                        rho = MRX.Par.tissue.rhoblood;
                        r1  = 1/MRX.Par.tissue.t1blood+MRX.Par.contrast.ca.la*MRX.Par.contrast.ry;
                    case 'bldplRA_act'
                        rho = MRX.Par.tissue.rhoblood;
                        r1  = 1/MRX.Par.tissue.t1blood+MRX.Par.contrast.ca.ra*MRX.Par.contrast.ry;
                    % ------------------------------------------------------------
                    %   Body
                    % ------------------------------------------------------------
                    case {'body_activity','pericardium_activity'}
                        rho = MRX.Par.tissue.rhofat;
                        r1  = 1/MRX.Par.tissue.t1fat;
                    case 'muscle_activity'
                        rho = MRX.Par.tissue.rhomuscle;
                        r1  = 1/MRX.Par.tissue.t1muscle;
                    case 'liver_activity'
                        rho = MRX.Par.tissue.rholiver;
                        r1  = 1/MRX.Par.tissue.t1liver;
                    case {'rib_activity','cortical_bone_activity','spine_activity','bone_marrow_activity'}
                        rho = MRX.Par.tissue.rhobone;
                        r1  = 1/MRX.Par.tissue.t1bone;
                    otherwise
                        rho = 0;
                end
                
                % ----------------------------------------------------------------
                %   Signal model (Spoiled GRE)
                % ----------------------------------------------------------------
                a  = cos(MRX.Par.scan.flip)*exp(-MRX.Par.scan.trep*r1);
                b  = 1-exp(-MRX.Par.scan.trep*r1);
                n  = MRX.Par.scan.nky0;
                TD = MRX.Par.scan.tsat;
                sig = rho*( (1-exp(-TD*r1)).*a.^(n-1) + b.*(1 - a.^(n-1))./(1 - a) );
                
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
        % Update contrast concentration at specific time
        function updateContrastConc(MRX,t,ca,cm)
            % Get contrast concentration in right atrium (ra), right ventricle (rv),
            % left atrium (la) and left ventricle (lv) at specific time t for blood 
            % pool and myocardium.
            % This method is used during the loop over MRXCAT time frames while 
            % creating the phantom. 
            
            
            % ----------------------------------------------------------------
            %   Update contrast conc indices for different compartments
            %   Contrast arrives first in RA, then RV, LA and finally LV
            % ----------------------------------------------------------------
            ra = min(round((t+4)*length(ca)/MRX.Par.scan.frames)+1,length(ca));
            rv = min(round((t+3)*length(ca)/MRX.Par.scan.frames)+1,length(ca));
            la = min(round((t+0)*length(ca)/MRX.Par.scan.frames)+1,length(ca));
            lv = min(round((t-1)*length(ca)/MRX.Par.scan.frames)+1,length(ca));
            
            MRX.Par.contrast.ca.ra = ca(ra);            % RA blood pool
            MRX.Par.contrast.ca.rv = ca(rv);            % RV blood pool
            MRX.Par.contrast.ca.la = ca(la);            % LA blood pool
            MRX.Par.contrast.ca.lv = ca(lv);            % LV blood pool
            
            MRX.Par.contrast.cm.ra = cm(ra);            % RA myocardium
            MRX.Par.contrast.cm.rv = cm(rv);            % RV myocardium
            MRX.Par.contrast.cm.la = cm(la);            % LA myocardium
            MRX.Par.contrast.cm.lv = cm(lv);            % LV myocardium
        end
        
        
        %=========================================================================
		% Calc Fermi fcn. using alpha, beta and time info
        function h = fermiFunction(MRX, t)
            alpha = MRX.Par.contrast.falpha;
            beta  = MRX.Par.contrast.fbeta;
            h     = (1+beta)./(1+beta*exp(alpha*t));
        end
        
        
        %=========================================================================
		% Time-shift AIF by time dT
        function y = shiftAIF(MRX, t, aif)
            dT = MRX.Par.contrast.tshift;
            y  = interp1( t+dT, aif, t, 'pchip','extrap' );
        end
        
        
    end
    
    methods( Static )
        %=========================================================================
        % Convolve c with h along t and scale with factor q
        function c = convolve(t, c, h, q)
            dt = (t(end)-t(1)) / (length(t)-1);
            
            c_ = zeros(length(c));
            for i=1:length(c)
                for j =1:length(c)
                    if(i+1-j>0)
                        c_(i,j) = c(i+1-j);
                    end
                end
            end
            if size(h,1) == 1
                h=h';
            end
            c = q*c_*h*dt;
        end
        
        
 
    end
    
end

