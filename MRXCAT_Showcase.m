classdef MRXCAT_Showcase < handle
    
    properties
        mrx             % MRXCAT object
        img             % image attributes
        filenames       % stores different filenames (e.g. XCAT .bin filename, MRXCAT .cpx filename)
        
        gui             % GUI handles
    end
    
    methods
        % Constructor
        function MRD = MRXCAT_Showcase()
            MRD.layout;
            
            MRD.initUpdateFields;
        end
        
        function layout( MRD )
            % Initialize GUI structure (sorted)
            MRD.gui.axes     = struct;
            MRD.gui.button   = struct;
            MRD.gui.checkbox = struct;
            MRD.gui.edit     = struct;
            MRD.gui.layout   = struct;
            MRD.gui.menu     = struct;
            MRD.gui.panel    = struct;
            MRD.gui.popup    = struct;
            MRD.gui.slider   = struct;
            MRD.gui.string   = struct;
            
            % Header: Define standard styles
            editStyle       = {'Style', 'edit', 'BackgroundColor', 'w', 'Enable', 'on'};               % edit box style
            sliderStyle     = {'Style', 'slider', 'BackgroundColor', [0.9 0.9 0.9], 'Enable', 'off'};  % slider style
            buttonStyle     = {'Enable', 'on'};                                                        % text button style
            popupStyle      = {'Style', 'popupmenu', 'BackgroundColor', 'w', 'Enable', 'on'};          % popup style
            checkboxStyle   = {'Style', 'checkbox', 'Enable', 'on'};                                   % checkbox style
            textStyle       = {'Style', 'text', 'HorizontalAlignment', 'left'};                        % static text style
            
            % Create GUI window
            screensz        = get(0,'ScreenSize');
            MRD.gui.fig      = figure( 'Position', [200, -screensz(2)+140, 1100, 800], 'MenuBar', 'none', ...
                'Name', 'MRXCAT Showcase GUI', 'NumberTitle', 'off' );
            
            MRD.gui.layout.top = uiextras.HBox( 'Parent', MRD.gui.fig );
            
            % Create viewing display and action buttons
            MRD.gui.layout.lay21 = uiextras.VBox( 'Parent', MRD.gui.layout.top );
            MRD.gui.axes.img     = axes( 'Parent', MRD.gui.layout.lay21, 'Position', [0 0 1 1] );
            
            % Create 3 sliders (time, slices, coils)
            MRD.gui.layout.lay22     = uiextras.VBox( 'Parent', MRD.gui.layout.top );
            MRD.gui.string.slice     = uicontrol( 'Parent', MRD.gui.layout.lay22, 'String', 'z', textStyle{:} );
            MRD.gui.slider.slice     = uicontrol( 'Parent', MRD.gui.layout.lay22, 'Callback', @MRD.updateSlider, sliderStyle{:} );
            MRD.gui.string.time      = uicontrol( 'Parent', MRD.gui.layout.lay22, 'String', 't', textStyle{:} );
            MRD.gui.slider.time      = uicontrol( 'Parent', MRD.gui.layout.lay22, 'Callback', @MRD.updateSlider, sliderStyle{:} );
            MRD.gui.string.coils     = uicontrol( 'Parent', MRD.gui.layout.lay22, 'String', 'c', textStyle{:} );
            MRD.gui.slider.coils     = uicontrol( 'Parent', MRD.gui.layout.lay22, 'Callback', @MRD.updateSlider, sliderStyle{:} );
            MRD.gui.layout.lay22.Sizes = [14, -1, 14, -1, 14, -1];

            % Create options panel
            MRD.gui.layout.lay23      = uiextras.VBox( 'Parent', MRD.gui.layout.top );
            
            % general options
            MRD.gui.panel.optsMRXCAT  = uipanel( 'Parent', MRD.gui.layout.lay23, 'Title', 'General MRXCAT options' );
            MRD.gui.layout.lay231     = uiextras.HBox( 'Parent', MRD.gui.panel.optsMRXCAT, 'Spacing', 5, 'Padding', 3 );
            MRD.gui.layout.lay2311    = uiextras.VBox( 'Parent', MRD.gui.layout.lay231, 'Spacing', 1, 'Padding', 1 );
            MRD.gui.string.mrxcatType = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'MRXCAT type', textStyle{:} );
            MRD.gui.string.tr         = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'TR [ms]', textStyle{:} );
            MRD.gui.string.te         = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'TE [ms]', textStyle{:} );
            MRD.gui.string.flipAngle  = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'Flip angle [deg]', textStyle{:} );
            MRD.gui.string.snr        = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'SNR/CNR', textStyle{:} );
            MRD.gui.string.coilsNo    = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', '# coils', textStyle{:} );
            MRD.gui.string.trajectory = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'Trajectory', textStyle{:} );
            MRD.gui.string.undersamp  = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'Undersampling factor', textStyle{:} );
            MRD.gui.checkbox.lowPass  = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'Low-pass Filter', 'Value', 1, 'Callback', @MRD.editParameters, checkboxStyle{:} );
            MRD.gui.string.lowPassStr = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'Filter strength', textStyle{:} );
            MRD.gui.string.bbox_low   = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'Box-Crop low (x|y|z)', textStyle{:} );
            MRD.gui.string.bbox_high  = uicontrol( 'Parent', MRD.gui.layout.lay2311, 'String', 'Box-Crop high (x|y|z)', textStyle{:} );
            MRD.gui.layout.lay2312    = uiextras.VBox( 'Parent', MRD.gui.layout.lay231, 'Spacing', 1, 'Padding', 1 );
            MRD.gui.popup.mrxcatType  = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', {'Cine';'Perfusion'}, 'Callback', @MRD.editParameters, popupStyle{:} );
            MRD.gui.edit.tr           = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', '2.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.te           = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', '1.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.flipAngle    = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', '50', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.snr          = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', '30', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.coils        = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', '8', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.popup.trajectory  = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', {'Cartesian';'Radial';'GoldenAngle'}, 'Callback', @MRD.editParameters, popupStyle{:} );
            MRD.gui.edit.undersamp    = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', '1', 'Callback', @MRD.editParameters, editStyle{:} );
            uiextras.Empty( 'Parent', MRD.gui.layout.lay2312 );
            MRD.gui.edit.lowPassStr   = uicontrol( 'Parent', MRD.gui.layout.lay2312, 'String', '1.2', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.layout.lay23121   = uiextras.HBox( 'Parent', MRD.gui.layout.lay2312 );
            MRD.gui.edit.bbox_x_low   = uicontrol( 'Parent', MRD.gui.layout.lay23121, 'String', '0.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.bbox_y_low   = uicontrol( 'Parent', MRD.gui.layout.lay23121, 'String', '0.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.bbox_z_low   = uicontrol( 'Parent', MRD.gui.layout.lay23121, 'String', '0.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.layout.lay23122   = uiextras.HBox( 'Parent', MRD.gui.layout.lay2312 );
            MRD.gui.edit.bbox_x_high  = uicontrol( 'Parent', MRD.gui.layout.lay23122, 'String', '1.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.bbox_y_high  = uicontrol( 'Parent', MRD.gui.layout.lay23122, 'String', '1.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.bbox_z_high  = uicontrol( 'Parent', MRD.gui.layout.lay23122, 'String', '1.0', 'Callback', @MRD.editParameters, editStyle{:} );
            
            % Cine options
            MRD.gui.panel.optsCine    = uipanel( 'Parent', MRD.gui.layout.lay23, 'Title', 'Cine options' );
            MRD.gui.layout.lay232     = uiextras.HBox( 'Parent', MRD.gui.panel.optsCine, 'Spacing', 5, 'Padding', 3 );
            MRD.gui.string.segments   = uicontrol( 'Parent', MRD.gui.layout.lay232, 'String', '# segments (cine)', textStyle{:} );
            MRD.gui.edit.segments     = uicontrol( 'Parent', MRD.gui.layout.lay232, 'String', '15', 'Callback', @MRD.editParameters, editStyle{:} );

            % Perfusion options
            MRD.gui.panel.optsPerf    = uipanel( 'Parent', MRD.gui.layout.lay23, 'Title', 'Perfusion options' );
            MRD.gui.layout.lay233     = uiextras.HBox( 'Parent', MRD.gui.panel.optsPerf, 'Spacing', 5, 'Padding', 3 );
            MRD.gui.layout.lay2331    = uiextras.VBox( 'Parent', MRD.gui.layout.lay233, 'Spacing', 1, 'Padding', 1 );
            MRD.gui.string.restStress = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Rest / Stress', textStyle{:} );
            MRD.gui.string.restMBF    = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Rest MBF [ml/g/min]', textStyle{:} );
            MRD.gui.string.stressMBF  = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Stress MBF [ml/g/min]', textStyle{:} );
            MRD.gui.string.FermiAlpha = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Fermi alpha', textStyle{:} );
            MRD.gui.string.FermiBeta  = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Fermi beta', textStyle{:} );
            MRD.gui.string.tShift     = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Time shift [s]', textStyle{:} );
            MRD.gui.string.RRCycleDur = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'RR duration [s]', textStyle{:} );
            MRD.gui.string.tSat       = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'SAT delay [ms]', textStyle{:} );
            MRD.gui.string.prof2k0    = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Profs -> k0', textStyle{:} );
            MRD.gui.string.dynamics   = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Time Frames', textStyle{:} );
            MRD.gui.string.respMotion = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Respiration', textStyle{:} );
            MRD.gui.string.dose       = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Gd dose [mmol/kg b.w.]', textStyle{:} );
            MRD.gui.string.relaxivity = uicontrol( 'Parent', MRD.gui.layout.lay2331, 'String', 'Gd relaxivity [l/mmol*s]', textStyle{:} );
            MRD.gui.layout.lay2332    = uiextras.VBox( 'Parent', MRD.gui.layout.lay233, 'Spacing', 1, 'Padding', 1 );
            MRD.gui.popup.stressRest  = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', {'Rest';'Stress'}, 'Callback', @MRD.editParameters, popupStyle{:} );
            MRD.gui.edit.restMBF      = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '1.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.stressMBF    = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '3.5', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.FermiAlpha   = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '0.15', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.FermiBeta    = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '0.15', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.tShift       = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '3', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.RRCycleDur   = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '1.0', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.tSat         = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '150', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.prof2k0      = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '40', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.dynamics     = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '30', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.popup.respiration = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', {'Off';'On'}, 'Callback', @MRD.editParameters, popupStyle{:} );
            MRD.gui.edit.dose         = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '0.05', 'Callback', @MRD.editParameters, editStyle{:} );
            MRD.gui.edit.relaxivity   = uicontrol( 'Parent', MRD.gui.layout.lay2332, 'String', '5.6', 'Callback', @MRD.editParameters, editStyle{:} );
            
            MRD.gui.panel.runLoad     = uipanel( 'Parent', MRD.gui.layout.lay23, 'Title', 'Create or Load' );
            MRD.gui.layout.lay234     = uiextras.HBox( 'Parent', MRD.gui.panel.runLoad );
            MRD.gui.button.runMRXCAT  = uicontrol( 'Parent', MRD.gui.layout.lay234, 'String', 'Create Dataset', 'Callback', @MRD.runMRXCAT, buttonStyle{:} );
            MRD.gui.button.loadMRXCAT = uicontrol( 'Parent', MRD.gui.layout.lay234, 'String', 'Load Dataset', 'Callback', @MRD.showMRXCAT, buttonStyle{:} );
            MRD.gui.layout.lay23.Sizes = [-12, 44,-13,60];
            
            % Resize whole GUI
            MRD.gui.layout.top.Sizes   = [-2, 24, 260]; % fixed width for slider (col 2), refresh buttons (col 4), fitting panel (col 5)
            
        end
        
        function runMRXCAT( MRD, hObject, eventdata )
            MRD.initUpdateFields;
            par = MRD.getParametersFromGUI;
            mrxcatTypes = get( MRD.gui.popup.mrxcatType, 'String' );
%             try
                if strcmpi( mrxcatTypes{get( MRD.gui.popup.mrxcatType, 'Value' )}, 'Cine' )
                    MRD.mrx = MRXCAT_CMR_CINE( '', 'demo_gui', par );
                else
                    MRD.mrx = MRXCAT_CMR_PERF( '', 'demo_gui', par );
                end
                MRD.filenames.cpx = MRD.mrx.Par.rec_file_name;
                MRD.filenames.bin = MRD.mrx.Filename;
                
                MRD.showMRXCAT;
%             catch
%                 warning('Something went wrong creating the MRXCAT phantom! Maybe you just didn''t select a .bin file!');
%             end
        end
        
        function showMRXCAT( MRD, hObject, eventdata )
            if nargin>1 && hObject == MRD.gui.button.loadMRXCAT
                [MRD.img.data, MRD.filenames.cpx] = DisplayMRXCAT( '', 0 );
            else
                [MRD.img.data, MRD.filenames.cpx] = DisplayMRXCAT( MRD.filenames.cpx, 0 );
            end
            if ~isempty( MRD.filenames.cpx )
                MRD.sendParameters2GUI;
                MRD.undersampleIfCartesian;
                MRD.updateImageDisplay( MRD.gui.axes.img );
                MRD.initUpdateFields;
                MRD.updateSlider( MRD.gui.slider.slice );
                MRD.updateSlider( MRD.gui.slider.time );
                MRD.updateSlider( MRD.gui.slider.coils );
            end
        end
        
        function editParameters( MRD, hObject, eventdata )
            switch hObject
                case MRD.gui.popup.mrxcatType
                    mrxcatTypes = get( MRD.gui.popup.mrxcatType, 'String' );
                    if strcmpi( mrxcatTypes{get( MRD.gui.popup.mrxcatType, 'Value' )}, 'Cine' )
                        set(MRD.gui.layout.lay232,'Enable','on');
                        set(MRD.gui.layout.lay233,'Enable','off');
                        set(MRD.gui.popup.trajectory,'String',{'Cartesian';'Radial';'GoldenAngle'});
                    else
                        set(MRD.gui.layout.lay232,'Enable','off');
                        set(MRD.gui.layout.lay233,'Enable','on');
                        set(MRD.gui.popup.trajectory,'String','Cartesian');
                    end                        
                case MRD.gui.checkbox.lowPass
                    if get(MRD.gui.checkbox.lowPass, 'Value')
                        set( MRD.gui.edit.lowPassStr, 'Enable', 'on' );
                    else
                        set( MRD.gui.edit.lowPassStr, 'Enable', 'off' );
                    end
            end
        end
        
        function par = getParametersFromGUI( MRD )
            
            % General options
            par.scan.trep               = str2double(get(MRD.gui.edit.tr,'String'));
            par.scan.te                 = str2double(get(MRD.gui.edit.te,'String'));
            par.scan.flip               = str2double(get(MRD.gui.edit.flipAngle,'String'))*pi/180;
            par.scan.snr                = str2double(get(MRD.gui.edit.snr,'String'));
            par.scan.coils              = str2double(get(MRD.gui.edit.coils,'String'));
            trajs                       = get( MRD.gui.popup.trajectory, 'String');
            if ~iscell(trajs), trajs    = {trajs}; end
            par.scan.trajectory         = trajs{get(MRD.gui.popup.trajectory,'Value')};
            par.scan.undersample        = str2double(get(MRD.gui.edit.undersamp,'String'));
            par.scan.lowpass            = get(MRD.gui.checkbox.lowPass,'Value');
            par.scan.lowpass_str        = str2double(get(MRD.gui.edit.lowPassStr,'String'));
            par.scan.bbox               = [ ... 
                str2double(get(MRD.gui.edit.bbox_x_low,'String')),str2double(get(MRD.gui.edit.bbox_x_high,'String')); ...
                str2double(get(MRD.gui.edit.bbox_y_low,'String')),str2double(get(MRD.gui.edit.bbox_y_high,'String')); ... 
                str2double(get(MRD.gui.edit.bbox_z_low,'String')),str2double(get(MRD.gui.edit.bbox_z_high,'String'));];
            
            % Cine parameters
            par.scan.segments           = str2double(get(MRD.gui.edit.segments,'String'));
            
            % Perfusion parameters
            par.contrast.rs             = get(MRD.gui.popup.stressRest,'Value');
            par.contrast.qr             = str2double(get(MRD.gui.edit.restMBF,'String'))/60;
            par.contrast.qs             = str2double(get(MRD.gui.edit.stressMBF,'String'))/60;
            par.contrast.falpha         = str2double(get(MRD.gui.edit.FermiAlpha,'String'));
            par.contrast.fbeta          = str2double(get(MRD.gui.edit.FermiBeta,'String'));
            par.contrast.tshift         = str2double(get(MRD.gui.edit.tShift,'String'));
            par.scan.trrc               = str2double(get(MRD.gui.edit.RRCycleDur,'String'));
            par.scan.tsat               = str2double(get(MRD.gui.edit.tSat,'String'));
            par.scan.nky0               = str2double(get(MRD.gui.edit.prof2k0,'String'));
            if str2double(get(MRD.gui.edit.dynamics,'String'))>0 && get(MRD.gui.popup.mrxcatType,'Value')>1 %&& ~RespMotion % only overwrite par.scan.frames, if Frames ~= 0
                par.scan.frames         = str2double(get(MRD.gui.edit.dynamics,'String'));
            end
            par.scan.resp               = get(MRD.gui.popup.respiration,'Value')-1;
            par.contrast.dose           = str2double(get(MRD.gui.edit.dose,'String'));
            par.contrast.ry             = str2double(get(MRD.gui.edit.relaxivity,'String'))/1000;   % [l/(mmol*ms)]
            
        end
        
        function sendParameters2GUI( MRD )
            load( [MRD.filenames.cpx(1:end-4) '_par.mat'] );
            
            % Cine or Perf
            if isfield(Par.scan,'phases')
                set( MRD.gui.popup.mrxcatType, 'Value', 1 );
            else
                set( MRD.gui.popup.mrxcatType, 'Value', 2 );
            end
            MRD.editParameters( MRD.gui.popup.mrxcatType );
            
            if isfield( Par.scan, 'segments')
                % Cine parameters
                set( MRD.gui.edit.segments, 'String', num2str(Par.scan.segments));
                set( MRD.gui.edit.te, 'String', num2str(Par.scan.te) );
            else
                % Perf parameters
                set( MRD.gui.popup.stressRest,  'Value',  Par.contrast.rs);
                set( MRD.gui.edit.restMBF,      'String', num2str(Par.contrast.qr*60));
                set( MRD.gui.edit.stressMBF,    'String', num2str(Par.contrast.qs*60));
                set( MRD.gui.edit.FermiAlpha,   'String', num2str(Par.contrast.falpha));
                set( MRD.gui.edit.FermiBeta,    'String', num2str(Par.contrast.fbeta));
                set( MRD.gui.edit.tShift,       'String', num2str(Par.contrast.tshift));
                set( MRD.gui.edit.RRCycleDur,   'String', num2str(Par.scan.trrc));
                set( MRD.gui.edit.tSat,         'String', num2str(Par.scan.tsat));
                set( MRD.gui.edit.prof2k0,      'String', num2str(Par.scan.nky0));
                if str2double(get(MRD.gui.edit.dynamics,'String'))>0 && get(MRD.gui.popup.mrxcatType,'Value')>1 %&& ~RespMotion % only overwrite par.scan.frames, if Frames ~= 0
                    par.scan.frames         = str2double(get(MRD.gui.edit.dynamics,'String'));
                end
                set( MRD.gui.popup.respiration, 'Value',  Par.scan.resp+1 );
                set( MRD.gui.edit.dose,         'String', num2str(Par.contrast.dose));
                set( MRD.gui.edit.relaxivity,   'String', num2str(Par.contrast.ry*1000));   % [l/(mmol*ms)]
            end
            % General options
            set( MRD.gui.edit.tr,           'String', num2str(Par.scan.trep));
            set( MRD.gui.edit.flipAngle,    'String', num2str(Par.scan.flip*180/pi));
            set( MRD.gui.edit.snr,          'String', num2str(Par.scan.snr) );
            set( MRD.gui.edit.coils,        'String', num2str(Par.ncoils) );
            trajs = get( MRD.gui.popup.trajectory, 'String');
            traj_ind = find(strcmpi(trajs,Par.scan.trajectory));
            set( MRD.gui.popup.trajectory,  'Value', traj_ind );
            if ~strcmpi( Par.scan.trajectory, 'Cartesian' )
                set(MRD.gui.edit.undersamp, 'String', num2str(Par.scan.undersample) );
            end
            set( MRD.gui.checkbox.lowPass,  'Value',  Par.scan.lowpass );
            set( MRD.gui.edit.lowPassStr,   'String', num2str(Par.scan.lowpass_str) );
            MRD.editParameters( MRD.gui.checkbox.lowPass );
            set( MRD.gui.edit.bbox_x_low,   'String', num2str(Par.scan.bbox(1,1)));
            set( MRD.gui.edit.bbox_x_high,  'String', num2str(Par.scan.bbox(1,2)));
            set( MRD.gui.edit.bbox_y_low,   'String', num2str(Par.scan.bbox(2,1)));
            set( MRD.gui.edit.bbox_y_high,  'String', num2str(Par.scan.bbox(2,2)));
            set( MRD.gui.edit.bbox_z_low,   'String', num2str(Par.scan.bbox(3,1)));
            set( MRD.gui.edit.bbox_z_high,  'String', num2str(Par.scan.bbox(3,2)));
            
        end
        
        % update image display on current axes
        function updateImageDisplay( MRD, cur_axes, eventdata )
            if nargin>1
                axes( cur_axes );
            else
                cur_axes = MRD.gui.axes.img;
            end
            cla( cur_axes );
            try
                imagesc(abs(MRD.img.data( :,:, MRD.img.slicesliderpos, MRD.img.dynsliderpos, MRD.img.coilsliderpos )));
            catch
                MRD.img.slicesliderpos = 1;
                MRD.img.dynsliderpos = 1;
                MRD.img.coilsliderpos = 1;
                imagesc(abs(MRD.img.data( :,:, MRD.img.slicesliderpos, MRD.img.dynsliderpos, MRD.img.coilsliderpos )));
            end
            tstr = sprintf('MRXCAT image, slice %2d, dyn %2d, coil %2d', MRD.img.slicesliderpos, MRD.img.dynsliderpos, MRD.img.coilsliderpos );
            title(tstr);
            axis off; axis equal; axis tight;
            colormap gray;
        end

        % slider callback
        function updateSlider( MRD, hObject, eventdata, slideTo ) % use eventdata to indicate where to slide to (if desired)
            % get previous slider value, min, max
            pcv  = get(hObject,'Value');
            pmin = get(hObject,'Min');
            pmax = get(hObject,'Max');
            
            % calculate new values
            cmin = 1;
            switch hObject
                case MRD.gui.slider.slice
                    cmax = MRD.img.zdim; 
                case MRD.gui.slider.time
                    cmax = MRD.img.tdim;
                case MRD.gui.slider.coils
                    cmax = MRD.img.cdim;
            end
            if cmax>1
                set(hObject,'Enable','on');
                if nargin>3
                    ccv = slideTo;
                else
                    ccv  = max(min(pcv*(cmax-cmin)/(pmax-pmin),cmax),cmin);
                end
                
                % set new values
                set(hObject,'Min',cmin,'Max',cmax,'SliderStep',[1/(cmax-cmin) 3/(cmax-cmin)],'Value',ccv);
                
                % update slider positions and update image display
                switch hObject
                    case MRD.gui.slider.slice
                        MRD.img.slicesliderpos = round(ccv);
                        set(MRD.gui.string.slice,'String',sprintf('z:%d',MRD.img.slicesliderpos));
                    case MRD.gui.slider.time
                        MRD.img.dynsliderpos = round(ccv);
                        set(MRD.gui.string.time,'String',sprintf('t:%d',MRD.img.dynsliderpos));
                    case MRD.gui.slider.coils
                        MRD.img.coilsliderpos = round(ccv);
                        set(MRD.gui.string.coils,'String',sprintf('c:%d',MRD.img.coilsliderpos));                        
                end
                MRD.updateImageDisplay;
            else
                switch hObject
                    case MRD.gui.slider.slice
                        MRD.img.slicesliderpos = 1;
                    case MRD.gui.slider.time
                        MRD.img.dynsliderpos = 1;
                    case MRD.gui.slider.coils
                        MRD.img.coilsliderpos = 1;
                end
                set(hObject,'Enable','off');
            end

        end

        % initialize MRXCAT Demo fields
        function initUpdateFields( MRD )
            
            if ~isfield( MRD.img, 'data' )
                MRD.img.data = '';
            end
            
            if ~isfield( MRD.img, 'slicesliderpos' ) %check if already initialized
                % initialize slider positions
                MRD.img.slicesliderpos  = 1;
                MRD.img.dynsliderpos    = 1;
                MRD.img.coilsliderpos   = 1;
            end
            if ~isfield( MRD.filenames, 'cpx' )
                MRD.filenames.cpx = '';
            end
            if ~isempty(MRD.img.data)
                MRD.img.xdim            = size(MRD.img.data,1);
                MRD.img.ydim            = size(MRD.img.data,2);
                MRD.img.zdim            = size(MRD.img.data,3);
                MRD.img.tdim            = size(MRD.img.data,4);
                MRD.img.cdim            = size(MRD.img.data,5);
                MRD.img.slicesliderpos  = round( MRD.img.zdim/2 );
                MRD.img.dynsliderpos    = round( MRD.img.tdim/2 );
                MRD.img.coilsliderpos   = round( MRD.img.cdim/2 );
            else
                MRD.img.xdim            = 1;
                MRD.img.ydim            = 1;
                MRD.img.zdim            = 1;
                MRD.img.tdim            = 1;
                MRD.img.cdim            = 1;
                MRD.img.slicesliderpos  = 1;
                MRD.img.dynsliderpos    = 1;
                MRD.img.coilsliderpos   = 1;
            end
                
        end
        
        function undersampleIfCartesian( MRD )
            strs=get(MRD.gui.popup.trajectory,'String');
            if( ( iscell(strs) && strcmpi( strs{get(MRD.gui.popup.trajectory,'Value')}, 'Cartesian' ) ) || ...
                strcmpi( strs, 'Cartesian' ) )
                R    = str2num( get( MRD.gui.edit.undersamp, 'String') );
                temp = MRXCAT.i2k(MRD.img.data,[1 2 3]);
                [a,b,c,d]=ind2sub( size(temp), find(max(abs(temp(:)))==abs(temp(:))));
                ind=1:R:size(temp,2);
                while ~any(ind==b)
                    ind = ind+1;
                end
                ind  = ind(ind<=size(temp,2));
                %temp = temp(:,ind,:,:,:,:);
                tempz = zeros(size(temp));
                for i=1:length(ind)
                    tempz(:,ind(i),:,:,:,:) = temp(:,ind(i),:,:,:,:);
                end
                MRD.img.data = MRXCAT.k2i(tempz,[1 2 3]);
            end
        end
    end
    
end