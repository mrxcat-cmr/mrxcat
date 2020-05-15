%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [data, fname] = DisplayMRXCAT( fname, plot_result )

if nargin < 2 || ~exist(fname,'file')
    % get MRXCAT output (.cpx file)
    [filename,pathname] = uigetfile({'*.cpx'});
    fname = [pathname filename];
end
if nargin < 2
    plot_result = 1;
end

if ~isequal(fname(2),0) % check if cancel was pressed!
    %-------------------------------------------------------------------------
    % load MRXCAT parameter file (acquisition matrix, coils, ...) and update
    %-------------------------------------------------------------------------
    load( [fname(1:end-4) '_par.mat'] );
    
    
    %-------------------------------------------------------------------------
    % get MRXCAT data
    %-------------------------------------------------------------------------
    fid  = fopen( fname );
    data = fread(fid,inf,'float','l');
    fclose(fid);
    %-------------------------------------------------------------------------
    % get sensitivity maps
    %-------------------------------------------------------------------------
    fid  = fopen( [fname(1:end-4) '.sen'] );
    sen  = fread(fid,inf,'float','l');
    fclose(fid);
    
    %-------------------------------------------------------------------------
    % reformat to complex data (image and sensitivities)
    %-------------------------------------------------------------------------
    data = reshape( data, 2, []);
    data = data(1,:)+1i*data(2,:);
    data = reshape( data, Par.acq_matrix(1), Par.acq_matrix(2), Par.acq_matrix(3), Par.ncoils, []);
    data = permute( data, [1 2 3 5 4] );
    sen  = reshape( sen, 2, []);
    sen  = sen(1,:)+1i*sen(2,:);
    sen  = reshape( sen, Par.acq_matrix(1), Par.acq_matrix(2), Par.acq_matrix(3), Par.ncoils, []);
    sen  = permute( sen, [1 2 3 5 4] );
    
    %-------------------------------------------------------------------------
    % "coil-combine" data
    %-------------------------------------------------------------------------
    sos = sum(data./sen,5);
    clear sen;
    
    %-------------------------------------------------------------------------
    % Display 4-panel figure showing time frame, slices, coil images and (for perf) example signal-time curves
    %-------------------------------------------------------------------------
    if plot_result
        figure;
        displayMovie( permute(squeeze(abs(sos(:,:,round(end/2),:))),[2,1,3]), [2 2 1],1,0.2 );
        title('time frames');
        displayMovie( permute(squeeze(abs(sos(:,:,:,round(end/2)))),[2,1,3]), [2 2 2],1,0.3 );
        title('slices');
        displayMovie( permute(squeeze(abs(data(:,:,round(end/2),round(end/2),:))),[2,1,3]), [2 2 3],1, 1 );
        title('coil maps');
        if isfield(Par.contrast,'aif') %is perf scan
            [sa,sm,sa_ind,sm_ind] = extractSignalTimeCurves( sos, filename, pathname);
            subplot(2,2,4); hold all;
            plot(1:length(sa),abs(sa),1:length(sm),abs(sm),1:length(sa),abs(sa_ind),1:length(sm),abs(sm_ind),'LineWidth',2);
            title('mean and single-voxel AIF and MYO signal');
            xlabel('time frame [heart beats]');
            ylabel('signal intensity [a.u.]');
            axis tight;
        end
    end
else
    warndlg('No file loaded');
    data = [];
    fname = '';
end

end

%=========================================================================
%	Display data as a movie
%=========================================================================
function displayMovie(data,subplotno,loops,dt,scale)

    numberimages = numel(data)/(size(data,1)*size(data,2));
    data = reshape(data,size(data,1),size(data,2),numberimages);
    imageposno = 1;
    
    if nargin>1 && numel(subplotno)>=2,
        gridsizex = subplotno(1);
        gridsizey = subplotno(2);
        if numel(subplotno)>=3,
            imageposno = subplotno(3);
        end
    else
        gridsizex = 1;
        gridsizey = 1;
    end
    
    if nargin<3, 
        defaultLoops = 3; 
    else
        defaultLoops = loops;
    end
    
    if nargin<4, 
        defaultDt = 0.1; 
    else
        defaultDt = dt;
    end
     
    if nargin<5, 
        defaultScale = 1; 
    else
        defaultScale = 0;
        minScale = scale(1);
        maxScale = scale(2);
    end

    if isreal(data),
        RealData = 1;       % real data
        maxint = max(data(:));
        minint = min(data(:));
    elseif isreal(data*i),
        RealData = -1;      % imaginary data
        maxint = max(imag(data(:)));
        minint = min(imag(data(:))); 
    else
        RealData = 0;       % complex data
        maxint = max(abs(data(:)));
        minint = min(abs(data(:))); 
    end
    
    if ~defaultScale
        minint = minScale;
        maxint = maxScale;
    end

    if (maxint==minint),
        maxint=maxint+1;
        minint=minint-1;
    end

    for loop=1:defaultLoops
        for imageno=1:numberimages
            subplot(gridsizex,gridsizey,imageposno);
            colormap gray;
            if RealData==1,         % real data
                imagesc(transpose(      data(:,:,imageno)  ), [minint, maxint]);
            elseif RealData==-1,    % imaginary data
                imagesc(transpose( imag(data(:,:,imageno)) ), [minint, maxint]);
            else                    % complex data
                imagesc(transpose(  abs(data(:,:,imageno)) ), [minint, maxint]);
            end
            graphtitle = sprintf('%2d / %2d',imageno,numberimages);
            title(graphtitle);
            axis image;
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            pause(defaultDt);
        end
    end
    
    imageno=bitshift(numberimages,-1)+1;
    subplot(gridsizex,gridsizey,imageposno);

    if RealData==1,         % real data
        imagesc(transpose(      data(:,:,imageno)  ), [minint, maxint]);
    elseif RealData==-1,    % imaginary data
        imagesc(transpose( imag(data(:,:,imageno)) ), [minint, maxint]);
    else                    % complex data
        imagesc(transpose(  abs(data(:,:,imageno)) ), [minint, maxint]);
    end
    graphtitle = sprintf('%2d / %2d',imageno,numberimages);
    title(graphtitle);
    axis image;
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');    
end


%=========================================================================
%	Extract Signal-time curves from MRXCAT using mask from .msk file
%   as segmentation
%=========================================================================
function [sa,sm,sa_ind,sm_ind] = extractSignalTimeCurves( data, filename, pathname )
    try
        filename = [filename(1:regexp(filename,'bh')+1) '.msk'];
        fmask = [pathname filename];
        fid = fopen( fmask );
        msk = fread(fid,inf,'uint8','l');
        siz = size(data);
        msk=reshape(msk,siz(1:4));
        fclose(fid);
    catch
        [filename,pathname] = uigetfile({'.msk'},'Select Mask file');
        fmask = [pathname filename];
        fid = fopen( fmask );
        msk = fread(fid,inf,'uint8','l');
        msk=reshape(msk,size(data));
        fclose(fid);
    end
    mym = msk==1;
    lvm = msk==5;
    datacc = sqrt( sum( abs( data ).^2, 5 ) )/sqrt(size(data,5));
    mym = datacc(:,:,round(end/2),:).* mym(:,:,round(end/2),:);
    lvm = datacc.* lvm(:,:,:,:);
    ind = find(squeeze(mym(:,:,1,1))>0);
    indl= find(squeeze(lvm(:,:,:,1))>0);
    for k=1:size(mym,4)
        temp = mym(:,:,1,k);
        sm(k) = mean(temp(ind));
        sm_ind(k) = temp(ind(1));
        temp = lvm(:,:,:,k);
        sa(k) = mean(temp(indl));
        sa_ind(k) = temp(indl(1));
    end
end
%=========================================================================    