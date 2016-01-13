function h = rastermap(GridSpecs,varargin)
% Display GeoEAS attribute raster in a Cartesian frame (up to 2D)
%
%% DESCRIPTION: rastermap.m
% Function to plot a GeoEAS-formatted attribute raster in a Cartesian
% coordinate system up to 3D. Raster origin is the lower left corner.
% Optionally, only the outline of a regular or irregular raster can be
% displayed.
% NOTE: Rasters can be regular, e.g., in 2D with rectangular pixels of 
%       fixed size, or irregular with variable pixel size, e.g., in 
%       2D rectangular pixels but with varying size
% NOTE: Pixels with missing values, flagged as NaNs, are color-coded white
% NOTE: This function uses imagesc for displaying regular rasters, and
%       pcolor for displaying irregular rasters
% NOTE: Accounting for NaNs using alphaData & the default OpenGL renderer
%       handles transparency, but distorts vector graphics.
%       set(gcf,'Renderer','Painters'); corrects the problem in GRAYSCALE
%       For COLOR, Renderer Painters does not handle transparency.
%       For COLOR, the default OpenGL rendered is used here;
%       this produces fine results on the screen, but not when
%       printing to an eps file; to overcome this print in GRAYSCALE
% NOTE: To convert COLOR image to GRAYSCALE outside this function, use:
%       colormap(1-colormap(gray)); set(gcf,'Renderer','Painters');
% NOTE: Exporting to an .eps file: print -deps (or -depsc) file.eps
%
%% SYNTAX:
%    rastermap(GridSpecs);
%    rastermap(GridSpecs,RasterIn,icolV,vType);
%    rastermap(GridSpecs,RasterIn,icolV,vType,colorLabels);
%    rastermap(GridSpecs,RasterIn,icolV,vType,colorLabels,plotPars);
%    h = rastermap(...);
%
%% INPUTS:
%    GridSpecs =            For REGULAR RASTERS: 
%                (nDim x 3) array with regular raster specifications:
%                e.g., in 2D: [nx xmin xsize ;  ny ymin ysize]
%                      nx    = # of pixels in x-direction
%                      xmin  = x-coordinate of lower-left pixel CENTER
%                      xsize = size of pixel in x-direction
%                      ny    = # of pixels in y-direction
%                      ymin  = y-coordinate of lower-left pixel CENTER
%                      ysize = size of pixel in y-direction
%                           For IRREGULAR RASTERS: 
%                (nDim x 3) CELL array with each cell containing the 
%                coordinates of pixel CENTERS along each direction;
%                coordinates must be monotonically increasing, e.g.,
%                in 2D: {0:100,0:50} or {[0:1:50 55:5:100],0:50}
%                NOTE: Dimensionality is determined by this array
%                NOTE: rastermap(GridSpecs) displays only raster outline
%    RasterIn  = (nPix x nCol) with input attribute rasters in GeoEAS
%                format; each column contains a different attribute raster
%                and all rasters have the same grid specifications.
%                NOTE: RasterIn can hold 1D, 2D, or 3D, attribute rasters 
%                      with missing values flagged as NaNs.
%    icolV     = scalar with column # (in RasterIn) of variable to display
%    vType      = string with flag for continuous ('cont')
%                 or categorical ('cat') attribute raster
%                 NOTE: This affects the color map used for each case
% colorLabels = Optional array with color limits OR category labels:
%                 For the continuous case (vType = 'cont'), 
%                       this is a (1 x 2) array [cmin cmax] with:  
%                       cmin = min color value for plotting
%                       cmax = max color value for plotting
%                 For the categorical case (vType = 'cat'),
%                       this is a (1 x nClass) cell array with strings for
%                       class labels, e.g., {'A','B','C'}
%                 NOTE: This could be an EMPTY array, in which case default
%                       color limits (min/max attribute values) and
%                       default class labels ('1','2',...) are used
%                 NOTE: For the categorical case, this array is used only 
%                       if a colorbar is requested (see below)
%  plotPars   = Optional (1 x 2) cell array {igray ibar} with flags
%                 for gray scale and bar legend to use for plotting:
%                       igray = 'color' -> color scale
%                       igray = 'gray' -> gray scale
%                       ibar  = 'nobar' -> no color/gray scale bar
%                       ibar  = 'vbar' -> vertical color/gray scale bar
%                       ibar  = 'hbar' -> horizontal color/gray scale bar
%                 NOTE: This could be an EMPTY array, in which case default
%                       ploting parameters, color & no bar, are used
%
%% OUTPUTS:
%    a MATLAB graph 
%  + handles to axes, surface plots and colorbar (if requested)
%
%% CUSTOM FUNCTIONS CALLED:
%    geoeas2matlab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%    rastermap(GridSpecs);
%    rastermap(GridSpecs,RasterIn,icolV,vType);
%    rastermap(GridSpecs,RasterIn,icolV,vType,colorLabels);
%    rastermap(GridSpecs,RasterIn,icolV,vType,colorLabels,plotPars);
%    h = rastermap(...);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2005                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixed parameters
VERS = 1;
NARGMIN = 1;
NARGDAT  = NARGMIN + 3;
NARGCOL  = NARGDAT + 1;
NARGPLOT = NARGCOL + 1;
YTOP1D = 0.05;         % Proportion of total x-length in 1D case
global MAXCLASS;
MAXCLASS = 10;         % max # of class labels to display

%% Error checking
if ~(   nargin == NARGMIN || nargin == NARGDAT || nargin == NARGCOL ...
     || nargin == NARGPLOT)
   tit=['You need ',num2str(NARGMIN),', ',num2str(NARGDAT),', ',...
        num2str(NARGCOL),' or ',num2str(NARGPLOT),...
        'input arguments for this function'];
   disp(tit);
   error('Check # of input arguments');
end
[colorLabels,plotPars] = checkerrors(GridSpecs,varargin{:});
% Note: plotPars = [igray ibar]
%       igray = 0, 1    -> color, gray scale 
%       ibar  = 0, 1, 2 -> nobar, vbar, hbar

%% Get grid specifications, and populate axes
if ~iscell(GridSpecs);  % Regular raster input
    iReg = 1;
    nDim = size(GridSpecs,1);
    nx = GridSpecs(1,1); xmin = GridSpecs(1,2); xsize = GridSpecs(1,3);
    x = (xmin:xsize:xmin+(nx-1)*xsize);
    ngTot = nx;
    if nDim > 1
        ny = GridSpecs(2,1); ymin = GridSpecs(2,2); ysize = GridSpecs(2,3);
        y = (ymin:ysize:ymin+(ny-1)*ysize);
        ngTot = nx*ny;
    end
    tit = 'Regular raster';
else                    % Irregular raster input
    iReg = 0;
    nDim = length(GridSpecs);
    x = GridSpecs{1};
    x = x(:)';
    nx = numel(x);
    ngTot = nx;
    if nDim > 1
        y = GridSpecs{2};
        y = y(:)';
        ny = numel(y);
        ngTot = nx*ny;
    end
    tit = 'Irregular raster';
end

%% Extract column to display
if nargin == NARGMIN
    iEmpty = 1;
    Data = nan(ngTot,1);
    iReg = 0;  % Use pcolor to display NaNs
    vType = 0;  % Empty array is treated as continuous
else
    iEmpty = 0;
    RasterIn = varargin{1};
    icolV = varargin{2};
    vType = varargin{3};
    if strcmp(vType,'cat')
        vType = 1;
    else
        vType = 0;
    end
    Data = double(RasterIn(:,icolV));
end
disp(' ');
disp(tit);

%% For irregular raster and change centers to edges
if iReg == 0
    switch nDim
        case 1
            % Convert to internal edges
            xE = x(1:end-1) + 0.5*diff(x);
            % Prepend and append 1st and last edges
            x = [x(1)-(xE(1)-x(1)) xE x(end)+(x(end)-xE(end))];
        case 2
            % Convert to internal edges
            xE = x(1:end-1) + 0.5*diff(x);
            yE = y(1:end-1) + 0.5*diff(y);
            % Prepend and append 1st and last edges
            x = [x(1)-(xE(1)-x(1)) xE x(end)+(x(end)-xE(end))];
            y = [y(1)-(yE(1)-y(1)) yE y(end)+(y(end)-yE(end))];
    end
end

%% Display axis limits used for plotting
switch nDim
    case 1
        if iReg == 0
            axisOut = [x(1) x(end)];
        else
            axisOut = [x(1)-0.5*xsize x(end)+0.5*xsize];
        end
    case 2
        if iReg == 0
            axisOut = [x(1) x(end) y(1) y(end)];
        else
            axisOut = [ x(1)-0.5*xsize x(end)+0.5*xsize ...
                        y(1)-0.5*ysize y(end)+0.5*ysize];
        end
end
disp('Axes limits');
tit=num2str(axisOut);
disp(tit);

%% Convert to categories from 1 to nClass
if vType == 1 
    classLabs = unique(Data); % Can return multiple NaNs
    if any(isnan(classLabs))
        classLabs(isnan(classLabs)) = []; % Delete NaNs
    end
    nClass = numel(classLabs);
    Data2 = nan(size(Data));
    for ic=1:nClass
        IndCode = Data == classLabs(ic);
        Data2(IndCode) = ic;
    end
    Data = Data2; clear Data2;
    disp(' ');
    tit = ['Found ',num2str(nClass),' categories'];
    disp(tit);
end

%% Display
newplot;
iHold = ishold;
colormap('jet');
switch nDim
    case 1
        % Row vector Data
        Data = Data';
        if iReg == 0
            % Append single node to data row
            Data = [Data 0];
            % Turn Data into a matrix
            Data = [Data; zeros(size(x))];
            % Pcolor
            hh = pcolor(x,(1:2),Data);
        else
            % Imagesc
            hh = imagesc(x,1,Data);
            if any(isnan(Data(:)))
                % Make NaNs transparent
                set(hh,'alphadata',~isnan(Data));
                if plotPars(1) == 1
                    % Handle vector graphics better
                    % Only transparent for grayscale
                    set(gcf,'Renderer','Painters');
                end
            end
        end
        % Erace y-tickmarks
        set(gca,'YTick',[]);
        % Control vertical exaggeration
        tmp = (x(end)-x(1))*YTOP1D;
        set(hh,'YData',[0 tmp]);
        axis('image'); box on;
    case 2
        % Convert to Matlab array
        Data = geoeas2matlab(Data,[nx ny]);
        if ny == 1
            Data = Data';
        end
        if nx == 1
            Data = flipud(Data);
        end
        % Flipud to account for Cartesian frame
        Data = flipdim(Data,1);
        if iReg == 0
            % Append extra column and row
            Data = [Data  zeros(size(Data,1),1)];
            Data = [Data; zeros(1,size(Data,2))];
            % Pcolor
            hh = pcolor(x,y,Data);
        else
            % Imagesc
            hh = imagesc(x,y,Data);
            axis('xy');
            if any(isnan(Data(:)))
                % Make NaNs transparent
                set(hh,'alphadata',~isnan(Data));
                if plotPars(1) == 1
                    % Handle vector graphics better
                    % Only transparent for grayscale
                    set(gcf,'Renderer','Painters');
                end
            end
        end
        axis('image'); box on;
end
if iHold == 0
    hold off;
end
% % Hide pixel boundaries for pcolor
if iEmpty == 0 && iReg == 0
    set(hh,'EdgeColor','none');
end

%% Define colormap
igray = plotPars(1);
ibar  = plotPars(2);
colormap('default');
if igray ~= 0;                      % grayscale colormap
    if vType == 1
        MapCol = gray(nClass);
    else
        MapCol = gray;
    end
    MapCol = 1 - MapCol; % invert to match GSLIB
else
    if vType == 1
        MapCol = jet(nClass);
    else
        MapCol = jet;
    end
end
colormap(MapCol);
% Set color limits for continuous case
if vType == 0
    caxis(colorLabels);
end

%% Add color bar for continuous case
if vType == 0
    switch ibar
      case 1  %%%%% vertical colorbar
        hbar = colorbar('EastOutSide');
      case 2  %%%%% horizontal colorbar
        hbar = colorbar('SouthOutSide');
      case 3  %%%%% vertical to the left
        hbar = colorbar('WestOutSide');
    end
end

%% Add colorbar for categorical case
if vType == 1 && ibar ~= 0
    % Default class labels
    classCodes = cell(1,nClass);
    for ic=1:nClass
        classCodes{ic} = [' ',num2str(classLabs(ic))];
    end
    if nargin >= NARGCOL
        tmp = varargin{NARGCOL - NARGMIN};
        if ~isempty(tmp)  %%% Read in user-specified class labels
            classCodes = tmp(:)';
        end
        % Add space in front of classes
        if ~isempty(tmp)
            for ic=1:length(tmp)
                classCodes{ic} = [' ',classCodes{ic}];
            end
        end
    end
    % Determine location of class labels in colorbar
    nPerClass = 64/nClass;
    if rem(64,nClass) == 0
        nPerClass = repmat(nPerClass,1,nClass);
    else
        nPerClass = floor(nPerClass);
        nPerClass=[repmat(nPerClass,1,nClass-1) 64-(nClass-1).*nPerClass];
    end
    ie = cumsum(nPerClass);
    is = ie - nPerClass + 1;
    iCentr = is + nPerClass./2;
    iCentr = iCentr./64;
    iCentr = 1 + (nClass - 1).*iCentr;
    % Add colorbar
    if ibar == 1 %%%% Vertical colorbar
        yTck = iCentr;
        hbar = colorbar('EastOutSide','YTick',yTck,'YTickLabel',classCodes);
        set (hbar,'YTickMode','manual','TickLength',[0 0]);
    else         %%%% Horizontal colorbar
        xTck = iCentr;
        hbar=colorbar('SouthOutSide','XTick',xTck,'XTickLabel',classCodes);
        set(hbar,'XTickMode','manual','TickLength',[0 0]);
    end
end

%% Return handles if requested
if nargout == 1
   hax = gca;
   if ibar == 0
      h = [hax; hh];
   else
      h = [hax; hh; hbar];
   end
end

%% FINISHED
disp(' ');
if nDim == 3
    tit = ['In 3D, use: axis ','"auto z"',...
           ' for automatic vertical exaggeration'];
    disp(tit);
    disp(' ');
end
disp(['Finished RASTERMAP: Version #',num2str(VERS)]);


%% Error checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for error checking
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [colorLabels,plotPars]=checkerrors(GridSpecs,varargin)
%    rastermap(GridSpecs);
%    rastermap(GridSpecs,RasterIn,icolV,vType);
%    rastermap(GridSpecs,RasterIn,icolV,vType,colorLabels);
%    rastermap(GridSpecs,RasterIn,icolV,vType,colorLabels,plotPars);
global MAXCLASS;

%% Check GridSpecs
if ~iscell(GridSpecs)
    if size(GridSpecs,2) ~= 3
       error('You need 3 columns in GridSpecs');
    end
    if size(GridSpecs,1) > 2
        error('You can only have up to 2 rows in GridSpecs');
    end
    nPixTot = prod(GridSpecs(:,1));
else
    nCells = length(GridSpecs);
    if nCells > 2
        error('This function only works up to 2D: Check GridSpecs');
    end
    nPixDir = zeros(nCells,1);
    for ii=1:nCells
        if any(diff(GridSpecs{ii})<=0)
            tit = ['Coordinates in each cell of GridSpecs',...
                   ' must be monotonically increasing'];
            disp(tit);
            error('Check GridSpecs');
        end
        nPixDir(ii) = numel(GridSpecs{ii});
    end
    nPixTot = prod(nPixDir);
end
% Return if only pixel outlines are plotted
if isempty(varargin)
    colorLabels = [0 1];
    plotPars = [0 0];
    return
end

%% Optional arguments
% Check data
RasterIn = varargin{1};
icolV = varargin{2};
[nPixTotRast,nVar] = size(RasterIn);
if nPixTotRast ~= nPixTot
    if ~iscell(GridSpecs)
        error('# of rows in RasterIn =/= prod(GridSpecs(:,1))');
    else
        error('# of rows in RasterIn =/= numel(cell2mat(GridSpecs))');
    end
end
% Check icolV
if ~isscalar(icolV)
    error('icolV must be a scalar');
end
if icolV > nVar
    error('icolV cannot be > # of columns in RasterIn');
end
if icolV < 1
    error('icolV cannot be < 1 ');
end
dataVals = RasterIn(:,icolV);
% Check vType
vType = varargin{3};
if ~ischar(vType)
    error('vType must be a string');
end
if ~(strcmp(vType,'cat') || strcmp(vType,'cont'))
    error('vType must be ''cat'' or ''cont''');
end
if strcmp(vType,'cat')
    vType = 1;
else
    vType = 0;
end
if vType == 1
    iNotNan = ~isnan(dataVals);
    if any( round(dataVals(iNotNan)) ~= dataVals(iNotNan) )
        error('For vType=''cat'', you need integer entries in DataIn(:,icolV)');
    end
    classLabs = unique(dataVals); % Can return multiple NaNs
    if any(isnan(classLabs))
        classLabs(isnan(classLabs)) = []; % Delete NaNs
    end
    nClass = numel(classLabs);
end
% Check colorLabels
if length(varargin) >= 4
    colorLabels = varargin{4};
    if ~isempty(colorLabels)
        if vType == 0 && ~isnumeric(colorLabels)
            error('For vType=''cont'', colorLabels must be a numeric array');
        end
        if vType == 0 && length(colorLabels) ~= 2 
            error('For vType=''cont'', you need 2 entries in colorLabels');
        end
        if vType == 0 && colorLabels(2) <= colorLabels(1)
            error('For vType=''cont'', colorLabels(2) cannot be <= colorLabels(1)');
        end
        if vType == 1
            if ~iscell(colorLabels)
                error('For vType=''cat'', colorLabels must be a cell array');
            end
            if length(colorLabels) ~= nClass
                error('vType=''cat'': # of cells in colorLabels ~= # of codes in DataIn'); 
            end
            for ic=1:nClass
                if ~ischar(colorLabels{ic})
                    error('vType=''cat'': Found non-string entry in colorLabels');
                end
            end
        end
    else
        if vType == 0
            colorLabels = [min(dataVals) max(dataVals)];
        else
            colorLabels = classLabs;
        end
    end
else
    if vType == 0
        colorLabels = [min(dataVals) max(dataVals)];
    else
        colorLabels = classLabs;
    end
end    
% Check plotPars
plotPars = [0 1];
if length(varargin) >= 5
    plotPars = varargin{5};
    if ~isempty(plotPars)
        if ~iscell(plotPars)
            error('plotPars must be a cell array');
        end
        if length(plotPars) ~= 2
            error('You need 2 entries in plotPars');
        end
        igray = plotPars{1};
        ibar  = plotPars{2};
        if ~(strcmp(igray,'gray') || strcmp(igray,'color'))
            error('plotPars(1), igray, must be ''color'' or ''gray');
        end
        if ~(strcmp(ibar,'nobar') || strcmp(ibar,'vbar') || strcmp(ibar,'hbar'))
            error('plotPars(2), ibar, must be ''nobar'', ''hbar'' or ''vbar');
        end
        if (vType == 1 && strcmp(ibar,'vbar')) || (vType == 1 && strcmp(ibar,'hbar'))
            if nClass > MAXCLASS
                disp('Error: Too many categories to display in colorbar');
                error('Check DataIn(:,icolV)');
            end
        end
        if strcmp(igray,'gray')
            igray = 1;
        else
            igray = 0;
        end
        if strcmp(ibar,'nobar')
            ibar = 0;
        elseif strcmp(ibar,'vbar')
            ibar = 1;
        else
            ibar = 2;
        end
        plotPars = [igray ibar];
    else
        plotPars = [0 1];
    end
end