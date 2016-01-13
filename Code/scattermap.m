function h=scattermap(DataIn,icolC,bullet,varargin)
% Display attribute values at scattered locations up to 3D
%
%% DESCRIPTION: scattermap.m
% Function to display attribute values at scattered sample locations.
% Attribute values are displayed with color-coded, filled cicles.
% NOTE: NaNs are NOT allowed in the coordinate columns,
%       whereas they are allowed in the attribute column
% NOTE: Sample locations with missing attribute values are displayed 
%       as circles with crosses
% NOTE: To change the color of bullet outline use:
%           set(h(2),'MarkerEdge','r')  % here, red
%       To change the line width of bullet outline use:
%           set(h(2),'LineWidth',2)     % here, 2
% NOTE: Exporting to an .eps file: print -deps (or -depsc) file.eps
% 
%% SYNTAX:
%    scattermap(DataIn,icolC,bullet);
%    scattermap(DataIn,icolC,bullet,icolV,vType);
%    scattermap(DataIn,icolC,bullet,icolV,vType,colorLabels);
%    scattermap(DataIn,icolC,bullet,icolV,vType,colorLabels,plotPars);
%    h = scattermap(...);
%
%% INPUTS:
%    DataIn     = (nDat x nCol) input arrray with sample coordinates 
%                 and data values
%    icolC      = (1 x nDim) array with column numbers (in DataIn) 
%                 with coordinates
%                 NOTE: No NaNs are allowed in coordinate columns
%                 NOTE: Dimensionality is determined by this array
%    bullet     = scalar to control size of color-filled circles
%                 drawn at sample locations
%                 NOTE: Use a value between 10 and 30
%    icolV      = Optional scalar with column # (in DataIn) with 
%                 attribute values to display
%    iType      = string with flag for continuous ('cont')
%                 or categorical ('cat') attribute raster
%                 NOTE: This affects the color map used for each case
%   colorLabels = Optional array with color limits OR category labels:
%                 For the continuous case (vType='cont'), 
%                       this is a (1 x 2) array [cmin cmax] with:  
%                       cmin = min color value for plotting
%                       cmax = max color value for plotting
%                 For the categorical case (vType='cat'),
%                       this is a (1 x nClass) cell array with strings for
%                       class labels, e.g., {'A','B','C'}
%                 NOTE: This could be an EMPTY array, in which case default
%                       color limits (min/max attribute values) and
%                       default class labels ('1','2',...) are used
%                 NOTE: For the categorical case, this array is used only 
%                       if a colorbar is requested (see below)
%    plotPars   = Optional (1 x 2) cell array {igray ibar} with flags
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
%  + handles to axes, patches and colorbar (if requested)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%    scattermap(DataIn,icolC,bullet);
%    scattermap(DataIn,icolC,bullet,icolV,vType);
%    scattermap(DataIn,icolC,bullet,icolV,vType,colorLabels);
%    scattermap(DataIn,icolC,bullet,icolV,vType,colorLabels,plotPars);
%    h = scattermap(...);

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
NARGMIN  = 3;
NARGDAT  = NARGMIN  + 2;
NARGCOL  = NARGDAT  + 1;
NARGPLOT = NARGCOL  + 1;
global MAXCLASS;
MAXCLASS = 10;
MARKCOLOR = 'k';      % Color of lines around circles
MARKLINE = 1;         % Width of lines around circles
MARKLINE3D = 0.1;     % Width of lines around circles in 3D
MARKMISS = 'k+';
VERS  = 1;

%% Error checking
if ~(   nargin == NARGMIN || nargin == NARGDAT || nargin == NARGCOL ...
     || nargin == NARGPLOT) 
   tit=['You need ',num2str(NARGMIN),', ',num2str(NARGDAT),', ',...
        num2str(NARGCOL),' or ',num2str(NARGPLOT),...
        ' input arguments for this function'];
   disp(tit);
   error('Check # of input arguments');
end
[colorLabels,plotPars]=checkerrors(DataIn,icolC,bullet,varargin{:});
% Note: plotPars = [igray ibar]
%       igray = 0, 1    -> color, gray scale 
%       ibar  = 0, 1, 2 -> nobar, vbar, hbar

%% Extract columns of coordinates
nDim = numel(icolC);
CoordsDat = single(DataIn(:,icolC));
nLocs = size(CoordsDat,1);
titLoc = ['Number of locations in input array = ', num2str(nLocs)];
disp(titLoc);

%% Handle case of bullet outlines
if nargin == NARGMIN
    switch nDim
        case 1
            y1 = zeros(nLocs,1);
            hscatt = scatter(CoordsDat,y1,bullet,...
                     'MarkerEdgeColor',MARKCOLOR,'LineWidth',MARKLINE);
            axis('image');
            set(gca,'YtickLabel',[],'YTick',[]);
            %xlabel('x-coordinates');
        case 2
            hscatt = scatter(CoordsDat(:,1),CoordsDat(:,2),bullet,...
                     'MarkerEdgeColor',MARKCOLOR,'LineWidth',MARKLINE);
            axis('image');
            %xlabel('x-coordinates');
            %ylabel('y-coordinates');
        case 3
            
            % In 3D, do not draw the outline of circles?
            hscatt = scatter3(CoordsDat(:,1),CoordsDat(:,2),CoordsDat(:,3),...
                     bullet,'MarkerEdgeColor',MARKCOLOR,...
                     'LineWidth',MARKLINE3D);
            grid off;
            %axis('image');  
            % Do not force axes in 3D to allow vertical exaggeration
            xlabel('x-coordinates');
            ylabel('y-coordinates');
            zlabel('z-coordinates');
    end
    box on;
    if nargout == 1
        h = [gca; hscatt];
    end
    disp(' ');
    disp(['Finished SCATTERMAP: Version #',num2str(VERS)]);
    return
end

%% Figure out # of data values in file and # finally plotted
icolV = varargin{1};
iType = varargin{2};
if strcmp(iType,'cat')
    iCat = 1;
else
    iCat = 0;
end
dataVals = DataIn(:,icolV);
in = find(isnan(dataVals));
nMiss = numel(in);
nDat = size(dataVals,1) - nMiss;
titDat = ['Number of attribute values plotted = ', num2str(nDat)];
disp(' ');
disp(titDat);

%% Convert to categories from 1 to nClass
if iCat == 1
    classLabs = unique(dataVals); % Can return multiple NaNs
    if any(isnan(classLabs))
        classLabs(isnan(classLabs)) = []; % Delete NaNs
    end
    nClass = numel(classLabs);
    dataVals2 = nan(size(dataVals));
    for ic=1:nClass
        IndCode = dataVals == classLabs(ic);
        dataVals2(IndCode) = ic;
    end
    dataVals = dataVals2; clear dataVals2;
    disp(' ');
    tit = ['Found ',num2str(nClass),' categories'];
    disp(tit);
end

%% Do the plotting
newplot;
iHold = ishold;
colormap('jet');
switch nDim
    case 1
        y1 = zeros(nLocs,1);
        hscatt = scatter(CoordsDat,y1,bullet,dataVals,...
                'filled','MarkerEdgeColor',MARKCOLOR,'LineWidth',MARKLINE);
        axis('image');
        set(gca,'YtickLabel',[],'YTick',[]);
        if nMiss > 0
            hold on;
            plot(CoordsDat(in,1),y1,MARKMISS);
        end
        %xlabel('x-coordinates');
    case 2
        hscatt = scatter(CoordsDat(:,1),CoordsDat(:,2),bullet,dataVals,...
                'filled','MarkerEdgeColor',MARKCOLOR,'LineWidth',MARKLINE);
        if nMiss > 0
            hold on;
            plot(CoordsDat(in,1),CoordsDat(in,2),MARKMISS);
        end
        axis('image');
        %xlabel('x-coordinates');
        %ylabel('y-coordinates');
    case 3
        % In 3D, do not draw the outline of circles?
        hscatt = scatter3(CoordsDat(:,1),CoordsDat(:,2),CoordsDat(:,3),...
                 bullet,dataVals,'filled','MarkerEdgeColor',MARKCOLOR,...
                 'LineWidth',MARKLINE3D);
        if nMiss > 0
            hold on;
            plot3(CoordsDat(in,1),CoordsDat(in,2),CoordsDat(in,3),...
                  MARKMISS);
        end
        grid off;
        %axis('image');  
        % Do not force axes in 3D to allow vertical exaggeration
        xlabel('x-coordinates');
        ylabel('y-coordinates');
        zlabel('z-coordinates');
end
if iHold == 0
    hold off;
end
box on;

%% Define colormap
igray = plotPars(1);
ibar  = plotPars(2);
colormap('default');
if igray ~= 0;                      % grayscale colormap
    if iCat == 1
        MapCol = gray(nClass);
    else
        MapCol = gray;
    end
    MapCol = 1 - MapCol; % invert to match GSLIB
else
    if iCat == 1
        MapCol = jet(nClass);
    else
        MapCol = jet;
    end
end
colormap(MapCol);
% Set color limits for continuous case
if iCat == 0
    caxis(colorLabels);
end

%% Add color bar for continuous case
if iCat == 0
    switch ibar
      case 1  %%%%% vertical colorbar
        hbar = colorbar('EastOutSide');
      case 2  %%%%% horizontal colorbar
        hbar = colorbar('SouthOutSide');
    end
end

%% Add colorbar for categorical case
if iCat == 1 && ibar ~= 0
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

%% Return handle if requested
if nargout == 1
   hax = gca;
   if ibar == 0
      h = [hax; hscatt];
   else
      h = [hax; hscatt; hbar];
   end
end

%% FINISHED
disp(' ');
disp(['Finished SCATTERMAP: Version #',num2str(VERS)]);

%% Function for error checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to perform error checking
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [colorLabels,plotPars]=checkerrors(DataIn,icolC,bullet,varargin)
%    scattermap(DataIn,icolC,bullet);
%    scattermap(DataIn,icolC,bullet,icolV,vType);
%    scattermap(DataIn,icolC,bullet,icolV,vType,colorLabels);
%    scattermap(DataIn,icolC,bullet,icolV,vType,colorLabels,plotPars);
global MAXCLASS;
nCols = size(DataIn,2);
% Check icolC
if numel(unique(icolC)) ~= numel(icolC)
    error('You need distinct entries in icolC');
end
if any(icolC > nCols )
    error('Found entry in icolC > # of columns in DataIn');
end
if any(icolC < 1)
    error('Found entry in icolC < 1');
end
CoordsDat = single(DataIn(:,icolC));
if any(isnan(CoordsDat))
    error('Found NaN in DataIn(:,icolC)');
end
nDim = numel(icolC);
if nDim > 3
    error('This function works only up to 3D');
end
% Check for duplicate points
CoordsDat = sortrows(CoordsDat);
if nDim == 1
    indS = diff(CoordsDat)==0;
else
    indS = all(diff(CoordsDat)'==0);
end
if any(indS)
    error('Found duplicate points in the data set');
end
clear CoordsDat indS;
% Check bullet
if numel(bullet) ~= 1
    error('bullet must be a scalar');
end
% Return if only bullets are plotted
if isempty(varargin)
    colorLabels = [0 1];
    plotPars = [0 0];
    return
end
%% Optional arguments
% Check icolV
icolV = varargin{1};
if numel(icolV) ~= 1
    error('icolV must be a scalar');
end
if any(icolC == icolV)
    error('Some entries in icolC are the same with that in icolV');
end
if icolV > nCols
    error('icolV > # of columns in DataIn');
end
if icolV < 1
    error('icolV < 1');
end
dataVals = DataIn(:,icolV);
% Check vType
vType = varargin{2};
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
if length(varargin) >= 3
    colorLabels = varargin{3};
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
if length(varargin) == 4
    plotPars = varargin{4};
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