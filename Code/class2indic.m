function Out = class2indic(DataIn,icolV,icolC)
% Indicator transform of class labels
%
%% DESCRIPTION: class2indic.m
% Function to compute the indicator transform of class labels
% NOTE: NaNs are allowed in any column of DataIn
% NOTE: Any NaN in the column of class labels is copied to the output
%
%% SYNTAX: 
%    Out = class2indic(DataIn,icolV);
%    Out = class2indic(DataIn,icolV,icolC);
%
%% INPUTS: 
%    DataIn = (nVals x nCols) input data array
%    icolV  = scalar with column # (in DataIn) containing class labels
%    icolC  = Optional (1 x nDim) array with column #s (in DataIn) with
%             coordinates of sample locations
%             NOTE: If icolC is provided, the sample coordinates are
%                   pre-pended to the output array of indicators
%
%% OUTPUTS:
%    Out    = (nVals x nClass) output array with indicators for each class
%             NOTE: In this case Out is single precision
%                                 OR
%             (nVals x nClass+nDim) output array with:
%             col#1:nDim: sample coordinates
%             col#nDim+1:nDim+nClass: class indicators
%             NOTE: In this case Out is double precision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%    Out = class2indic(DataIn,icolV);
%    Out = class2indic(DataIn,icolV,icolC);

%% CREDITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                               May 2005                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixed parameters
NARGMIN   = 2;
NARGMAX   = NARGMIN + 1;

%% Error checking
if ~( nargin == NARGMIN || nargin == NARGMAX )
    tit = ['You need ',num2str(NARGMIN),' or ',num2str(NARGMAX),...
           ' input arguments for this function'];
    disp(' ');
    disp(tit);
    error('Check # of input arguments');
end
if ~isscalar(icolV)
    error('icolV must be a scalar');
end
if icolV < 1
    error('icolV must be >= 1');
end
[nValsIn,nCols] = size(DataIn);
if icolV > nCols
    error('icolV must be <= # of columns in DataIn');
end
if nargin == NARGMAX
    if isempty(icolC)
        error('You need some entries in icolC');
    end
    if numel(unique(icolC)) ~= numel(icolC)
        error('Entries in icolC must be unique');
    end
    if any(icolC<1)
        error('Found entry in icolC < 1');
    end
    if any(icolC>nCols)
        error('Found entry in icolC > # of colums in DataIn');
    end
    if any(icolC==icolV)
        error('Entries in icolC must be different than icolV');
    end
end

%% Get some input parameters
classValsIn = DataIn(:,icolV);
iNotNan = ~isnan(classValsIn);
nVals = sum(iNotNan);
if sum(iNotNan) == nValsIn
    iMiss = 0;
    classVals = classValsIn;
else
    iMiss = 1;
    classVals = classValsIn(iNotNan);
end
classLabs = unique(classVals);
nClass = numel(classLabs);

%% Declare some arrays
if nargin == NARGMIN 
    Out = zeros(nValsIn,nClass,'single');
else
    Out = zeros(nValsIn,nClass);
end
prop = zeros(1,nClass);

%% Loop over # of classes and compute indicators
for ic=1:nClass
    classIndic = classVals == classLabs(ic);
    prop(ic) = sum(classIndic)./nVals;
    Out(iNotNan,ic) = classIndic;
end

%% Handle missing values
if iMiss == 1
    Out(~iNotNan,:) = NaN;
end

%% Account for coordinates columns
if nargin == NARGMAX
    Out = [DataIn(:,icolC) Out];
end

%% Finished
disp(' ');
disp('Finished CLASS2INDIC: Version #1');
tit = ['Found ',num2str(nClass),' categories, with proportions:'];
disp(' ');
disp(tit);
tit = num2str(prop);
disp(tit);