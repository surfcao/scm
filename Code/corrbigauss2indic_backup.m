function Out = corrbigauss2indic(ModPars,lagDist,azims,dips,globProps,MeasSpecs)
% Indicator covariograms/semivariograms from correlogram of Gaussian RF
%
%% DESCRIPTION: corrbigauss2indic.m
% Function to compute model (theoretical) measures of spatial structure 
% (joint probabilities/covariograms/semivariograms/transition probabilities) 
% at a set of specified directions and lags (up to 3D), given an input 
% correlogram model of a Gaussian random function, i.e., a set of
% correlation coefficients for a set of h-specific bivariate Gaussian PDFs
%
%% SYNTAX:
%    Out=corrbigauss2indic(ModPars,lagDist,azims,dips,globProps,MeasSpecs);
%
%% INPUTS:
%    ModPars  = (nStruct x 4,6,9) array with correlogram parameters for
%               a standard Gaussian random function
%               Examples with 2 nested structures:
%               1D: [mtype(1) msill(1) rng1(1) extrapar(1);
%                    mtype(2) msill(2) rng1(2) extrapar(2)]
%               2D: [mtype(1) msill(1) ang1(1) rng1(1) rng2(1) extrapar(1);
%                    mtype(2) msill(2) ang1(2) rng1(2) rng2(2) extrapar(2)]
%               3D: [mtype msill ang1 ang2 ang3 rng1 rng2 rng3 extrapar;
%                    mtype msill ang1 ang2 ang3 rng1 rng2 rng3 extrapar]
%               mtype = functional form (e.g., spherical, gaussian, etc).
%               msill = variance of nested structure
%               ang1, ang2, ang3 = anisotropy rotation angles (see Notes)
%               rng1, rng2, rng3 = ranges of anisotropy ellipsoid
%               extrapar = extra parameter used for certain models,
%                          e.g., dampened cosine model.
%               See function dist2msemivar for model specifications
%               NOTE: Dimensionality is determined by this array
%    lagDist  = (nDir x 1) cell array, with each cell containing the 
%               (nLag(idir)) lag distances for each direction.
%    azims    = (nDir x 1) array with lag azimuth angles for each direction
%                   in 1D, azims must be an EMPTY array
%                   in 2D & 3D, azims must contain as many entries as
%                         # of cells in array lagDist
%               NOTE: These are NOT the same angles as those specifying 
%                     the orientation of the anisotropy ellipse or 
%                     ellipsoid in array ModPars
%    dips     = (nDir x 1) array with lag dip angles for each direction
%                   in 1D & 2D, dips must be an EMPTY array
%                   in 3D, dips must contain as many entries as # of cells
%                         in array lagDist
%   globProps = (1 x nClass) array with global (stationary) class
%               proportions
%               NOTE: The actual class labels do not matter
%   MeasSpecs = (nMeas x 3) array with tail & head class (order i and j
%               in array globProps) & indicator structural measure type 
%               to compute for each direction:
%               e.g.: [ iclassT(1) iclassH(1) iMeas(1);
%                       iclassT(2) iclassH(2) iMeas(2) ]
%               iclassT(im) = tail category
%               iclassH(im) = head category
%               iMeas(im)  = indicator structural measure type:
%                          = 0 -> non-centered indicator covariogram
%                          = 1 -> indicator semivariogram
%                          = 2 -> indicator cross-semivariogram
%                          = 3 -> indicator covariogram
%                          = 4 -> indicator correlogram
%                          = 5 -> transition probability diagram
%                e.g.: [ 1 1 1 ] -> indicator semivariogram of 
%                                   class #1 and #2 (corresponding to 
%                                   first 2 entries in globProps)
%                      [ 2 2 3 ] -> indicator covariogram of class #2
%                                   (2nd enty in globProps)
%                      [ 1 2 2 ] -> indicator cross-semivariogram between
%                                   class #1 and class #2
%                      [ 1 1 3 ; 2 2 3; 1 2 3] -> 
%                              indicator covariogram of class #1
%                              indicator covariogram of class #2,
%                              indicator cross-covariogram of classes #1,#2
%
%% OUTPUTS:
%    Out      = (6 x 1) output structure array:
%               Out.title  = text description of array contents
%               Out.model  = ModPars array for Gaussian field copied here
%               Out.azim   = [azims(idir)] azimuth angle per direction
%                            NOTE: This array is EMPTY in 1D cases
%               Out.dip    = [dips(idir)] dip angle per direction
%                            NOTE: This array is empty in 1D & 2D cases
%               Out.meas   = array MeasSpecs copied from input
%               Out.data   = (nDir x 1) cell array, with each cell
%                            containing a (nLag(idir) x 5 x nMeas) array:
%                            col#1:   lag distances (from array lagDist)
%                            col#2:   indicator structural measure value
%                            col#3:4: tail and head classes
%                            col#5:   imeas = structural measure type
%                            Same column sequence is repeated for different
%                            combinations of tail and head classes and
%                            measure types, thus building the 3rd dimension
%                            of this array
%
%% NOTES:        
%               All angles are measured with the observer standing
%               at the end of the axis and looking towards its origin
%               positive angle = clockwise rotation
%               ang1: rotation (around z-axis) of the xy plane
%                     positive = clockwise rotation from NS azimuth
%               ang2: dip rotation (around x-axis) of the yz plane
%                     positive = dipping down from y-axis
%               ang3: plunge rotation (around y-axis) of the xz plane
%                     positive = clockwise down from x-axis
%
%% CUSTOM FUNCTIONS CALLED:
%    checkmodelstruct
%    rotationmatrix
%    dist2mcov
%    bigausscdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% SYNTAX:
%    Out=corrbigauss2indic(ModPars,lagDist,azims,dips,globProps,MeasSpecs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                            September 2007                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixed parameters
VERS = 1;
NARGMIN = 6;

%% Error checking
if ~( nargin == NARGMIN )
    tit = ['You need ',num2str(NARGMIN),...
           ' input arguments for this function'];
    disp(' ');
    disp(tit);
    error('Check # of input arguments');
end
checkerrors(ModPars,lagDist,azims,dips,globProps,MeasSpecs);

%% Get dimensionality
switch size(ModPars,2)
  case 4
      nDim = 1;
  case 6
      nDim = 2;
  case 9
      nDim = 3;
end

%% Check for consistency in entries of ModPars
checkmodelstruct(ModPars,nDim);

%% Compute rotation matrix for all structures
ModPars(:,2) = ModPars(:,2)./sum(ModPars(:,2));
switch nDim
  case 1
      [RotMat,iAnis] = rotationmatrix(ModPars(:,3));
  case 2
      [RotMat,iAnis] = rotationmatrix(ModPars(:,[3 4 5]));
  case 3
      [RotMat,iAnis] = rotationmatrix(ModPars(:,[3 4 5 6 7 8]));
end

%% Get some parameters
nStruct  = size(ModPars,1);
nDir = length(lagDist);
nLag = zeros(nDir,1);
for id=1:nDir
    nLag(id) = numel(lagDist{id});
end
nMeas = size(MeasSpecs,1);
nClass = numel(globProps);
switch nDim
    case 1
        azimsOut = [];
        dipsOut  = [];
    case 2
        azimsOut = azims;
        azims = azims(:).*(pi/180);
        dips = zeros(nDir,1);
        dipsOut = [];
    case 3
        azimsOut = azims;
        azims = azims(:).*(pi/180);
        dipsOut = dips;
        dips    = dips(:).*(pi/180);
end

%% Initialize output structure array
Out.title  = 'Gaussian-based indicator measures of spatial structure';
Out.model  = ModPars;
Out.azim   = azimsOut(:)';
Out.dip    = dipsOut(:)';
Out.meas   = MeasSpecs;
Out.data   = cell(nDir,1);
for id=1:nDir
    Out.data{id} = nan(numel(lagDist{id}),5,nMeas);
end

%% Compute model correlation values along nDir directions
corrGauss = cell(nDir,1);
switch nDim
 case 1   %%%%%%%%%%%%%%%%% 1D case
   lagDistCurr = lagDist{1};
   lagDistCurr = lagDistCurr(:);
   covDir = dist2mcov(abs(lagDistCurr),[nLag(1) 1],ModPars,nStruct);
   corrGauss{1} = covDir;

 case 2   %%%%%%%%%%%%%%%%% 2D case
   %ang1=ModPars(:,3); rng1=ModPars(:,4); rng2=ModPars(:,5);
   %extrapar=ModPars(:,6); 
   for id = 1:nDir
       lagDistCurr = lagDist{id};
       lagDistCurr = lagDistCurr(:);
       xcomp = sin(azims(id))*cos(dips(id));  % x-comp of input unit lag
       ycomp = cos(azims(id))*cos(dips(id));  % y-comp of input unit lag
       dx = lagDistCurr*xcomp;
       dy = lagDistCurr*ycomp;
       covDir = zeros(nLag(id),1);
       for ist=1:nStruct
           % Account for anisotropy
           if iAnis(ist) == 1
                % Rotate and scale vector components
                dAnis = [dx dy]*RotMat(:,:,ist);
           else
                dAnis = [dx dy];
           end
           % Compute anisotropic distance
           dAnis = sqrt(sum(dAnis.^2,2));
           % Compute contribution of current nested structure
           p = ModPars(ist,[1 2 4 end]);
           covDir  = covDir + dist2mcov(dAnis,[nLag(id) 1],p,1);
       end
       corrGauss{id} = covDir;
   end

   case 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D case
   %ang1=ModPars(:,3); ang2=ModPars(:,4); ang3=ModPars(:,5);
   %rng1=ModPars(:,6); rng2=ModPars(:,7); rng3=ModPars(:,8);
   %extrapar=ModPars(:,9);
   for id = 1:nDir
       lagDistCurr = lagDist{id};
       lagDistCurr = lagDistCurr(:);
       xcomp = sin(azims(id))*cos(dips(id));  % x-comp of input unit lag
       ycomp = cos(azims(id))*cos(dips(id));  % y-comp of input unit lag 
       zcomp = sin(dips(id));                 % z-comp of input unit lag
       dx = lagDistCurr*xcomp;
       dy = lagDistCurr*ycomp;
       dz = lagDistCurr*zcomp;
       covDir = zeros(nLags(id),1);
       for ist=1:nStruct
           % Account for anisotropy
           if iAnis(ist) == 1
                % Rotate and scale vector components
                dAnis = [dx dy dz]*RotMat(:,:,ist);
           else
                dAnis = [dx dy dz];
           end
           % Compute anisotropic distance
           dAnis = sqrt(sum(dAnis.^2,2));
           % Compute contribution of current nested structure
           p = ModPars(ist,[1 2 6 end]);
           covDir  = covDir + dist2mcov(dAnis,[nLag(id) 1],p,1);
       end
       corrGauss{id} = covDir;
   end
end        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch nDim

%% Compute Gaussian thresholds
globProps = double(globProps(:));
globProps(globProps > (1-eps)) = 1-eps;
globProps(globProps < eps) = eps;
globProps = globProps./sum(globProps);
thresh = cumsum(globProps(1:nClass-1));
thresh = -sqrt(2).*erfcinv(2.*thresh);
thresh = [-inf; thresh; inf];
varInd = globProps.*(1-globProps);
stdInd = sqrt(varInd);

%% Loop over directions and compute structural measures
for id=1:nDir
    for il=1:nLag(id)
        % Covariance matrix (Gaussian) for current lag
        %S = [1 covGauss{id}(il); covGauss{id}(il) 1];
        % Current lag correlation coefficient
        corrCurr = corrGauss{id}(il);
        for im=1:nMeas
            % Determine lower and upper integration limits
            iT = MeasSpecs(im,1);
            iH = MeasSpecs(im,2);
            yLowT = thresh(iT);
            yLowH = thresh(iH);
            yLow = [yLowT yLowH];
            yHighT = thresh(iT+1);
            yHighH = thresh(iH+1);
            yHigh = [yHighT yHighH];
            % Handle 0 distance
            if lagDist{id}(il) < eps
                if iT == iH
                    jointProb = globProps(iT);
                else
                    jointProb = 0;
                end
            else
                % Integrate to compute joint probabilities
                jointProb = bigausscdf(yLow,yHigh,corrCurr);
            end
            % Compute requested measure
            switch MeasSpecs(im,3)
            case 0 % Non-centered covariance
                % Nothing to change
                outMeas = jointProb;
            case 1 % Semivariogram
                outMeas = globProps(iT) - jointProb;
            case 2 % Cross-semivariogram
                outMeas = -jointProb;
            case 3 % Centered covariance
                pT = globProps(iT);
                pH = globProps(iH);
                outMeas = jointProb - pT.*pH;
            case 4 % Correlogram
                pT = globProps(iT);
                pH = globProps(iH);
                outMeas = jointProb - pT.*pH;
                outMeas = outMeas./(stdInd(iT).*stdInd(iH));
            case 5 % Transition probability
                pT = globProps(iT);
                outMeas = jointProb./pT;
            end
            % Populate output cell array
            Out.data{id}(il,1,im) = lagDist{id}(il);
            Out.data{id}(il,2,im) = outMeas;
            Out.data{id}(il,3,im) = MeasSpecs(im,1);
            Out.data{id}(il,4,im) = MeasSpecs(im,2);
            Out.data{id}(il,5,im) = MeasSpecs(im,3);
        end
    end
end

%% FINISHED
disp(' ');
disp(['Finished CORRBIGAUSS2INDIC: Version #',num2str(VERS)]);

%% Error checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function for error checking
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkerrors(ModPars,lagDist,azims,dips,globProps,MeasSpecs)

if ~iscell(lagDist)
    error('lagDist must be a cell array');
end
if size(lagDist,2) ~= 1
    error('lagDist must be a column cell array');
end
nDir = length(lagDist);
nClass = numel(globProps);
if any(globProps<0) || any(globProps>1)
    error('Entries of globProps must be in [0 1]');
end
if sum(globProps) ~= 1
    disp(' ');
    disp('WARNING: Sum of class proportions ~= 1; will standardize');
    disp(' ');
end
if size(MeasSpecs) ~= 3
    error('You need 3 columns in MeasSpecs');
end
nMeas = size(MeasSpecs,1);
for im=1:nMeas
    if MeasSpecs(im,1) < 1
        error('Found entry in MeasSpecs(:,1) < 1');
    end
    if MeasSpecs(im,1) > nClass
        error('Found entry in MeasSpecs(:,1) > numel(globProps)');
    end
    if MeasSpecs(im,2) < 1
        error('Found entry in MeasSpecs(:,2) < 1');
    end
    if MeasSpecs(im,2) > nClass
        error('Found entry in MeasSpecs(:,2) > numel(globProps)');
    end
    iMeas = MeasSpecs(im,3);
    if ~( iMeas == 0 || iMeas == 1 || iMeas == 2 || iMeas == 3 || ...
          iMeas == 4 || iMeas == 5)
        error('Entries in MeasSpecs(:,3) must be 0, 1, 2, 3, 4 or 5');
    end
    if MeasSpecs(im,3) == 1
        if MeasSpecs(im,1) ~= MeasSpecs(im,2)
            error('Semivariogram between different tail and head classes');
        end
    end
    if MeasSpecs(im,3) == 2
        if MeasSpecs(im,1) == MeasSpecs(im,2)
            error('Cross-semivariogram between same tail and head classes');
        end
    end
end
switch size(ModPars,2)
  case 4 % 1D
      if nDir > 1
         error('In 1D, there is only one direction');
      end
      if length(lagDist) ~= 1
          error('In 1D, lagDist must contain only 1 cell');
      end
      if ~isempty(azims)
          error('In 1D, azims must be an empty array');
      end
      if numel(dips) ~= 0
          error('In 1D, dips must be an empty array');
      end
  case 6 % 2D
      if numel(dips) ~= 0
          error('In 2D, dips must be an empty array');
      end
      if length(azims) ~= nDir
          error('# of entries in azims =/= # of cells in lagDist');
      end
  case 9 % 3D
      if length(azims) ~= nDir
          error('# of rows in azims =/= # of cells in lagDist');
      end
      if length(dips) ~= nDir
          error('# of rows in dips =/= # of cells in lagDist');
      end
  otherwise
      error('Invalid specification of ModPars'); 
end
if sum(ModPars(:,2)) ~= 1
    disp('WARGNING: Sill of inpute model ~= 1: Will standardize...');
end