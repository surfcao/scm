function Out=raster2structdir(RastIn,GridSpecs,DirOffs,nLag,MeasSpecs,iOmni)
% Compute directional measures of spatial structure from attribute rasters
%
%% DESCRIPTION: raster2structdir.m
% Function to calculate alternative directional measures of spatial
% structure (auto- or cross- semivariograms/covariograms/correlograms) 
% from multiple rasters (up to 3D). The raster cell size must be regular 
% (fixed), but it need not be square. Optionally, the directional measures
% can be averaged into a pseudo omni-directional measure, which is 
% appended at the end of the other directional measures.
% NOTE: This function is a Matlab translation of the GSLIB function gam.
% NOTE: This function can accommodate indicator (0/1) variables,
%       but no error checking is performed to ensure that the columns
%       of RastIn contain indeed such 0/1 entries.
% NOTE: For the case of a cross-semivariogram, the output tail and
%       head mean and variance are those of v1(s+h)-v1(s) and 
%       v2(s+h)-v2(s), respectively. This is NOT the same as the
%       corresponding output of the GSLIB function gam, which outputs
%       the mean of 0.5*(v1(s+h)+v1(s)) and the mean of 0.5*(v2(s+h)+v2(s))
% NOTE: The following convention is used for tail and head variables:
%       A lag of h=+1 entails that vrT(i=1) is paired with vrH(j=2),
%       hence tail indices are 1:nx-1 and head indices are 1+1:nx
%       This is NOT the convention of GSLIB function gam, whose output
%       is reversed in terms of direction. In other words, the output
%       for a lag h=+1 in the x direction corresponds to a lag -1 0 0
%       in GSLIB. The convention used in this function is CONSISTENT
%       with function samplestruct and rasters2structgrid, whereas this
%       is NOT the case for GSLIB's functions gam and gamv.
% NOTE: In 3D, a positive izoff value is the z-component of a lag vector
%       pointing DOWNWARDS; this is consistent with the definition of
%       the dip angle which increases (positively) downwards.
% NOTE: For iOmni=1, no error cheking is performed to ensure that the
%       distances along the corresponding directions are similar...
%
%% SYNTAX:
%   Out = raster2structdir(RastIn,GridSpecs,DirOffs,nLag,MeasSpecs,iOmni);
%
%% INPUTS:
%   RastIn    = (nPix x nCol) array with input attribute rasters in 
%               GeoEAS format; each column contains a different attribute
%               raster and all rasters have the same grid specs.
%               NOTE: RastIn can hold 1D, 2D, or 3D, attribute rasters 
%                     with missing values flagged as NaNs.
%   Gridspecs = (nDim x 3) array with grid specifications for RastIn
%               e.g., in 2D: [nx xmin xsize ;  ny ymin ysize]
%               nx = number of nodes in x-direction
%               xmin = Cartesian x-coordinate of grid origin
%                      (USE AN ARBITRARY VALUE)
%               xsize = size of grid cell in x-direction
%               ny = number of nodes in x-direction
%               ymin = Cartesian y-coordinate of grid origin
%                      (USE AN ARBITRARY VALUE)
%               ysize = size of grid cell in y-direction
%   DirOffs   = (nDir x nDim) array [ixoff iyoff izoff] with # of pixels
%               along each direction specifying the magnitude and 
%               direction of each lag vector:
%               ixoff = x-component of lag vector
%               iyoff = y-component of lag vector
%               1D example: [  1 ] -> tail -> head (WE)
%                           [ -1 ] -> head <- tail (EW)
%                           [ 1 ; -1] -> both of the above
%                     NOTE: DirOffs must be a column vector in 1D
%               2D example: [1 0] -> WE-direction
%                           [0 1] -> SN-direction
%                           [1 1] -> N45E-direction (for square pixels)
%                           [1 0 ; 0 1] -> WE and SN directions
%               3D example: [0 0  1] -> vertical downward direction
%                           [0 0 -1] -> vertical upward direction
%   nLag      = (1 x nDir) vector with number of lags per direction,
%               i.e., # of mulitples of lag vector specified in DirOffs
%   MeasSpecs = (nMeas x 3) array with column #s for tail & head 
%               attribute rasters & structural measure specs to compute
%               for each direction:
%               e.g.: [ icolT(1) icolH(1) iMeas(1);
%                       icolT(2) icolH(2) iMeas(2) ]
%               icolT(im) = tail variable raster
%               icolH(im) = head variable raster
%               iMeas(im)  = structural measure type:
%                          = 0 -> non-centered covariogram
%                          = 1 -> semivariogram
%                          = 2 -> cross-semivariogram
%                          = 3 -> covariogram
%                          = 4 -> correlogram
%                          = 5 -> transition probability diagram
%                          = 6 -> general relative semivariogram
%                          = 7 -> pairwise relative semivariogram
%                          = 8 -> semimadogram
%                e.g.: [ 1 1 1 ] -> semivariogram of raster in col#1
%                      [ 2 2 3 ] -> covariogram of raster in col#2
%                      [ 1 2 2 ] -> cross-semivariogram between
%                                   raster in col#1 and raster in col#2
%                      [ 1 1 3 ; 2 2 3; 1 2 3] -> 
%                                covariogram of raster in col#1,
%                                covariogram of raster in col#2,
%                                cross-covariogram of rasters in cols#1,2
%   iOmni     = scalar with flag to average directional measures into
%               a single pseudo omni-directional measure (=1) or not;
%               this measure is appended below the other directional ones.
%
%% OUTPUTS:
%   Out       = Output structure array with 5 fields:
%               Out.title  = text description of array contents
%               Out.datafl = name of input array with attribute rasters
%                            from which structural measures were computed
%               Out.dir    = (nDir x 1,3,6) array with direction specs:
%                             1D:  90 -> right shift
%                                 -90 -> left shift
%                             2D: [azim azimTol horBand],
%                                 azim = azimuth angle (degrees from NS)
%                                 azimTol = azimuth tolerange (degrees)
%                                 horBand = horizontal bandwidth
%                                 NOTE: azimTol and horBand are set to 
%                                        a small number 0.001
%                             3D: [azim azimTol horBand dip dipTol verBahd]
%                                 dip = dip angle (degrees from NS)
%                                 dipTol = dip tolerange (degrees)
%                                 verBand = vertical bandwidth
%                                 NOTE: dipTol and verBand are set to 
%                                       a small number 0.001
%               Out.meas   = array MeasSpecs copied from input
%               Out.struct = (nDir x 1) cell array, with each cell
%                            containing a (nLag(idir) x 10 x nMeas) array:
%                            col#1:  lag distances
%                            col#2:  corresponding sample structure values
%                            col#3:  number of pairs per lag
%                            col#4:  column # for tail variable raster
%                            col#5:  column # for head variable raster
%                            col#6:  type of structural measure
%                            col#7:  mean of tail values
%                            col#8:  mean of head values
%                            col#9:  variance of tail values 
%                            col#10: variance of head values
%                            Same column sequence is repeated for different
%                            combinations of tail and head variables and
%                            measure types, thus building the 3rd dimension
%                            of this array
%                  NOTE: if iOmni = 1, Out.struct has 1 additional cell
%                        corresponding to the pseudo omni-directional
%                        measure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% SYNTAX:
%   Out = raster2structdir(RastIn,GridSpecs,DirOffs,nLag,MeasSpecs,iOmni);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2005                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixed parameters
NARGMIN = 6;
TINY = 1/10^10;
VERS = 1;

%% Error checking
if ~(nargin == NARGMIN)
   tit = ['You need ',num2str(NARGMIN),...
          ' input arguments for this function'];
   disp(tit);
   error('Check # of input arguments');
end
[nPix,nCol] = size(RastIn);
checkerrors(nPix,nCol,GridSpecs,DirOffs,nLag,MeasSpecs,iOmni);

%% Get some parameters
icolT = MeasSpecs(:,1);
icolH = MeasSpecs(:,2);
iMeas = MeasSpecs(:,3);
nMeas = numel(iMeas);
nDim  = size(GridSpecs,1);
nDir  = size(DirOffs,1);

%% Determine offsets, grid cell size, and # of grid nodes in rasters
ng = single(GridSpecs(:,1));
switch nDim
  case 1
    ixoff = DirOffs(:,1);
    xsize = GridSpecs(end);
  case 2 
    ixoff = DirOffs(:,1);
    iyoff = DirOffs(:,2);
    xsize = GridSpecs(1,end);
    ysize = GridSpecs(2,end);
  case 3
    ixoff = DirOffs(:,1);
    iyoff = DirOffs(:,2);
    izoff = DirOffs(:,3);
    xsize = GridSpecs(1,end);
    ysize = GridSpecs(2,end);
    zsize = GridSpecs(3,end);
    ng12  = ng(1)*ng(2);
end

%% Compute azimuth and dip angles for the nDir offsets
switch nDim
    case 1
        DirOut = 90*ones(nDir,1);
        DirOut(ixoff<0) = -DirOut(ixoff<0);
        DirOut(ixoff==0) = 0;
        if iOmni == 1
            DirOut = [DirOut; NaN];
        end
    case 2
        azimOut = zeros(nDir,1);
        for id=1:nDir
            % Unit distance along current direction
            %unitdist =   (ixoff(id)*xsize)^2 + (iyoff(id)*ysize)^2;
            %unitdist = sqrt(max(unitdist,0));
            % Compute azimuth and dip angles
            azdir = rad2deg( atan2( ixoff(id)*xsize,iyoff(id)*ysize ) );
            azimOut(id) = round(azdir);
        end
        DirOut = [azimOut 0.001*ones(nDir,1) 0.001*ones(nDir,1)];
        if iOmni == 1
            DirOut = [DirOut; nan(1,3)];
        end
    case 3
        azimOut = zeros(nDir,1);
        dipOut  = zeros(nDir,1);
        for id=1:nDir
            % Unit distance along current direction
            %unitdist =   (ixoff(id)*xsize)^2 + (iyoff(id)*ysize)^2 ...
            %                                 + (izoff(id)*zsize)^2;
            %unitdist = sqrt(max(unitdist,0));
            % Compute azimuth and dip angles
            azdir = rad2deg( atan2( ixoff(id)*xsize,iyoff(id)*ysize ) );
            azimOut(id) = round(azdir);
            xydist  = sqrt((ixoff(id)*xsize)^2 + (iyoff(id)*ysize)^2);
            zdist   = izoff(id)*zsize;
            if xydist > TINY
                  ddir = rad2deg( atan( zdist/xydist ) );
            else
                if izoff == 0
                   error('Found all zero entries in a row of DirOffs');
                end
                if izoff < 0 
                   ddir = -90;
                else
                   ddir = 90;
                end
            end
            dipOut(id) = round(ddir);
        end
        azimOut = [azimOut 0.001*ones(nDir,1) 0.001*ones(nDir,1)];
        dipOut  = [dipOut 0.001*ones(nDir,1) 0.001*ones(nDir,1)];
        DirOut = [azimOut dipOut];
        if iOmni == 1
            DirOut = [DirOut; nan(1,6)];
        end
end %%%%%%%% end switch ndim

%% Declare output structure array
Out.title  = 'Directional sample measures of spatial structure';
inName = inputname(1);  % Report the name of the 1st input argument
Out.structfl = ['Input array with attribute rasters: ',inName];
Out.dir    = single(DirOut);
Out.meas   = single(MeasSpecs);
Out.struct = cell(nDir,1);

%% Loop over # of directions, and compute spatial structure
for id = 1:nDir
  % Declare various arrays
  IvTail   = nan(nLag(id),nMeas,'single');
  IvHead   = nan(nLag(id),nMeas,'single');
  IvMeas   = nan(nLag(id),nMeas,'single');
  Struct   = nan(nLag(id),nMeas,'single');
  Npairs   = zeros(nLag(id),nMeas,'single');
  Tmn      = nan(nLag(id),nMeas,'single');
  Hmn      = nan(nLag(id),nMeas,'single');
  Tvar     = nan(nLag(id),nMeas,'single');
  Hvar     = nan(nLag(id),nMeas,'single');
  ilagsdir = (single(0):single(nLag(id)))';
  onesdir  = ones(size(ilagsdir),'single');
  % Compute starting and ending indices of tail and head rasters 
  % along all directions
  % NOTE: isx_sh -> starting head index along X direction (shifted)
  %       iex_sh -> ending head index along X direction
  %       isx_or -> starting tail index along X direction (original)
  %       iex_or -> ending tail index along X direction
  switch nDim
      case 1
          ix = single(abs(ixoff(id))*ilagsdir);
          if ixoff(id) >= 0
              isx_sh = ix + 1;
              iex_sh = onesdir*ng(1);
              isx_or  = onesdir;
              iex_or  = onesdir*ng(1)-ix;
          else
              isx_sh = onesdir; 
              iex_sh = onesdir*ng(1)-ix;
              isx_or  = ix + 1;  
              iex_or  = onesdir*ng(1);
          end
          % Unit distance along current direction
          unitdist =   (ixoff(id)*xsize)^2;
          unitdist = sqrt(max(unitdist,0));
      case 2
          ix = single(abs(ixoff(id))*ilagsdir);
          if ixoff(id) >= 0
              isx_sh = ix + 1;
              iex_sh = onesdir*ng(1);
              isx_or  = onesdir;
              iex_or  = onesdir*ng(1)-ix;
          else
              isx_sh = onesdir; 
              iex_sh = onesdir*ng(1)-ix;
              isx_or  = ix + 1;  
              iex_or  = onesdir*ng(1);
          end
          iy = single(abs(iyoff(id))*ilagsdir);
          if iyoff(id) >= 0
              isy_sh = iy + 1;
              iey_sh = onesdir*ng(2);
              isy_or  = onesdir;
              iey_or  = onesdir*ng(2)-iy;
          else
              isy_sh = onesdir; 
              iey_sh = onesdir*ng(2)-iy;
              isy_or = iy + 1;  
              iey_or = onesdir*ng(2);
          end
          % Unit distance along current direction
          unitdist =   (ixoff(id)*xsize)^2 + (iyoff(id)*ysize)^2;
          unitdist = sqrt(max(unitdist,0));
      case 3
          ix = single(abs(ixoff(id))*ilagsdir);
          if ixoff(id) >= 0
              isx_sh = ix + 1;
              iex_sh = onesdir*ng(1);
              isx_or  = onesdir;
              iex_or  = onesdir*ng(1)-ix;
          else
              isx_sh = onesdir; 
              iex_sh = onesdir*ng(1)-ix;
              isx_or  = ix + 1;  
              iex_or  = onesdir*ng(1);
          end
          iy = single(abs(iyoff(id))*ilagsdir);
          if iyoff(id) >= 0
              isy_sh = iy + 1;
              iey_sh = onesdir*ng(2);
              isy_or  = onesdir;
              iey_or  = onesdir*ng(2)-iy;
          else
              isy_sh = onesdir; 
              iey_sh = onesdir*ng(2)-iy;
              isy_or = iy + 1;  
              iey_or = onesdir*ng(2);
          end
          iz = single(abs(izoff(id))*ilagsdir);
          if izoff(id) >= 0
              isz_sh = iz + 1;
              iez_sh = onesdir*ng(3);
              isz_or  = onesdir;
              iez_or  = onesdir*ng(3)-iz;
          else
              isz_sh = onesdir; 
              iez_sh = onesdir*ng(3)-iz;
              isz_or = iz + 1;  
              iez_or = onesdir*ng(3);
          end
          % Unit distance along current direction
          unitdist =   (ixoff(id)*xsize)^2 + (iyoff(id)*ysize)^2 ...
                     + (izoff(id)*zsize)^2;
          unitdist = sqrt(max(unitdist,0));
  end %% end switch nDim
  % Distances corresponding to all lags along this direction
  Lagdist  = repmat( ilagsdir.*unitdist, 1, nMeas );
  % Loop over # of lag shifts along this direction
  for il = 1:nLag(id)+1
     % Figure out tail and head indices for current lag
     switch nDim
     case 1
       indTcurr = isx_or(il):iex_or(il);
       indHcurr = isx_sh(il):iex_sh(il);
     case 2
       % Compute all pairs of tail and head indices for X and Y
       [IndTXcurr,IndTYcurr] = meshgrid(isx_or(il):iex_or(il),...
                                        isy_or(il):iey_or(il));
       [IndHXcurr,IndHYcurr] = meshgrid(isx_sh(il):iex_sh(il),...
                                        isy_sh(il):iey_sh(il));
       % Compute indices in GeoEAS array
       indHcurr = (IndHYcurr(:) - 1).*ng(1) + IndHXcurr(:);
       indTcurr = (IndTYcurr(:) - 1).*ng(1) + IndTXcurr(:);
     case 3
       [IndTXcurr,IndTYcurr,IndTZcurr] = meshgrid(isx_or(il):iex_or(il),...
                                                  isy_or(il):iey_or(il),...
                                                  isz_or(il):iez_or(il));
       [IndHXcurr,IndHYcurr,IndHZcurr] = meshgrid(isx_sh(il):iex_sh(il),...
                                                  isy_sh(il):iey_sh(il),...
                                                  isz_sh(il):iez_sh(il));
       indHcurr = (IndHZcurr(:) - 1).*ng12  + ...
                  (IndHYcurr(:) - 1).*ng(1) + IndHXcurr(:);
       indTcurr = (IndTZcurr(:) - 1).*ng12  + ...
                  (IndTYcurr(:) - 1).*ng(1) + IndTXcurr(:);
     end
     for im = 1:nMeas
        % disp(num2str(il))
        IvTail(il,im) = icolT(im);
        IvHead(il,im) = icolH(im);
        IvMeas(il,im) = iMeas(im);
        % vrHead(indH,H)
        headCurr = single(RastIn(indHcurr,icolH(im)));
        % vrTail(indT,T)
        tailCurr = single(RastIn(indTcurr,icolT(im)));
        if iMeas(im) == 2  %%% Cross-semivariogram
            % vrHead(indH,H) - vrHead(indT,H)
            headCurr = headCurr - single(RastIn(indTcurr,icolH(im)));
            % vrTail(indH,T) - vrTail(indT,T)
            tailCurr = single(RastIn(indHcurr,icolT(im))) - tailCurr;
        end
        % Compute structureal measure
        tmp = vrtvrh2struct(tailCurr,headCurr,iMeas(im));
        Struct(il,im) = tmp(1);
        Npairs(il,im) = tmp(2);
        Tmn(il,im)    = tmp(3);
        Hmn(il,im)    = tmp(4);
        Tvar(il,im)   = tmp(5);
        Hvar(il,im)   = tmp(6);
     end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over measures
  end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over lags
  % Construct output array for this direction
  Out.struct{id} = cat(3,Lagdist,Struct,Npairs,IvTail,IvHead,IvMeas,...
                       Tmn,Hmn,Tvar,Hvar);
  Out.struct{id} = permute(Out.struct{id},[1 3 2]);
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over directions
clear Lagdist Struct Npairs IvTail IvHead IvMeas Tmn Hmn Tvar Hvar;

%% Account for omni-directional averaging
if iOmni == 1
    nLagAll = nLag(1)+1;
    % 1st pass to compute total # of pairs
    NpairsTot = zeros(nLagAll,nMeas);
    for im=1:nMeas
        for id=1:nDir
            %size( Out.struct{id}(:,3,im) )
            NpairsTot(:,im) = NpairsTot(:,im) + Out.struct{id}(:,3,im);
        end
    end
    % 2nd pass to compute weighted sums
    Tmp = zeros(nLagAll,10,nMeas,'single');
    for im=1:nMeas
        Tmp(:,3,im) = NpairsTot(:,im);
        Tmp(:,4,im) = icolT(im);
        Tmp(:,5,im) = icolH(im);
        Tmp(:,6,im) = iMeas(im);
        for id=1:nDir
            % Weights = direction # of pairs / total # of pairs
            w =  Out.struct{id}(:,3,im)./NpairsTot(:,im);
            % Distance
            Tmp(:,1,im) = Tmp(:,1,im) + w.*Out.struct{id}(:,1,im);
            % Tail mean
            Tmp(:,7,im) = Tmp(:,7,im) + w.*Out.struct{id}(:,7,im);
            % Head mean
            Tmp(:,8,im) = Tmp(:,8,im) + w.*Out.struct{id}(:,8,im);
            % Tail variance
            mT = Out.struct{id}(:,7,im); % tail mean
            vT = Out.struct{id}(:,9,im); % tail variance
            % Mean of squares = variance + mean^2
            Tmp(:,9,im) = Tmp(:,9,im) + w.*(vT + mT.^2);
            % Head variance
            mH = Out.struct{id}(:,8,im); % head mean
            vH = Out.struct{id}(:,10,im); % head variance
            % Mean of squares = variance + mean^2
            Tmp(:,10,im) = Tmp(:,10,im) + w.*(vH + mH.^2);
            % Structure
            switch iMeas(im);
                case 0 % non-centered covariance = mean of cross-products
                    Tmp(:,2,im) = Tmp(:,2,im) + w.*Out.struct{id}(:,2,im);
                case 1 % semivariogram
                    Tmp(:,2,im) = Tmp(:,2,im) + w.*Out.struct{id}(:,2,im);
                case 2 % cross-semivariogram
                    Tmp(:,2,im) = Tmp(:,2,im) + w.*Out.struct{id}(:,2,im);
                case 3 % centered covariance
                    % non-centered covariance
                    s = Out.struct{id}(:,2,im) + mT.*mH;
                    Tmp(:,2,im) = Tmp(:,2,im) + w.*s;
                case 4 % correlogram
                    % non-centered covariance
                    s = Out.struct{id}(:,2,im).*sqrt(vT.*vH) + mT.*mH;
                    Tmp(:,2,im) = Tmp(:,2,im) + w.*s;
                case 5 % transition probability
                    % joint prob = cond prop * tail prob = non-centered cov
                    Tmp(:,2,im) = w.*Out.struct{id}(:,2,im).*Out.struct{id}(:,7,im);
            end
        end
    end
    % 3rd pass to adjust mean of squares to variances
    for im=1:nMeas
        % Tail variance
        Tmp(:,9,im) = Tmp(:,9,im) - Tmp(:,7,im).^2;
        % Head variance
        Tmp(:,10,im) = Tmp(:,10,im) - Tmp(:,8,im).^2;
        % Structure
        switch iMeas(im)
            case 3 % Centered covariance
                Tmp(:,2,im) = Tmp(:,2,im) - Tmp(:,7,im).*Tmp(:,8,im);
            case 4 % Correlation
                s = Tmp(:,2,im) - Tmp(:,7,im).*Tmp(:,8,im);
                Tmp(:,2,im) = s./sqrt(Tmp(:,9,im).*Tmp(:,10,im));
            case 5 % transition probability 
                % = joint prob / tail mean
                Tmp(:,2,im) = Tmp(:,2,im)./Tmp(:,7,im);
        end
    end
    Out.struct = [Out.struct; {Tmp}];
end

%% FINISHED
disp(' ');
tit = ['Number of directions considered: ',num2str(nDir)];
disp(tit);
tit = ['Number of structural measures per direction: ',num2str(nMeas)];
disp(tit);
if iOmni == 1
    disp('Directional measure(s) also averaged into');
    disp('   pseudo omni-directional measure(s)');
end
disp(' ');
titVers = ['Finished RASTER2STRUCTDIR: Version #',num2str(VERS)];
disp(titVers);


%% Function for computing structural measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute structural measure from head and tail data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = vrtvrh2struct(vrTail,vrHead,imeas)

%% Return, if there are no pairs for this lag
if isempty(vrHead) || isempty(vrTail)
   out = [NaN 0 NaN NaN NaN NaN];
   return
end

%% Fixed parameter
TINY   = 1/(10^10);

%% Handle missing values
Tmp = vrTail + vrHead;
in = isnan(Tmp);
if ~isempty(in)
    vrTail(in) = [];
    vrHead(in) = [];
end

%% Return, if there are no pairs in current lag
nPairs  = numel(vrTail);
if nPairs == 0
    out = [NaN 0 NaN NaN NaN NaN];
    return
end

%% Compute tail and head mean, and variance
tmean = mean(vrTail);
hmean = mean(vrHead);
tvar = var(vrTail,1);
hvar = var(vrHead,1);

%% Compute appropriate structural measure
switch imeas(1)
    case 0                                  % non-centered covariance
       gam = vrHead.*vrTail;
    case 1                                  % semivariance
       gam = 0.5*( vrHead-vrTail).^2;
    case 2                                  % cross-semivariance
      gam    = 0.5*vrHead.*vrTail;
    case 3                                  % covariance
      gam = (vrHead-hmean).*(vrTail-tmean);
    case 4                                  % correlation coefficient
      if hvar < TINY || tvar < TINY
         gam = zeros(nPairs,1);
      else
         gam = (vrHead-hmean).*(vrTail-tmean);
         gam = gam./(sqrt(hvar)*sqrt(tvar));
      end
    case 5                                  % transition probability
      if tmean < TINY
          gam = zeros(nPairs,1);
      else
          gam = vrHead.*vrTail;
          gam = gam./tmean;
      end
    case 6                                   % general relative
      gam = 0.5*( vrHead-vrTail ).^2;
    case 7                                   % pairwise relative
      gam = 0.5*( vrHead-vrTail ).^2;
      gam = gam./( mean([vrTail vrHead],2).^2 );
    case 8                                   % semimadogram
      gam = 0.5*abs( vrHead-vrTail );
    otherwise
      error('Invalid structural measure specification');
end                                          % end measure identification

%% Report average
gam = mean(gam);

%% Account for definition of general relative variogram
if imeas == 6          % general relative
   denom =  (0.5*mean([tmean hmean])^2);
   gam = gam / denom;
end

%% Construct output array
out = [gam nPairs tmean hmean tvar hvar];


%% Function for error checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to perform error checking
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkerrors(nPix,nCol,GridSpecs,DirOffs,nLag,MeasSpecs,iOmni)

% Check tail and head columns
icolCodes = unique(MeasSpecs(:,1:2));
if any(icolCodes > nCol)
    error('Found tail/head index in MeasSpecs > # of columns in RastIn');
end
if any(icolCodes < 1)
    error('Found tail/head index in MeasSpecs < 1');
end
% Check GridSpecs
if size(GridSpecs,2) ~= 3
   error('You need 3 entries in GridSpecs');
end
% Check dimensions of tail and head rasters
if nPix ~= prod( GridSpecs(:,1) ); 
    error('# of rows in RastIn =/= grid size from GridSpecs');
end
% Check DirOffs
nDim = size(GridSpecs,1);
if size(DirOffs,2) ~= nDim
   error('# of columns in DirOffs =/= # of rows in GridSpecs');
end
nDir = size(DirOffs,1);
if length(nLag) ~= nDir;
   error('Length of nLag =/= # of rows in DirOffs');
end
if nDim == 1 && nDir > 2
   error('In 1D, you cannot have > 2 directions');
end
if nDim == 1 && any(DirOffs == 0)
    error('In 1D, entries of DirOffs cannot be 0');
end
% Check MeasSpecs
if any(MeasSpecs(:,3)>8) || any(MeasSpecs(:,3)<0)
    error('MeasSpecs(:,3), iMeas, must be between 0 and 8');
end
% Check whether requested measures of spatial structure make sense
nMeas = size(MeasSpecs,1);
for im=1:nMeas
    tv = MeasSpecs(im,1);
    hv = MeasSpecs(im,2);
    it = MeasSpecs(im,3);
    titm = ['MeasSpecs row #', num2str(im)];
    if tv == hv && it == 2;
       disp(titm);
       error('Cross-semivariogram between the same attribute');
    end
    if tv ~= hv && it > 5
       disp(titm);
       error('Ambiguous choise of head/tail/structural measure');
    end
end
if ~isscalar(iOmni)
    error('iOmni must be a scalar');
end
if ~(iOmni == 0 || iOmni == 1)
    error('iOmni must be 0 or 1');
end
if numel(unique(nLag)) ~= 1 && iOmni == 1
    error('For omni-directional averaging, nLag must have all same entries');
end
if iOmni == 1 && any(MeasSpecs(:,3)>5)
    error('For iOmni=1, MeasSpecs(:,3), iMeas, must be between 0 and 5');
end