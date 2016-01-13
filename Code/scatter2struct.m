function Out=scatter2struct(DataIn,icolC,DirSpecs,LagMidTol,MeasSpecs)
% Sample measures of spatial structure from scattered data
%
%% DESCRIPTION: scatter2struct.m
% Function to calculate sample directional or omni-directional measures of 
% spatial structure (auto-, cross-semivariograms/covariograms/correlograms) 
% from scattered data on possibly multiple attributes. The different
% attributes need not be sampled at the same locations, and hence NaNs
% are allowed in the attribute columns, but not in the coordinate columns.
% NOTE: This function is based on the GSLIB program gamv.
% NOTE: This function can accommodate indicator (0/1) variables,
%       but no error checking is performed to ensure that the columns
%       of RastIn contain indeed such 0/1 entries.
% NOTE: GSLIB's gamv has a BUG, in that a maximum squared distance dismxs
%       is used to quickly discard data pairs; that dismxs is defined as:
%       ((nLags + 0.5)*lagMid)^2. This statement, however, assumes that 
%       the lag tolerance is 0.5 times the lag mid point and hence ignores 
%       the actual input lag tolerance.
%       As a result, gamv returns fewer pairs for large distance lags
%       (close to nLag*lagMid) when a lag tolerance is > than 0.5*lagMid.
% NOTE: Distance classes are defined with CLOSED limits, i.e.,
%       as [lagMid-lagTol lagMid+lagTol]
% NOTE: If any lag distance has lower limit <=0, this is set to a very
%       small number and the 0 distance is treated separately
% NOTE: An azimuth tolerance of, say, 22.5 degrees, indicates a 
%       45 degree sector; and similarly for the dip tolerance
%       In general, all tolerance values are HALF tolerances.
% NOTE: If the azimuth tolerance is set to a value >=90, an
%       omni-directional structural measure is computed. Although the
%       output is flagged as omni-directional, no test is conducted
%       to ensure that the horizontal bandwidth is set to a large number.
%       To compute a truly omni-directional measure set the horizontal
%       bandwidth to a very large number, OR specify an EMPTY array 
%       DirSpecs. Similar arguments hold for the 3D case. 
% NOTE: In 3D, dip increases positively DOWNWARDS from horizontal plane.
%
%% SYNTAX:
%   Out = scatter2struct(DataIn,icolC,DirSpecs,LagMidTol,MeasSpecs);
%
%% INPUTS:
%   DataIn   = (nDat x nCol) array with sample coordinates and data values
%              NOTE: Missing values, flagged as NaNs, are not allowed
%                    in coordinate columns (see below)
%   icolC    = (1 x nDim) array with column #s (in DataIn) of sample
%              coordinates
%              NOTE: Dimensionality is determined by this array
%   DirSpecs = (nDir x 1,3,6) array with direction specifications:
%              1D case: DirSpecs = azim1D
%                     azim1D =   0 -> omni-directional case
%                     azim1D =  90 -> Tail -> Head computation
%                     azim1D = -90 -> Head <- Tail computation
%                     e.g., [90 0 -90]'
%              2D case: DirSpecs = [azim azimTol horBand]
%                     azim    = azimuth angle in degrees 
%                              (positive = clockwise from North;0=NS,90=EW)
%                     azimTol = azimuth tolerance in degrees
%                               NOTE: azimTol is HALF the sector tolerance
%                               NOTE: omni-directional measures computed
%                                     by setting azimTol > 90
%                     horBand = horizontal bandwidth
%              3D case: DirSpecs=[azim azimTol horBand dip dipTol verBand]
%                     dip     = dip angle in degrees 
%                              (positive = down from horizontal plane)
%                     dipTol  = dip tolerance in degrees
%                               NOTE: dipTol is HALF the v-sector tolerance
%                               NOTE: The definition of an omni-directional 
%                                     measure depends only on azimTol
%                     verBand = vertical bandwidth
%              NOTE: This array can be EMPTY, in which case an
%                    omni-directional mode is selected
%  LagMidTol = (nDir x 1) cell array with each cell containing a
%              (nLag(idir) x 2) array [lagMid lagTol] with:
%              lagMid = midpoint of distance class
%              lagTol = tolerance for each class
%              NOTE: Distance classses are thus defined as:
%                    [ lagMid-lagTol, lagMid+lagTol ]
%              NOTE: This could also be a regular array NOT a cell array
%                    In this case, the # of directions nDir is determined
%                    by array DirSpecs and the same distances & tolerances 
%                    are used for ALL directions
%  MeasSpecs = (nMeas x 3) array with column #s (in DataIn) for attributes
%              to consider as tail & head variables, and structural 
%              measure type to compute FOR EACH direction:
%              e.g.: [ icolT(1) icolH(1) iMeas(1);
%                      icolT(2) icolH(2) iMeas(2) ]
%                icolT(im) = column # for tail variable
%                icolH(im) = column # for head variable 
%                iMeas(im) = structural measure type
%                      = 0 -> non-centered covariogram
%                      = 1 -> semivariogram
%                      = 2 -> cross-semivariogram
%                      = 3 -> covariogram
%                      = 4 -> correlogram
%                      = 5 -> transition probability diagram (indicator)
%                      = 6 -> general relative semivariogram
%                      = 7 -> pairwise relative semivariogram
%                      = 8 -> semimadogram
%                e.g.: [ 4 4 1 ] -> semivariogram of attribute in col#4
%                      [ 5 5 3 ] -> covariogram of attribute in col#5
%                      [ 4 5 2 ] -> cross-semivariogram between
%                                   attributes in col#1 and col#2
%                      [ 4 4 3 ; 5 5 3; 4 5 3] -> 
%                              covariogram of attribute in col#4,
%                              covariogram of attribute in col#5,
%                              cross-covariogram of attributes in cols#4,5
%
%% OUTPUTS:
%   Out      = Output structure array with 6 fields:
%              Out.title  = text description of array contents
%              Out.datafl = name of input array with attribute data
%                           from which structural measures were computed
%              Out.dim    = scalar with dimensionality
%              Out.dir    = array DirSpecs copied from input
%              Out.meas   = array MeasSpecs copied from input
%              Out.struct = (nDir x 1) cell array, with each cell
%                           containing a (nLag(idir) x 10 x nMeas) array:
%                           col#1:  average distance in distance class
%                                   specified by LagMidTol
%                           col#2:  corresponding sample structure values
%                           col#3:  number of pairs per lag
%                           col#4:  column (in DataIn) with tail variable
%                           col#5:  column (in DataIn) with head variable
%                           col#6:  type of structural measure
%                           col#7:  mean of tail values
%                           col#8:  mean of head values
%                           col#9:  variance of tail values 
%                           col#10: variance of head values
%                           Same column sequence is repeated for different
%                           combinations of tail and head variables and
%                           measure types, thus building the 3rd dimension
%                           of this array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% SYNTAX:
%   Out = scatter2struct(DataIn,icolC,DirSpecs,LagMidTol,MeasSpecs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2006                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixed parameters
TINY = 1/10^20;
NARGMIN = 5;
VERS = 1;
ICOLDIST = 1; ICOLSTRUCT = 2; ICOLPAIRS = 3; ICOLTAIL = 4; ICOLHEAD = 5;
ICOLTYPE = 6; ICOLMEANT = 7; ICOLMEANH = 8; ICOLVART = 9; ICOLVARH = 10;
ICOLUPDATE = [ICOLDIST ICOLSTRUCT ICOLPAIRS ICOLMEANT ICOLMEANH ...
              ICOLVART ICOLVARH];

%% Error checking
if ~(nargin == NARGMIN)
    tit = ['You need ',num2str(NARGMIN),...
           ' input arguments for this function'];
    disp(tit);
    error('Check # of input arguments');
end
checkerrors(DataIn,icolC,DirSpecs,LagMidTol,MeasSpecs);

%% Get some parameters
icolT = MeasSpecs(:,1);
icolH = MeasSpecs(:,2);
iMeas = MeasSpecs(:,3);
nMeas = numel(iMeas);
icolC = icolC(:);
nDim  = numel(icolC);
nDat  = single(size(DataIn,1));

%% Find # of directions
if iscell(LagMidTol)
    nDir = size(LagMidTol,1);
    for id=1:nDir
        LagMidTol{id} = single(LagMidTol{id});
    end
else
    if ~isempty(DirSpecs)
        nDir = size(DirSpecs,1);
    else
        nDir = 1;
        switch nDim
            case 1
                DirSpecs = 0;
            case 2
                DirSpecs = [0 inf inf];
            case 3
                DirSpecs = [0 inf inf 0 inf inf];
        end
    end
    LagMidTol2 = cell(nDir,1);
    for id=1:nDir
        LagMidTol2{id} = single(LagMidTol);
    end
    LagMidTol = LagMidTol2;
    clear LagMidTol2;
end

%% Check for NaNs in coordinates and duplicate points
CoordsDat = single(DataIn(:,icolC));
if any(isnan(CoordsDat))
    error('Found NaN in coordinate columns');
end
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
% Single precision coordinates
CoordsDat = single(DataIn(:,icolC));

%% Transform LagMidTol to distance class limits
lagLow = cell(nDir,1);
lagHig = cell(nDir,1);
nLagsDir = zeros(nDir,1);
for id=1:nDir
    lagLow{id} = LagMidTol{id}(:,1) - LagMidTol{id}(:,2);
    lagHig{id} = LagMidTol{id}(:,1) + LagMidTol{id}(:,2);
    nLagsDir(id) = single(numel(lagHig{id}));
end

%% Account for 0 or negative lower limits
for id=1:nDir
    lagLow{id}(lagLow{id}<eps) = TINY;
    if LagMidTol{id}(1,1) == 0
        lagLow{id}(1) = -TINY;
        lagHig{id}(1) = TINY;
    else
        lagLow{id} = [-TINY; lagLow{id}];
        lagHig{id} = [ TINY; lagHig{id}];
        nLagsDir(id) = nLagsDir(id)+1;
    end
end

%% Figure out direction specs
DirSpecs = single(DirSpecs);
DirSpecsOut = DirSpecs;
if nDim == 1
    azim1D = DirSpecs(:,1);
    iOmni = azim1D == 0;
end
if nDim > 1
    azim = DirSpecs(:,1);
    azimTol = DirSpecs(:,2);
    horBand = DirSpecs(:,3);
    iOmni = azimTol >= 90;
    if nDim > 2
        dip = DirSpecs(:,4);
        dipTol = DirSpecs(:,5);
        verBand = DirSpecs(:,6);
    end
end

%% Compute azimuth and dip stuff
if nDim == 1
    azimuth1 = (pi/180)*(90-azim1D);
    % set to 90 the 0 azimuth
    azimuth1(azim1D == 0) = 90;
    uvxazm1 = cos(azimuth1);
    % uvyazm1 = sin(90-azim1D) = 0 for 1D cases
else
    azmuth = (pi/180)*(90-azim); % azimuth from EW
    uvxazm = cos(azmuth); % x-component of azimuth
    uvyazm = sin(azmuth); % y-component of azimuth
    atol = (pi/180)*(azimTol);
    csatol = cos(atol); % x-component of azimTol
    if nDim > 2
        declin = (pi/180)*(90-dip); % declination
        uvzdec = cos(declin); % Vertical component of declination
        uvhdec = sin(declin); % Horizontal component of declination
        dtol = (pi/180)*(dipTol);
        csdtol = cos(dtol);
    end
end

%% Declare output structure array
Out.title  = 'Sample (experimental) measures of spatial structure';
inName = inputname(1);  % Report the name of the 1st input argument
Out.datafl = ['Input array with scattered attribute data: ',inName];
Out.dim    = nDim;
Out.dir    = DirSpecsOut;
Out.meas   = MeasSpecs;
Out.struct = cell(nDir,1);

%% Initialize output array
for id=1:nDir
    Out.struct{id} = zeros(nLagsDir(id),10,nMeas,'single');
end

%% Display progress
disp(' ');
disp('Computing requested measure(s) of spatial structure...');

%% 1D case
if nDim == 1
% Loop over # of sample data and compute
for ii=1:nDat
    % Distance between current datum and all others > i
    CoordsCurr = CoordsDat(ii:nDat,1);
    coords1 = CoordsDat(ii,1);
    % DX is defined as x(j) - x(i), since head location x(j) is
    % defined at a vector h away from tail location x(i)
    dx = CoordsCurr - coords1;
    distCurr = abs( dx );
    i0 = find(distCurr < TINY);
    % Loop over # of directions
    for id=1:nDir
        distCurr(i0) = 1;
        dcazm = dx.*uvxazm1(id)./distCurr;
        dcazm(i0) = 1;
        distCurr(i0) = 0;
        % Loop over # of lags
        for il=1:nLagsDir(id)
            jj = find(   distCurr >= lagLow{1}(il) ...
                       & distCurr <= lagHig{1}(il) );
            if isempty(jj)
                continue
            end
            distLag = distCurr(jj);
            dcazmLag = dcazm(jj);
            iNeg = dcazmLag < 0;
            jj = ii + jj - 1;
            ii2 = repmat(ii,numel(distLag),1);
            % Loop over # of measures
            for im=1:nMeas
                vrHead = DataIn(jj,icolH(im));  % vr2(s+h)
                vrTail = DataIn(ii2,icolT(im)); % vr1(s)
                % Exchange ii and jj to account for -dx
                vrHead(iNeg) = DataIn(ii2(iNeg),icolH(im));
                vrTail(iNeg) = DataIn(jj(iNeg),icolT(im));
                Tmp = vrTail + vrHead;
                inonan = find(~isnan(Tmp));
                if isempty(inonan)
                    continue
                end
                vrHead = vrHead(inonan);
                vrTail = vrTail(inonan);
                 % Account for cross-variogram OR omni-directional
                vrTailPr = [];
                vrHeadPr = [];
                if iMeas(im) == 2 || iOmni(id) == 1
                    % Reverse ii and jj to get what was left out
                    vrTailPr = DataIn(jj,icolT(im));  % vr1(s+h)
                    vrHeadPr = DataIn(ii2,icolH(im)); % vr2(s)
                    vrTailPr(iNeg) = DataIn(ii2(iNeg),icolT(im));
                    vrHeadPr(iNeg) = DataIn(jj(iNeg),icolH(im));
                    Tmp = vrTailPr + vrHeadPr;
                    inonan = find(~isnan(Tmp));
                    if isempty(inonan)
                        continue
                    end
                    vrTailPr = vrTailPr(inonan);
                    vrHeadPr = vrHeadPr(inonan);
                end
                % Compute structural measure
                outCurr = vrTH2struct(distLag,vrTail,vrHead,...
                          iMeas(im),iOmni(id),vrTailPr,vrHeadPr,TINY);
                Out.struct{id}(il,ICOLUPDATE,im) = ...
                Out.struct{id}(il,ICOLUPDATE,im) + outCurr;
            end
        end
    end
end
end  %%%%%% 1D case

%% 2D case
if nDim == 2
% Loop over # of sample data and compute
for ii=1:nDat
    % Distance between current datum and all others > i
    CoordsCurr = CoordsDat(ii:nDat,:);
    coords1 = CoordsDat(ii,:);
    % dx = x(j) - x(i)
    dx = CoordsCurr(:,1) - coords1(1);
    dy = CoordsCurr(:,2) - coords1(2);
    distCurr = sqrt(dx.^2 + dy.^2);
    %distCurr = max(distCurr,0);
    i0 = find(distCurr < TINY);
    % Loop over # of directions
    for id=1:nDir
        distCurr(i0) = 1;
        % Inner product between unit directional vector and dx dy
        % yields angle between vector with dx,dy and current direction 
        dcazm = ( dx.*uvxazm(id) + dy.*uvyazm(id) )./distCurr;
        dcazm(i0) = 1;
        distCurr(i0) = 0;
        dcazm( abs(dcazm) < csatol(id) ) = NaN;
        if all(isnan(dcazm)) == 1
            %disp('test aTol');
            continue
        end
        % Outter product -> perpendicular vector
        band = uvxazm(id).*dy - uvyazm(id).*dx;
        dcazm( abs(band) > horBand(id) ) = NaN;
        clear band;
        if all(isnan(dcazm)) == 1
            %disp('test hband');
            continue
        end
        % Finished azimuth checks
        distDir = distCurr;
        distDir(isnan(dcazm)) = NaN;
        % Loop over # of lags
        for il=1:nLagsDir(id)
            jj = find(   distDir >= lagLow{id}(il) ...
                       & distDir <= lagHig{id}(il) );
            if isempty(jj)
                continue
            end
            distLag = distDir(jj);
            dcazmLag = dcazm(jj);
            iNeg = dcazmLag < 0;
            jj = ii + jj - 1;
            ii2 = repmat(ii,numel(distLag),1);
            % Loop over # of measures
            for im=1:nMeas
                vrHead = DataIn(jj,icolH(im));
                vrTail = DataIn(ii2,icolT(im));
                vrHead(iNeg) = DataIn(ii2(iNeg),icolH(im));
                vrTail(iNeg) = DataIn(jj(iNeg),icolT(im));
                Tmp = vrTail + vrHead;
                inonan = find(~isnan(Tmp));
                if isempty(inonan)
                    continue
                end
                vrHead = vrHead(inonan);
                vrTail = vrTail(inonan);
                % Account for cross-variogram OR omni-directional
                vrTailPr = [];
                vrHeadPr = [];
                if iMeas(im) == 2 || iOmni(id) == 1
                    vrTailPr = DataIn(jj,icolT(im));
                    vrHeadPr = DataIn(ii2,icolH(im));
                    vrTailPr(iNeg) = DataIn(ii2(iNeg),icolT(im));
                    vrHeadPr(iNeg) = DataIn(jj(iNeg),icolH(im));
                    Tmp = vrTailPr + vrHeadPr;
                    inonan = find(~isnan(Tmp));
                    if isempty(inonan)
                        continue
                    end
                    vrTailPr = vrTailPr(inonan);
                    vrHeadPr = vrHeadPr(inonan);
                end
                % Compute structural measure
                outCurr = vrTH2struct(distLag,vrTail,vrHead,...
                          iMeas(im),iOmni(id),vrTailPr,vrHeadPr,TINY);
                Out.struct{id}(il,ICOLUPDATE,im) = ...
                Out.struct{id}(il,ICOLUPDATE,im) + outCurr;
            end
        end
    end
end
end  %%%%%% 2D case

%% 3D case
if nDim == 3
% Loop over # of sample data and compute
for ii=1:nDat
    % Distance between current datum and all others > i
    CoordsCurr = CoordsDat(ii:nDat,:);
    coords1 = CoordsDat(ii,:);
    % dx = x(j) - x(i)
    dx = CoordsCurr(:,1) - coords1(1);
    dy = CoordsCurr(:,2) - coords1(2);
    dz = CoordsCurr(:,3) - coords1(3);
    distCurr = sqrt(dx.^2 + dy.^2 + dz.^2);
    %distCurr = max(distCurr,0);
    i0H = find(distCurr < TINY);
    dxy = sqrt(dx.^2 + dy.^2);
    %dxy = max(dxy,0);
    i0XY = find(dxy < TINY);
    % Loop over # of directions
    for id=1:nDir
        distCurr(i0H) = 1;
        dxy(i0XY) = 1;
        % Azimuth computations
        dcazm = ( dx.*uvxazm(id) + dy.*uvyazm(id) )./dxy;
        dcazm(i0XY) = 1;
        dxy(i0XY) = 0;
        dcazm( abs(dcazm) < csatol(id) ) = NaN;
        if all(isnan(dcazm)) == 1
            %disp('test aTol');
            continue
        end
        band = uvxazm(id).*dy - uvyazm(id).*dx;
        dcazm( abs(band) > horBand(id) ) = NaN;
        clear band;
        if all(isnan(dcazm)) == 1
            %disp('test hband');
            continue
        end
        % Dip computations
        dxyDir = dxy;
        dxyDir(dcazm<0) = -dxyDir(dcazm<0);
        dcdec = ( dxyDir.*uvhdec(id) + dz.*uvzdec(id) )./distCurr;
        dcdec(i0H) = 0;
        distCurr(i0H) = 0;
        dcdec( abs(dcdec) < csdtol(id) ) = NaN;
        if all(isnan(dcdec)) == 1
            %disp('test dTol');
            continue
        end
        band = uvhdec(id).*dz - uvzdec(id).*dxyDir;
        dcdec( abs(band) > verBand(id) ) = NaN;
        clear band;
        if all(isnan(dcdec)) == 1
            %disp('test vband');
            continue
        end
        % Finished azimuth and dip checks
        distDir = distCurr;
        distDir(isnan(dcazm)) = NaN;
        distDir(isnan(dcdec)) = NaN;
        % Loop over # of lags
        for il=1:nLagsDir(id)
            jj = find(   distDir >= lagLow{id}(il) ...
                       & distDir <= lagHig{id}(il) );
            if isempty(jj)
                continue
            end
            distLag = distDir(jj);
            dcazmLag = dcazm(jj);
            dcdecLag = dcdec(jj);
            iNeg = dcazmLag < 0 | dcdecLag < 0;
            jj = ii + jj - 1;
            ii2 = repmat(ii,numel(distLag),1);
            % Loop over # of measures
            for im=1:nMeas
                vrHead = DataIn(jj,icolH(im));
                vrTail = DataIn(ii2,icolT(im));
                vrHead(iNeg) = DataIn(ii2(iNeg),icolH(im));
                vrTail(iNeg) = DataIn(jj(iNeg),icolT(im));
                Tmp = vrTail + vrHead;
                inonan = find(~isnan(Tmp));
                if isempty(inonan)
                    continue
                end
                vrHead = vrHead(inonan);
                vrTail = vrTail(inonan);
                % Account for cross-variogram OR omni-directional
                vrTailPr = [];
                vrHeadPr = [];
                if iMeas(im) == 2 || iOmni(id) == 1
                    vrTailPr = DataIn(jj,icolT(im));
                    vrHeadPr = DataIn(ii2,icolH(im));
                    vrTailPr(iNeg) = DataIn(ii2(iNeg),icolT(im));
                    vrHeadPr(iNeg) = DataIn(jj(iNeg),icolH(im));
                    Tmp = vrTailPr + vrHeadPr;
                    inonan = find(~isnan(Tmp));
                    if isempty(inonan)
                        continue
                    end
                    vrTailPr = vrTailPr(inonan);
                    vrHeadPr = vrHeadPr(inonan);
                end
                % Compute structural measure
                outCurr = vrTH2struct(distLag,vrTail,vrHead,...
                          iMeas(im),iOmni(id),vrTailPr,vrHeadPr,TINY);
                Out.struct{id}(il,ICOLUPDATE,im) = ...
                Out.struct{id}(il,ICOLUPDATE,im) + outCurr;
            end
        end
    end
end
end %%%% 3D case

%% Loop over # of directions, and # of measures and report averages
for id=1:nDir
    for im=1:nMeas
        Out.struct{id}(:,ICOLTAIL,im) = icolT(im);
        Out.struct{id}(:,ICOLHEAD,im) = icolH(im);
        Out.struct{id}(:,ICOLTYPE,im) = iMeas(im);
        nPrs = Out.struct{id}(:,ICOLPAIRS,im);
        i0 = find(nPrs == 0);
        nPrs(i0) = 1;
        % Distance
        meanD = Out.struct{id}(:,ICOLDIST,im);
        meanD = meanD./nPrs;
        meanD(i0) = NaN;
        Out.struct{id}(:,ICOLDIST,im) = meanD;
        % Structure
        gam = Out.struct{id}(:,ICOLSTRUCT,im);
        gam = gam./nPrs;
        % Tail mean
        meanT = Out.struct{id}(:,ICOLMEANT,im);
        meanT = meanT./nPrs;
        meanT(i0) = NaN;
        Out.struct{id}(:,ICOLMEANT,im) = meanT;
        % Head mean
        meanH = Out.struct{id}(:,ICOLMEANH,im);
        meanH = meanH./nPrs;
        meanH(i0) = NaN;
        Out.struct{id}(:,ICOLMEANH,im) = meanH;
        % Tail variance
        varT = Out.struct{id}(:,ICOLVART,im);
        varT = varT./nPrs;
        varT = varT - meanT.^2;
        varT(i0) = NaN;
        Out.struct{id}(:,ICOLVART,im) = varT;
        % Head variance
        varH = Out.struct{id}(:,ICOLVARH,im);
        varH = varH./nPrs;
        varH = varH - meanH.^2;
        varH(i0) = NaN;
        Out.struct{id}(:,ICOLVARH,im) = varH;
        switch iMeas(im)
            %%% Nothing to do for non-centered covariogram
            case 1 %%% Semivariogram
                gam = 0.5*gam;
            case 2 %%% Cross-semivariogram
                gam = 0.5*gam;
            case 3 %%% Covariogram
                gam = gam - meanT.*meanH;
            case 4 %%% Correlogram
                varTH = varT.*varH;
                ind0 = find( varTH < TINY );
                varTH(ind0) = 1;
                gam = (gam - meanT.*meanH)./sqrt(varTH);
                gam(ind0) = 0;
            case 5 %%% Transition probability diagram
                ind0 = find( meanT < TINY );
                meanTtmp = meanT;
                meanTtmp(ind0) = 1;
                gam = gam./meanTtmp;
                gam(ind0) = 0;
            case 6 %%% General relative semivariogram
                meanTH = 0.5*(meanT+meanH);
                meanTH = meanTH.*meanTH;
                ind0 = find( meanTH < TINY);
                meanTH(ind0) = 1;
                gam = gam./meanTH;
                gam(ind0) = 0;
            case 7 %%% Pairwise relative semivariogram
                gam = 0.5*gam;
            case 8 %%% Semimadogram
                gam = 0.5*gam;
        end
        gam(i0) = NaN;
        % Report in output
        Out.struct{id}(:,ICOLSTRUCT,im) = gam;
    end
end

%% FINISHED
disp(' ');
tit = ['Number of directions considered: ',num2str(nDir)];
disp(tit);
tit = ['Number of structural measures per direction: ',num2str(nMeas)];
disp(tit);
disp(' ');
titVers = ['Finished SCATTER2STRUCT: Version #',num2str(VERS)];
disp(titVers);


%% Function for computing structural measure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute structural measure from head and tail data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = vrTH2struct(distIn,vrTail,vrHead,iMeas,iOmn,...
               vrTailPr,vrHeadPr,TINY)

%% Compute appropriate structural measure
switch iMeas
    case 0                                  % non-centered covariance
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = sum( vrHead.*vrTail ); % vr2(s+h)*vr1(s)
        if iOmn == 1
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            gam = gam + sum(vrHeadPr.*vrTailPr); % vr2(s)*vr1(s+h)
        end
    case 1                                  % semivariogram
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = sum( (vrHead - vrTail).^2 ); % vr2(s+h)-vr1(s)
        if iOmn == 1
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            gam = gam + sum( (vrHeadPr - vrTailPr).^2 ); % vr2(s)-vr1(s+h)
        end
    case 2                                  % cross-semivariogram
        nP = numel(vrTail);
        meanD = sum(distIn);
        vrTail2 = vrTail+vrTailPr;
        vrHead2 = vrHead+vrHeadPr;
        meanT = 0.5*sum(vrTail2);
        meanH = 0.5*sum(vrHead2);
        varT = 0.25*sum(vrTail2.^2);
        varH = 0.25*sum(vrHead2.^2);
        gam = sum( (vrHeadPr - vrHead).*(vrTail-vrTailPr) );
        % The above is: vr2(s)-vr2(s+h) * vr1(s)-vr1(s+h)
        %      same as: vr2(s+h)-vr2(s) * vr1(s+h)-vr1(s)
        % For positive dcazm >=0, the above equation is:
        % (vrHead(ii,H) - vrHead(jj,H))*(vrTail(ii,T)-vrTail(jj,T))
        % For negative dcazm <0, indices ii and jj are reversed
        % NOTE: Using (vrHead - vrHeadPr)*(vrTailPr-vrTail)
        %       i.e.,  ( jj - ii ) & ( jj - ii ) yields SAME result
    case 3                                  % covariance
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = sum( vrHead.*vrTail ); % vr2(s+h)*vr1(s)
        if iOmn == 1
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            gam = gam + sum(vrHeadPr.*vrTailPr); % vr2(s)*vr1(s+h)
        end
    case 4                                  % correlogram
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = sum( vrHead.*vrTail );
        if iOmn == 1
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            gam = gam + sum(vrHeadPr.*vrTailPr);
        end
    case 5                                  % transition probability
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = sum( vrHead.*vrTail ); % vr2(s+h)*vr1(s)
        if iOmn == 1
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            gam = gam + sum(vrHeadPr.*vrTailPr); % vr2(s)*vr1(s+h)
        end
    case 6                                  % general relative
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = sum( ( vrHead - vrTail ).^2 );
        if iOmn == 1
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            gam = gam + sum( (vrHeadPr - vrTailPr).^2 );
        end
    case 7                                  % pairwise relative
        ind = find(abs(vrTail+vrHead)<TINY);
        vrTail(ind) = [];
        vrHead(ind) = [];
        distIn(ind) = [];
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = 2*(vrTail-vrHead)./(vrTail+vrHead);
        gam = sum(gam.*gam);
        if iOmn == 1
            ind = find(abs(vrTailPr+vrHeadPr)<TINY);
            vrTailPr(ind) = [];
            vrHeadPr(ind) = [];
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            tmp = 2*(vrTailPr-vrHeadPr)./(vrTailPr+vrHeadPr);
            gam = gam + sum(tmp.*tmp);
        end
    case 8                                   % semimadogram
        nP = numel(vrTail);
        meanD = sum(distIn);
        meanT = sum(vrTail);
        meanH = sum(vrHead);
        varT = sum(vrTail.^2);
        varH = sum(vrHead.^2);
        gam = sum( abs(vrHead - vrTail) );
        if iOmn == 1
            nP = nP + numel(vrTailPr);
            meanD = 2*meanD;
            meanT = meanT + sum(vrTailPr);
            meanH = meanH + sum(vrHeadPr);
            varT = varT + sum(vrTailPr.^2);
            varH = varH + sum(vrHeadPr.^2);
            gam = gam + sum( abs(vrHeadPr - vrTailPr) );
        end
end

%% Construct output array
out = [meanD gam nP meanT meanH varT varH];



%% Function for error checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to perform error checking
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkerrors(DataIn,icolC,DirSpecs,LagMidTol,MeasSpecs)

TINY = 1/10^6;
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
nDim = numel(icolC);
% Check DirSpecs
if ~isempty(DirSpecs)
    if nDim == 1
        if size(DirSpecs,2) ~= 1
            error('In 1D, you need 1 column in DirSpecs');
        end
        for ii=1:numel(DirSpecs)
            if ~(DirSpecs(ii)==-90 || DirSpecs(ii)==0 || DirSpecs(ii)==90)
                error('DirSpecs must contain -90, 0 or 90');
            end
        end
        nDir = numel(DirSpecs);
    end
    if nDim > 1
        if nDim == 2
            if size(DirSpecs,2) ~= 3
                error('In 2D, you need 3 columns in DirSpecs');
            end
        else
            if size(DirSpecs,2) ~= 6
                error('In 3D, you need 6 columns in DirSpecs');
            end
        end
        if any(DirSpecs(:,2)<TINY)
            error('You need positive entries in DirSpecs(:,2), azimTol');
        end
        if any(DirSpecs(:,3)<TINY)
            error('You need positive entries in DirSpecs(:,3), horBand');
        end
        if nDim > 2
            if any(DirSpecs(:,5)<TINY)
                error('You need positive entries in DirSpecs(:,5), dipTol');
            end
            if any(DirSpecs(:,6)<TINY)
                error('You need positive entries in DirSpecs(:,6), verBand');
            end
        end
        nDir = size(DirSpecs,1);
    end
else
    nDir = 1;
end
% Check LagMidTol
if iscell(LagMidTol)
    if length(LagMidTol) ~= nDir
        error('# of cells in LagMidTol incompatible with DirSpecs');
    end
    for id=1:length(LagMidTol);
        LagCurr = LagMidTol{id};
        if size(LagCurr,2) ~= 2
            error('You need 2 columns in the cells of LagMidTol');
        end
    end
else
    if size(LagMidTol,2) ~= 2
        error('You need 2 columns in LagMidTol');
    end
end
% Check MeasSpecs
if size(MeasSpecs,2) ~= 3
    error('You need 3 columns in MeasSpecs');
end
icolT = MeasSpecs(:,1);
icolH = MeasSpecs(:,2);
if any(icolT > nCols )
    error('Found entry in MeasSpecs(:,1) > # of columns in DataIn');
end
if any(icolT < 1)
    error('Found entry in MeasSpecs(:,1) < 1');
end
if any(icolH > nCols )
    error('Found entry in MeasSpecs(:,2) > # of columns in DataIn');
end
if any(icolH < 1)
    error('Found entry in MeasSpecs(:,2) < 1');
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
for im=1:nMeas
    tmp = MeasSpecs(im,1:2);
    for ii=1:nDim
        if any(tmp == icolC(ii))
            error('Found entry in MeasSpecs(:,1:2) equal to one in icolC');
        end
    end
end