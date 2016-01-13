function [CoeffMat,ModParsMulti] = fitcrosstrans(distIn,GammaHat,...
                                   ModPars,CoeffMat0,iWeight,optimPars)
% Fit auto- & cross-transiograms by iterative weighted least squares
%
%% DESCRIPTION:
%  NOTE: ModPars contains the common structure specs
%        CoeffMat0 contains initial partial sill estimates, assumed
%        to sum to the respective class proportions
%        CoeffMat contains improved versions of CoeffMat0
%        ModParsMulti containts all auto- and cross-transiogram models
%        with auto-sills summing to 1-p, and cross-sill summing to p
%        This is consistent with the specification adopted in 
%        function modelstruct2dir
%
%% SYNTAX:
% [CoeffMat,ModParsMulti] =
% fitcrosstrans(distIn,GammaHat,ModPars,CoeffMat0,iWeight,optimPars)

%% Fixed parameters
VERS = 1;
%TINY = 1/1000000;
NREPMAX = 100;

%% Get some parameters
nLag = numel(distIn);
nStruct = size(ModPars,1);
nVar = size(CoeffMat0,1);
nIterMax = optimPars(1);
tolerMin = optimPars(2);
% Set all partial sills of ModPars to 1s
ModPars(:,2) = 1;
nReport  = max([ 1 min([nIterMax/10 NREPMAX]) ]);
switch size(ModPars,2)
    case 4
        nDim = 1;
    case 6
        nDim = 2;
    case 9
        nDim = 3;
end

%% Compute isotropic parameters
switch nDim
    case 1
        ModParsIso = ModPars;
    case 2
        ModParsIso = ModPars(:,[1 2 4 end]);
    case 3
        ModParsIso = ModPars(:,[1 2 6 end]);
end

%% Pre-compute some arrays which do not change through the iterations
GammaBasic = zeros(nLag,nStruct);
sumGammaBasicSq = zeros(nStruct,1);
for is=1:nStruct
    modCurr = ModParsIso(is,:);
    GammaBasic(:,is) = dist2msemivar(distIn,[nLag 1],modCurr,1);
    sumGammaBasicSq(is) = sum(GammaBasic(:,is).^2);
end
nStruct2Keep = nStruct - 1;
if iWeight == 1
    w = (nLag:-1:1)';
else
    w = ones(nLag,1);
end
indDiag = 1:(nVar+1):nVar*nVar;

%% Iterate
wssOld = inf;
CoeffMatOld = CoeffMat0;
it = 0;
while it <= nIterMax && wssOld > tolerMin
    it = it + 1;
    CoeffMatCurr = CoeffMatOld;
    % Loop over # of structures
    for idrop=1:nStruct
        % Rows to keep from ModPars
        rows2keep = 1:nStruct;
        rows2keep(idrop) = [];
        % Compute cropped LMC for all lags
        GammaCurr = zeros(nVar,nVar,nLag);
        GammaDeficitCurr = zeros(nVar,nVar,nLag);
        for ii=1:nVar
            for jj=1:nVar
                % Cross- (off-diagonal) terms
                if ii~=jj
                    sillsIJ = CoeffMatCurr(ii,jj,:);
                    structIJ = zeros(nLag,1);
                    for ik=1:nStruct2Keep
                        ikeep = rows2keep(ik);
                        % Add contribution of current nested structure
                        structIJ = structIJ + ...
                            sillsIJ(ikeep).*GammaBasic(:,ikeep);
                    end
                    GammaCurr(ii,jj,:) = structIJ; 
                    % structIJ is now forward CROPPED LMC-model for ii,jj
                    % Compute Gamma deficit matrix due to cropped LMC
                    structHatIJ = GammaHat(ii,jj,:);
                    GammaDeficitCurr(ii,jj,:) = ...
                        (structHatIJ(:)- structIJ).*GammaBasic(:,idrop);
                end
            end
        end
        % Auto- (diagonal) terms
        for ii=1:nVar
            % Diagonal terms of GammaCurr were intializaed to 0 above,
            % so including them too, is OK
            tmp = squeeze(GammaCurr(ii,1:nVar,:));
            % Replace ii,ii terms of GammaCurr by sum of cross-terms
            structII = 1-sum(tmp);
            GammaCurr(ii,ii,:) = structII;
            % Compute deficit for ii,ii terms
            structHatII = GammaHat(ii,ii,:);
            GammaDeficitCurr(ii,ii,:) = ...
           (structHatII(:)- structII(:)).*GammaBasic(:,idrop);
        end
        % Compute difference matrix over all lags
        GammaDeficitCurr = sum(GammaDeficitCurr,3);
        % Spectral decomposition
        %[Q,L] = eig(GammaDeficitCurr);
        %eigVals = diag(L);
        %eigVals(eigVals<0) = TINY;
        %GammaDeficitCurr = Q*diag(eigVals)*Q';
        % Current coregionalization matrix for dropped structure
        Ccurr = GammaDeficitCurr./sumGammaBasicSq(idrop);
        % Guard against negative entries
        autoSills = Ccurr(indDiag);
        Ccurr(indDiag) = 10;
        Ccurr(Ccurr<0) = 0;
        Ccurr(indDiag) = autoSills;
        % Update coregionalization matrices for FULL forward LMC
        CoeffMatCurr(:,:,idrop) = Ccurr;
        
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End loop over structures
    % We now have nStruct proposed PSD coregionalization matrices
    
    % Compute misfit from full forward LMC
    wssCurr = 0;
    for ii=1:nVar
        for jj=1:nVar
            % Cross-terms
            if ii ~= jj
                sillsIJ = CoeffMatCurr(ii,jj,:);
                structIJ = zeros(nLag,1);
                for is=1:nStruct
                    % Add contribution of current nested structure
                    structIJ = structIJ + ...
                        sillsIJ(is).*GammaBasic(:,is);
                end
                % Over-write ii,jj portion of GammaCurr
                GammaCurr(ii,jj,:) = structIJ;
                % NOTE: structIJ is now forward FULL LMC-model for ii,jj
                % Compute current weighted sum of squared errors
                structHatIJ = GammaHat(ii,jj,:);
                wssCurr = wssCurr + sum( w.*(structHatIJ(:) - structIJ).^2 );
            end
        end
    end
    % Auto-terms
    for ii=1:nVar
        structHatII = GammaHat(ii,ii,:);
        % Extract all current terms, apart from ii,ii
        i2kp = 1:nVar;
        i2kp(ii) = [];
        tmp = squeeze(GammaCurr(ii,i2kp,:));
        % Replace ii,ii terms by sum of all other ones
        structII = 1-sum(tmp);
        GammaCurr(ii,ii,:) = structII;
        % Update wss
        wssCurr = wssCurr + sum( w.*(structHatII(:) - structII(:)).^2 );
    end
    if rem(it,nReport) == 0
        tit = ['Iteration #',num2str(it)];
        disp(tit);
        tit = ['Weighted sum of squares = ',num2str(wssCurr)];
        disp(tit);
    end
    CoeffMatOld = CoeffMatCurr;
    wssOld = wssCurr;
end

%% Scale auto-coefficients to match sills of cross-coefficients
CoeffMat = CoeffMatCurr;
clear CoeffMatCurr;
SillsOut = sum(CoeffMat,3);
autoSills = SillsOut(indDiag);
autoSills = autoSills(:);
Tmp = SillsOut;
Tmp(indDiag) = 0;
autoSillsTarg = 1 - sum(Tmp,2);
clear Tmp;
deltaAutoSills = autoSillsTarg - autoSills;
for ii=1:nVar
    propSill = CoeffMat(ii,ii,:)/autoSills(ii);
    sillNew = CoeffMat(ii,ii,:) + propSill*deltaAutoSills(ii);
    CoeffMat(ii,ii,:) = sillNew;
end
% Now sills of auto-terms sum to p

%% Scale auto-sills of ModParsAll to sum to 1-p
SillsOut = CoeffMat;
for ii=1:nVar
    sillCurr = CoeffMat(ii,ii,:);
    classPropII = sum(sillCurr);
    propVarCurr = sillCurr./classPropII;
    sillNew = sillCurr + (1-2*classPropII).*propVarCurr;
    SillsOut(ii,ii,:) = sillNew;
end

%% Now construct output ModPars
ModParsMulti = cell(nVar,nVar);
for ii=1:nVar
    for jj=1:nVar
        sillsIJ = SillsOut(ii,jj,:);
        tmp = ModPars;
        tmp(:,2) = sillsIJ(:);
        ModParsMulti{ii,jj} = tmp;
    end
end
disp(' ')
disp('In array CoeffMat: ');
disp('All partial sills sum to p');
disp(' ');
disp('In array ModParsMulti: ');
disp('Auto-(partial)sills sum to 1-p');
disp('Cross-(partial)sills sum to p');

%% FINISHED
disp(' ');
disp(['Finished FITCROSSTRANS: Version #',num2str(VERS)]);