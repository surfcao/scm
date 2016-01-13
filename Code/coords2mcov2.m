function Out = coords2mcov2(Coords1,ModPars,Coords2)

%% Get dimensionality
nDim=0;
switch size(ModPars,2)
  case 4
      nDim = 1;
  case 6
      nDim = 2;
  case 9
      nDim = 3;
end
if nDim <=0
    disp('ModPars is not correctly specified!');
end


%% Check for consistency in entries of ModPars
%checkmodelstruct(ModPars,nDim);

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

nStruct  = size(ModPars,1);
nD1= size(Coords1,1);

%% Single set of input coordinates
if nargin == 2
    Out = zeros(nD1,nD1);
    for ii=1:nD1
        % Loop over # of dimensions and compute vector components
        ndCurr = nD1-ii+1;
        Vcomps = zeros(ndCurr,nDim);
        for id = 1:nDim
            Vcomps(:,id) = Coords1(ii,id) - Coords1(ii:end,id);
        end
        % Loop over # of nested structures
        covTot = zeros(ndCurr,1);
        for ist=1:nStruct
            % Account for anisotropy
            if iAnis(ist) == 1
                % Compute rotated vector components for current structure
                dAnis = Vcomps*RotMat(:,:,ist);
                % Compute anisotropic distance
                dAnis = sqrt( sum(dAnis.^2,2) );
            else
                % Compute regular distance
                dAnis = sqrt( sum(Vcomps.^2,2) );
            end
            % Compute nested covariance contribution
            covTot = covTot + dist2mcov2(dAnis,[ndCurr 1],...
                              ModPars(ist,1),ModPars(ist,2),...
                              ModPars(ist,4),ModPars(ist,end));
        end
        % Fill-in output covariance array
        Out(ii,ii:end) = covTot;
        Out(ii:end,ii) = covTot';
    end
    return
end

%% Two sets of input coordinates
nD2= size(Coords2,1);
Out = zeros(nD1,nD2);
if nD1 <= nD2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop over all locations of set #1 (row-wise loop)
    for ii=1:nD1
        % Loop over # of dimensions and compute vector components
        Vcomps = zeros(nD2,nDim);
        for id = 1:nDim
            Vcomps(:,id) = Coords1(ii,id) - Coords2(:,id);
        end
        % Loop over # of nested structures
        covTot = zeros(nD2,1);
        for ist=1:nStruct
            % Account for anisotropy
            if iAnis(ist) == 1
                % Compute rotated vector components for current structure
                dAnis = Vcomps*RotMat(:,:,ist);
                % Compute anisotropic distance
                dAnis = sqrt( sum(dAnis.^2,2) );
            else
                % Compute regular distance
                dAnis = sqrt( sum(Vcomps.^2,2) );
            end
            % Compute nested covariance contribution
            covTot = covTot + dist2mcov2(dAnis,[nD2 1],...
                              ModPars(ist,1),ModPars(ist,2),...
                              ModPars(ist,4),ModPars(ist,end));
        end
        % Fill-in output covariance array
        Out(ii,:) = covTot;
    end
    
else  %%% nD1 > nD2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop over all locations of set #2 (column-wise loop)
    for jj=1:nD2
        % Compute vector components
        Vcomps = zeros(nD1,nDim);
        for id = 1:nDim
            Vcomps(:,id) = Coords1(:,id) - Coords2(jj,id);
        end
        % Loop over # of nested structures
        covTot = zeros(nD1,1);
        for ist=1:nStruct
            % Account for anisotropy
            if iAnis(ist) == 1
                % Compute rotated vector components for current structure
                dAnis = Vcomps*RotMat(:,:,ist);
                % Compute anisotropic distance
                dAnis = sqrt( sum(dAnis.^2,2) );
            else
                % Compute regular distance
                dAnis = sqrt( sum(Vcomps.^2,2) );
            end
            % Compute nested covariance contribution
            covTot = covTot + dist2mcov2(dAnis,[nD1 1],...
                              ModPars(ist,1),ModPars(ist,2),...
                              ModPars(ist,4),ModPars(ist,end));
        end
        % Fill-in output covariance array
        Out(:,jj) = covTot;
    end
    
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch nD1 <= nD2


function C = dist2mcov2(D,sizeD,ctype,csill,crange,cextra)
% Isotropic model covariance values for given distances 
% SIMILAR to dist2mcov, copied here for speed

%% Declare output
C = zeros(sizeD);

%% Evaluate model covariogram
switch ctype
  case 0 %Nugget model
      C = C + csill*(D==0);
  case 1 %Spherical model
      Dr = min(D/crange,1);
      C = C + csill*(1.0 - (1.5*Dr-0.5*Dr.^3) );
  case 2 %Exponential model
      C = C + csill*( exp(-3*D/crange) );
  case 3 %Gaussian model
      C = C + csill*( exp ( -3*(D/crange).^2 ) );
  case 4 %Power model
      C = C + 10^10 - crange*(D.^cextra);
  case 5 %Cosine model
      C = C + csill*( cos(2*pi*D/crange) );
  case 6 %Dampened cosine model
      C = C + csill*( exp(-3.0*D./cextra) .* cos(2*pi*D/crange) );
  case 7 %Cardinal sine model
      Dr = 2*pi*D/crange;
      Cr = csill*( sin(Dr)./max(1/10^10,Dr) );
      Cr(Dr <= 1/10^10) = csill;
      C    = C + Cr;
  case 8 %Logarithmic model
      X = sqrt(D.^2 + cextra^2)./crange;
      C = C + 10^10 - csill*(log(X));
  case 9 %Whittle's correlation
      Dr = D/crange;
      Ctmp = csill.*Dr.*besselk(1,Dr);
      Ctmp(isnan(Ctmp)) = csill;
      C = C + Ctmp;
  case 10 %Cubic model
     Dr = min(D/crange,1);
     C = C + csill*(1.0 - 7*Dr.^2 + 35/4*Dr.^3 ...
                             - 7/2*Dr.^5 + 3/4*Dr.^7  );
  case 11 %Generalized Cauchy model
     C = C + csill*( 1 + (D/crange).^2 ).^(-cextra);
  case 12 %stable model
     C = C + csill*( exp ( -3*(D/crange).^cextra ) );
end