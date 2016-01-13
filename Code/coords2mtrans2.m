function Out = coords2mtrans2(Coords1,ModPars,globProps,Coords2)

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

nStruct  = size(ModPars,1);
nClass = numel(globProps);

%% Compute Gaussian thresholds
globProps = double(globProps(:));
globProps(globProps > (1-eps)) = 1-eps;
globProps(globProps < eps) = eps;
globProps = globProps./sum(globProps);
thresh = cumsum(globProps(1:nClass-1));
thresh = -sqrt(2).*erfcinv(2.*thresh);
thresh = [-inf; thresh; inf];

nD1= size(Coords1,1);
if nargin == 3
    tpmat = cell(nD1,nD1);
    for ii=1:nD1
        for jj=1:nD1
            Vcomps = zeros(1,nDim);
            for id = 1:nDim
                Vcomps(1,id) = Coords1(ii,id) - Coords1(jj,id);
            end
            covTot = 0;
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
                p = ModPars(ist,[1 2 4 end]);
                covTot = covTot + dist2mcov(dAnis,[1 1],p,1);
            end
            for iT=1:nClass
                for iH=1:nClass
                    yLowT = thresh(iT);
                    yLowH = thresh(iH);
                    yLow = [yLowT yLowH];
                    yHighT = thresh(iT+1);
                    yHighH = thresh(iH+1);
                    yHigh = [yHighT yHighH];
                    jointProb = bigausscdf(yLow,yHigh,covTot);
                    pT = globProps(iT);
                    tpmat{ii,jj}(iT,iH)= jointProb./pT;
                end
            end
        end
    end
    Out=tpmat;
else

%% Two sets of input coordinates
nD2= size(Coords2,1);
tpmat = cell(nD1,nD2);
for ii=1:nD1
     for jj=1:nD2
        Vcomps = zeros(nD2,nDim);
        for id = 1:nDim
            Vcomps(:,id) = Coords1(ii,id) - Coords2(jj,id);
        end
        covTot = 0;
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
            p = ModPars(ist,[1 2 4 end]);
            covTot = covTot + dist2mcov(dAnis,[nD2 1],p,1);
        end
       for iT=1:nClass
           for iH=1:nClass
                yLowT = thresh(iT);
                yLowH = thresh(iH);
                yLow = [yLowT yLowH];
                yHighT = thresh(iT+1);
                yHighH = thresh(iH+1);
                yHigh = [yHighT yHighH];
                jointProb = bigausscdf(yLow,yHigh,covTot);
                pT = globProps(iT);
                tpmat{ii,jj}(iT,iH)= jointProb./pT;
           end
        end
    
     end
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch nD1 <= nD2
  Out=tpmat;
end









