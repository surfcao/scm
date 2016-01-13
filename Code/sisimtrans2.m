function Out = sisimtrans2(CoordsSim,DataIn,icolC,icolV,classLabs,...
      ModPars,globProps,srchSpecs,srchData,simPars,iWeight)
  
% Sequential indicator simulation with combination of transition
% probabilities
% By Guofeng Cao
%% DESCRIPTION: sisimtrans.m
% Function to perform sequential indicator simulation (SISIM) at scattered
% locations or grid nodes of a raster, using various types of transition
% probability fusion. All supports are assumed to be quasi-punctual.

%% INPUTS:
%  CoordsSim = [nx xmin xsiz; ny xmin ysiz]; 
%              nx, xmin, xsiz: # of nodes, minimum x-coordinate,
%              and fixed cell size in x-direction
%              NOTE: Dimensionality is determined by this array, not
%                    by array icolC (see below), since the latter could
%                    be empty in the case of unconditional simulation
%   DataIn   = (nDat x nCol) array with sample data (conditioning data)
%              Missing values, flagged as NaNs in any column, are ignored
%   icolC    = (1 x nDim) array with column # (in DataIn) of sample
%              coordinates
%   icolV    = (1 x 1,1+nPrior) array with column #s (in DataIn) for data
%              col#1: sample class labels (category codes at sample locs)
%              col#2:end: prior probabilities (if any) at sample locations
%  classLabs = (1 x nClass) array with class labels, e.g., [0 1 2]
%              NOTE: The # of classes is determined by this array
%              NOTE: This array must contain the same # of classes as
%                    in DataIn(:,icolV) for conditional simulation
%   ModPars  = (nStruct x 6) array with covariogram parameters for
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
%  globProps = (1 x nClass) array with global (regional) class proportions
%  srchSpecs = (1 x 1,3,6) array with search ellipsoid parameters:
%               1D: [srad1]
%               2D: [sang1 srad1 srad2]
%               3D: [sang1 sang2 sang3 srad1 srad2 srad3]
%               sang1, sang2, sang3 = search ellipsoid angles
%               srad1, srad2, srad3 = search ellipsoid radii
%  srchData  = (1 x 3,4) array [ndMin ndMax nsMax (ndOctMax)] with 
%              data parameters for local search:
%              ndMin = minimum # of nearby sample data for simulation
%                      if fewer than ndMin data are found in the
%                      search neighborhood, no prediction is performed
%                      and the corresponding output value is set to NaN
%              ndMax = maximum # of nearby sample data for simulation
%                      if more than ndMax data are found in the
%                      search neighborhood, only the ndMax CLOSEST
%                      data are used for prediction
%   simPars  = (1 x 2,3) array [nSim seed (nMulti)] with simulation 
%               parameters:
%               nSim  = # of simulated realizations to generate
%               seed  = random number seed
%               nMulti = # of multiple nested grids to consider
%                        NOTE: This is only used for regular grids;
%                              nMulti=0 -> purely random path
%                              nMulti>0 -> stratified random path
%  iweights = different transtion probablity combination methods.
%             0->Naive Bayes
%             1->Naive Bayes with OK (ordinary krigint)weights
%             2->Consensus therory 
%             3-> Tau model with weights=1 (Permanence of odds ratios)
%             4-> Tau model with weights [Krishnan(2008)]
%             5-> Nu model with weights
%             6-> Tau model with OK weights, assuming the nearest neight is unknown [Cao2009]
%             7-> Nu model with OK weights, assuming the nearest neight is unknown 
%             8-> Tau model with correlation coefficient as weights.
%             9-> Nu model using conditional probability as nu weights
%             10->SuperDad method to relax CI assumption.[]
%             11->SupreDad with odds of ratios
%             12->Full rank, truncated multivariate gaussian
%             13->Same as method 4, but using average of the weights.
%             14->tau model with mutual information approximation as
%             weights. 
%   
%% OUTPUTS:
%   Out      = output SINGLE precision array with simulated class labels; 
%              each column contains simulated labels for one realization.
%              Out has size (nx*ny*nSim)
%

%%
espsilon = 10^(-10);
nDim = size(CoordsSim,1);
nClass = numel(classLabs);

%checksearch(srchSpecs,nDim,1,srchData);
iCond = 1;
if isempty(DataIn)
   iCond = 0;
end

%% Determine # of simulation locations
GridSpecs = CoordsSim;
ng    = GridSpecs(:,1);
nLocs = prod(CoordsSim(:,1));

%% Get some other parameters
classLabs = classLabs(:)';
nSim  = simPars(1);
seed  = simPars(2);
if numel(simPars) == 3
   nMulti = simPars(end);
else
   nMulti = 0;
end
ndMin = srchData(1);
ndMax = srchData(2);


%% Construct global cumulative PMF
globProps = globProps(:);
globProps(globProps > 1) = 1;
globProps(globProps < 0) = 0;
globProps = globProps./sum(globProps);

%%  Intialize output array
Out = nan(nLocs,nSim,'single');
if ~isempty(DataIn)
   for i = 1:nSim
      Out(:,i) = DataIn(:,icolV);
   end
end

%% Simulation at grid nodes of a raster
   %% Loop over simulation grid nodes
   %     for n= 1:nSim
   %        titC = ['Simulation', num2str(n)];
   %        disp(titC);
   rand('twister',sum(100*clock));
   seed = 1000*rand(1);
   indPath = randompath(1,ng,seed,nMulti);
   for in = 1:nLocs
      titC = ['Working at node ', num2str(in)];
      if mod(in, 10) == 0
        disp(titC);
      end 
      % Index of current node in GeoEAS array
      indLoc = indPath(in);
      % Skip this location if it is a sample location
      if iCond == 1 && ~isnan(Out(indLoc,1))
         continue
      end
      IndexClose = searchnodes2(Out(:,1),CoordsSim,indLoc,srchSpecs,ndMax);
      ValClose = Out(IndexClose,:);
      nCount = size(IndexClose,1);
      if nCount < ndMin
         Out(indLoc,:) = pmf2simclass(globProps,classLabs,nClass,nSim);
         continue;
      end
      %get the coords of neighbour nodes.
      neighCoords = indgeoeas2coords(IndexClose,CoordsSim);

      %get the coords of center node.
      centerCoords = indgeoeas2coords(indLoc,CoordsSim);
      tpcell = coords2mtrans2(neighCoords,ModPars,globProps,centerCoords);
      for n = 1:nSim
         Cpmf = ones(nClass,1);
         switch (iWeight)
         case 0 %  Naive Bayes
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  Cpmf(c) = tpmat(c,ValClose(i,n))*Cpmf(c);
               end
            end
            Cpmf = Cpmf.*globProps;
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
         case 1 % Naive Bayes with OK weights
            covmat = ones(nCount+1,nCount+1);
            covmat(nCount+1,nCount+1) = 0;
            covd2u = ones(nCount+1,1);
            covmat(1:nCount,1:nCount) = coords2mcov2(neighCoords,ModPars);
            covd2u(1:nCount,1) = coords2mcov2(neighCoords,ModPars,centerCoords);
            weights = inv(covmat)*covd2u;
            % readjust the weights if there is negative weights.
            weights = weights(1:nCount);
            if(min(weights) < 0)
               weights = weights - min(weights);
               weights = weights./sum(weights);
            end
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  Cpmf(c) = tpmat(c,ValClose(i,n))*(weights(i))*Cpmf(c);
               end
            end
            Cpmf = Cpmf.*globProps;
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
         case 2 % Consensus
            covmat = coords2mcov2(neighCoords,ModPars);
            covd2u = coords2mcov2(neighCoords,ModPars,centerCoords);
            weights = inv(covmat)*covd2u;
            rang = max(weights)-min(weights);
            if abs(rang) > espsilon
               weights = (weights-min(weights))/rang;
            end
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  Cpmf(c) = (tpmat(ValClose(i,n),c)/globProps(c))^(weights(i))*Cpmf(c);
               end
            end
            Cpmf = Cpmf.*globProps;
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
         case 3 %Tau model with \tau =1 Krishnan(2008)
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = (x_i/x_0(c))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
         case 4 %Tau model with weights Krishnan(2008)
            covmat = ones(nCount+1,nCount+1);
            covmat(nCount+1,nCount+1) = 0;
            covd2u = ones(nCount+1,1);
            covmat(1:nCount,1:nCount) = coords2mcov2(neighCoords,ModPars);
            covd2u(1:nCount,1) = coords2mcov2(neighCoords,ModPars,centerCoords);
            weights = inv(covmat)*covd2u;
            % readjust the weights if there is negative weights.
            weights = weights(1:nCount);
            if(min(weights) < 0)
               weights = weights - min(weights);
               weights = weights./sum(weights);
            end
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = ((x_i/x_0(c))^weights(i))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;

         case 5 %Nu model with weights
            covmat = ones(nCount+1,nCount+1);
            covmat(nCount+1,nCount+1) = 0;
            covd2u = ones(nCount+1,1);
            covmat(1:nCount,1:nCount) = coords2mcov2(neighCoords,ModPars);
            covd2u(1:nCount,1) = coords2mcov2(neighCoords,ModPars,centerCoords);
            weights = inv(covmat)*covd2u;
            % readjust the weights if there is negative weights.
            weights = weights(1:nCount);
            if(min(weights) < 0)
               weights = weights - min(weights);
               weights = weights./sum(weights);
            end
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = weights(i)*(x_i/x_0(c))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;

         case 6 %Tau model with OK weights Krishnan(2008) cross-validation
            weights = ones(nCount,1);
            tmp = neighCoords(1,:);
            tmpneigh = neighCoords(2:end,:);
            covmat = ones(nCount,nCount);
            covmat(nCount,nCount) = 0;
            covd2u = ones(nCount,1);
            covmat(1:nCount-1,1:nCount-1) = coords2mcov2(tmpneigh,ModPars);
            covd2u(1:nCount-1,1) = coords2mcov2(tmpneigh,ModPars,tmp);
            w = inv(covmat)*covd2u;
            weights(2:nCount) = w(1:nCount-1);
            if(min(weights(2:nCount)) < 0)
               weights(2:nCount) = weights(2:nCount)- min(weights(2:nCount));
               weights(2:nCount) = weights(2:nCount)./sum(weights(2:nCount));
            end
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = ((x_i/x_0(c))^weights(i))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;

         case 7 %Nu model with OK weights Krishnan(2008) cross-validation
            weights = ones(nCount,1);
            tmp = neighCoords(1,:);
            tmpneigh = neighCoords(2:end,:);
            covmat = ones(nCount,nCount);
            covmat(nCount,nCount) = 0;
            covd2u = ones(nCount,1);
            covmat(1:nCount-1,1:nCount-1) = coords2mcov2(tmpneigh,ModPars);
            covd2u(1:nCount-1,1) = coords2mcov2(tmpneigh,ModPars,tmp);
            w = inv(covmat)*covd2u;
            weights(2:nCount) = w(1:nCount-1);
            if(min(weights(2:nCount)) < 0)
               weights(2:nCount) = weights(2:nCount)- min(weights(2:nCount));
               weights(2:nCount) = weights(2:nCount)./sum(weights(2:nCount));
            end
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = weights(i)*(x_i/x_0(c))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
         case 8 % Tau model with correlation coefficient
            weights = ones(nCount,1);
            tmp = neighCoords(1,:);
            tmpneigh = neighCoords(2:end,:);
            weights(2:nCount,1) = coords2mcov2(tmpneigh,ModPars,tmp);
            if(min(weights(2:nCount)) < 0)
               weights(2:nCount) = weights(2:nCount)- min(weights(2:nCount));
               weights(2:nCount) = weights(2:nCount)./sum(weights(2:nCount));
            end
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = ((x_i/x_0(c))^weights(i))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;   
         case 9 % NB with conditional probability
             for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  Cpmf(c) = (tpmat(ValClose(i,n),c)/globProps(c))*Cpmf(c);
               end
            end
            Cpmf = Cpmf.*globProps;
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
         case 10 % SuperDad with Naive Bayes
            globProps = double(globProps(:));
            globProps(globProps > (1-eps)) = 1-eps;
            globProps(globProps < eps) = eps;
            globProps = globProps./sum(globProps);
            thresh = cumsum(globProps(1:nClass-1));
            thresh = -sqrt(2).*erfcinv(2.*thresh);
            thresh = [-inf; thresh; inf];   
            for c = 1:nClass
                Cpmf(c)=tpcell{1,1}(c,ValClose(1,n));
                covTot2=coords2mcov2([centerCoords;neighCoords(1,:)],ModPars);
                yLow2=[thresh(c),thresh(ValClose(1,n))];
                yHigh2=[thresh(c+1),thresh(ValClose(1,n)+1)];
                biprob=mvncdf(yLow2,yHigh2,[0,0],covTot2);
                Cpmf(c)=biprob;
                for i = 2:nCount
                 covTot3= coords2mcov2([centerCoords;neighCoords(1,:);neighCoords(i,:)],ModPars);
                 yLow3 =[thresh(c),thresh(ValClose(1,n)),thresh(ValClose(i,n))];
                 yHigh3=[thresh(c+1),thresh(ValClose(1,n)+1),thresh(ValClose(i,n)+1)];
                 triprob=mvncdf(yLow3,yHigh3,[0,0,0],covTot3);
                 cprob=triprob/biprob;
                 Cpmf(c) = cprob*Cpmf(c);
               end
            end
           % Cpmf = Cpmf.*globProps;
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
         
         case 11 % SuperDad with ratio
            globProps = double(globProps(:));
            globProps(globProps > (1-eps)) = 1-eps;
            globProps(globProps < eps) = eps;
            globProps = globProps./sum(globProps);
            thresh = cumsum(globProps(1:nClass-1));
            thresh = -sqrt(2).*erfcinv(2.*thresh);
            thresh = [-inf; thresh; inf]; 
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
                covTot2=coords2mcov2([centerCoords;neighCoords(1,:)],ModPars);
                yLow2=[thresh(c),thresh(ValClose(1,n))];
                yHigh2=[thresh(c+1),thresh(ValClose(1,n)+1)];
                biprob=mvncdf(yLow2,yHigh2,[0,0],covTot2);
             %  bipmf(c)=(biprob/(globProps(ValClose(1,n))-biprob))/(((globProps(ValClose(1,n))-biprob)/biprob)*globProps(c)/(1-globProps(c)));
                Cpmf(c)=((globProps(ValClose(1,n))-biprob)/biprob)*(globProps(c)/(1-globProps(c)));
                for i = 2:nCount     
                    covTot3= coords2mcov2([centerCoords;neighCoords(1,:);neighCoords(i,:)],ModPars);
                    yLow3 =[thresh(c),thresh(ValClose(1,n)),thresh(ValClose(i,n))];
                    yHigh3=[thresh(c+1),thresh(ValClose(1,n)+1),thresh(ValClose(i,n)+1)];
                    triprob=mvncdf(yLow3,yHigh3,[0,0,0],covTot3);
                    covTot2=coords2mcov2([neighCoords(1,:);neighCoords(i,:)],ModPars);
                    yLow2=[thresh(ValClose(1,n)),thresh(ValClose(i,n))];
                    yHigh2=[thresh(ValClose(1,n)+1),thresh(ValClose(i,n)+1)];
                    biprob_nn=mvncdf(yLow2,yHigh2,[0,0],covTot2);
                    cprob=((biprob_nn-triprob)/triprob)*(biprob/(globProps(ValClose(i,n))-biprob));
                    Cpmf(c) = cprob*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;      
         case 12 % full rank
            globProps = double(globProps(:));
            globProps(globProps > (1-eps)) = 1-eps;
            globProps(globProps < eps) = eps;
            globProps = globProps./sum(globProps);
            thresh = cumsum(globProps(1:nClass-1));
            thresh = -sqrt(2).*erfcinv(2.*thresh);
            thresh = [-inf; thresh; inf];   
            covmat_all = ones(nCount+1,nCount+1);
            covmat_dd=ones(nCount,nCount);
            covmat_all = coords2mcov2([neighCoords;centerCoords],ModPars);
            covmat_dd = coords2mcov2(neighCoords,ModPars);
            yLow=thresh(ValClose(:,n));
            yHigh=thresh(ValClose(:,n)+1);
            mu=zeros(1,nCount+1);
            biprob=mvncdf(yLow',yHigh',mu(1:nCount),covmat_dd);
            for c = 1:nClass
                yLow2=thresh([ValClose(:,n);c]);
                yHigh2=thresh([ValClose(:,n)+1;c+1]);
                triprob=mvncdf(yLow2',yHigh2',mu,covmat_all);
                Cpmf(c)=triprob/biprob;
            end
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
            
            case 13 %Tau model with OK weights Krishnan(2008) cross-validation, it uses average of the weights.      
                
            sumweights= zeros(nCount,1);
            for j = 1:nCount
                weights = ones(nCount,1);
                tmp = neighCoords(j,:);
                tmpneigh = neighCoords([1:j-1,j+1:end],:);
                covmat = ones(nCount,nCount);
                covmat(nCount,nCount) = 0;
                covd2u = ones(nCount,1);
                covmat(1:nCount-1,1:nCount-1) = coords2mcov2(tmpneigh,ModPars);
                covd2u(1:nCount-1,1) = coords2mcov2(tmpneigh,ModPars,tmp);
                w = inv(covmat)*covd2u;
                weights([1:j-1,j+1:end]) = w(1:nCount-1);
               if(min(weights(1:nCount)) < 0)
                   weights(1:nCount) = weights(1:nCount)- min(weights(1:nCount));
               end
                sumweights=sumweights+weights;
            end
            weights(1:nCount)=sumweights./sum(sumweights);
%            sumweights
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = ((x_i/x_0(c))^weights(i))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;
            
            case 14 %tau model with mutual information approximation.
            weights = ones(nCount,1); 
            for j = 2:nCount
               covmat = eye(j);
               covmat2 = ones(j,1);
               tmpneigh=[centerCoords;neighCoords(1:j-1,:)];
               % sigma_22;
               covmat(1:j,1:j) = coords2mcov2(tmpneigh,ModPars);
               % sigma_12
               covmat2(1:j)= coords2mcov2(neighCoords(j,:),ModPars,tmpneigh);
               rho=coords2mcov2(centerCoords,ModPars,neighCoords(j,:));
               rho2=det(covmat-(covmat2*covmat2'));
               if rho~=0&&rho2~=0
                  weights(j)=-log(det(covmat)/rho2)/log(1-rho^2);
               end
            end
            comp = 1-globProps;
            x_0 = comp./globProps;
            for c = 1:nClass
               for i = 1:nCount
                  tpmat = tpcell{i,1};
                  x_i = (1-tpmat(ValClose(i,n),c))/tpmat(ValClose(i,n),c);
                  Cpmf(c) = ((x_i/x_0(c))^weights(i))*Cpmf(c);
               end
            end
            x = Cpmf.*x_0;
            Cpmf = 1./(1+x);
            Cpmf = Cpmf./sum(Cpmf);
            Out(indLoc,n) = pmf2simclass(Cpmf,classLabs,nClass,1);
            continue;            
         end % end switch
             
      end
   end


%% FIXPMF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PmfOut,ViolNum,ViolAve,ViolMax] = fixpmf(PmfIn,nClass,...
      ViolNum,ViolAve,ViolMax)
% Correct PMF values & keep track of PMF violations
%
%% DESCRIPTION: fixpmf.m
% Function to post-process a set of input PMF values (pmfIn) and enforce
% its validity, along with recording any violations.
%
%% SYNTAX:
%   [PmfOut,ViolNum,ViolAve,ViolMax] = fixpmf(PmfIn,nClass,...
%                                             ViolNum,ViolAve,ViolMax);
%
%% INPUTS:
%   PmfIn     = (nClass x nSim) array with original PMF values
%   nClass    = scalar with # of rows in PmfIn
%   ViolNum   = (nSim x nClass) array with # of violating sites
%   ViolAve   = (nSim x nClass) array with average magnitude of violations
%   ViolMax   = (nSim x nClass) array with max magnitude of violations
%
%% OUTPUTS:
%   PmfOut    = (nClass x 1) array with corrected PMF values
%   ViolNum   = updated version of input ViolNum
%   ViolAve   = updated version of input ViolAve
%   ViolMax   = updated version of input ViolMax
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% SYNTAX:
%   [PmfOut,ViolNum,ViolAve,ViolMax] = fixpmf(PmfIn,nClass,...
%                                             ViolNum,ViolAve,ViolMax);

%% CREDITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                            November 2005                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Declare some arrays
PmfOut = PmfIn; % nCat x nSim array

%% Enforce [0,1] interval
PmfOut(PmfIn < 0) = 0;
PmfOut(PmfIn > 1) = 1;

%% Row vector with sum of each column of PmfOut
sumPmf = sum(PmfOut);
sumPmf(sumPmf <= 0) = 1;

%% Accumulate violation stats across realizations
for ic = 1:nClass
   PmfOut(ic,:) = PmfOut(ic,:)./sumPmf; % Enforse sum = 1
   indSim = find( PmfOut(ic,:) ~= PmfIn(ic,:) );
   if ~isempty(indSim)
      ViolNum(indSim,ic) = ViolNum(indSim,ic) + 1;
      viol = abs(PmfOut(ic,indSim)-PmfIn(ic,indSim))';
      ViolAve(indSim,ic) = ViolAve(indSim,ic) + viol;
      ViolMax(indSim,ic) = max(ViolMax(indSim,ic),viol);
   end
end

%% PMF2SIMCLASS
function out = pmf2simclass(PmfIn,classLabs,nClass,nSim)
% Simulated class labels from PMF values
%
%% DESCRIPTION: pmf2simclass.m
% Function to generate a set of (1 x nSim) class labels from a set of
% (nClass x nSim) input PMF values
% NOTE: No error checking is performed
% NOTE: It is assumed that seed of rand is fixed outside this function
%
%% SYNTAX:
%   out = pmf2simclass(PmfIn,classLabs,nClass,nSim);
%
%% INPUTS:
%   PmfIn     = (nClass x 1) or (nClass x nSim) array with PMF values
%   classLabs = (1 x nClass) array with class labels (category codes)
%   nClass    = scalar with # of rows in PmfIn
%   nSim      = scalar with # of realizations
%
%% OUTPUTS:
%   out       = (1 x nSim) array with simulated class labels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% SYNTAX:
%   out = pmf2simclass(PmfIn,classLabs,nClass,nSim);

%% CREDITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                            November 2006                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ensure that input PMF has nSim columns
if size(PmfIn,2) == 1
   PmfIn = PmfIn(:, ones(nSim,1) );
end
% Proceed with computation
PmfIn = cumsum(PmfIn);           % Cumulative version of PmfIn
% Generate simulated values from cumulative PMFs
r = rand(1,nSim);              % Array of uniform deviates
r = r( ones(nClass,1) ,:);   % Replicate nClass rows
% Subscripts satisfying random # <= cumulative PMFs
[ii,jj] = find(r <= PmfIn);
% 1st instance of 1 in each column of matrix r<=Cpmf
t = logical(diff([0; jj]));
ii = ii(t);
% Extract appropriate class labels
out = classLabs(ii);
