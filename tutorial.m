%% Prediction and Simulation in Categorical Field: A Transition Probability Combination Approach
% This tutorial illustrates the application of a geostatistical approach
% for modeling of categorical spatial data, such as land
% use/land cover class and soil type data, using different methods 
% for fusing spatial transition probabilities. These transition 
% probabilities, which quantify the likelihood of observing a particular 
% class label given that another class label (of the same or different class) has been observed elsewhere, 
% are modeled as functions of distance (transiograms or
% transition probability diagrams), and are used as measures of spatial 
% structure (texture) in categorical fields (in analogy to
% variograms/covariograms in continuous fields).
%
%% Preparation
% In what follows, examples will be given after a brief introduction to
% relevant theory. Data for this tutorial are provided in the 'Data' directory 
% of this package. Please set up the path and load the tutorial
% data before running the example code.
load sample_data.mat
who
%% Categorical Spatial Data
% Categorical spatial data are a common source of information in many 
% scientific disciplines. Examples include land use classes in geography, 
% lithofacies in geology and socioenomoic survey data. Categorical data
% are typically represented by mutually exclusive and collectively 
% exhaustive (MECE) classes and visualized as area-class maps. The class
% labels of such variables may be distinguished only by their attributes
% (nominal, e.g., land use and land cover classes) or ranked orders 
% indicating whether an observation has a higher or lower value than 
% another one (ordinal, e.g., best, better).
%
%% Displaying categorical spatial data: function rastermap
% Example use of the the GeostatsToolbox function |rastermap| for
% displaying categorical atributes on a raster (up to 2D).
%
% GeostatsToolbox function |rastermap| for displaying raster data
%
% Inputs:
%
% * GridSpecs   = [nx xmin xsize; ny ymin ysize] raster grid specifications 
% * Rasterin    = data array with Geo-EAS formatted attribute rasters
% * icolV       = column # (in RasterIn) for attribute raster to display
% * vType       = string with flag for continuous ('cont') or categorical ('cat') attribute raster
% * colorLabels = Optional array with color limits OR category labels
% * plotPars    = Optional (1 x 2) cell array {igray ibar} with flags for gray scale and bar legend to use for plotting
%%
% Use default color, class labels ('1','2','3'), and vertical legend bar to display the categorical map
gridspecs=[300,0.5,1;300,0.5,1];
figure;
rastermap(gridspecs,simcat,1,'cat');
% Use gray, vertical legend bar to display the categorical map with
% categories 'A','B','C'
figure;
rastermap(gridspecs,simcat,1,'cat',{'A','B','C'},{'gray','vbar'});

%% Displaying categorical spatial data: function scattermap
% Example use of the the GeostatsToolbox function |scattermap| for
% displaying categorical attributes at scatter locations(up to 2D).
%
% GeostatsToolbox function: scattermap
%
% Inputs
%
% * DataIn      = (nDat x nCol) input arrray with sample coordinates 
%                 and data values
% * icolC       = (1 x nDim) array with column numbers (in DataIn) with coordinates
% * bullet      = scalar to control size of color-filled circles drawn at sample locations
% * icolV       = Optional scalar with column # (in DataIn) with attribute values to display
% * vType       = string with flag for continuous ('cont') or categorical ('cat') attribute raster
% * colorLabels = Optional array with color limits OR category labels
% * plotPars    = Optional (1 x 2) cell array {igray ibar} with flags
%                 for gray scale and bar legend to use for plotting
figure;
scattermap(coords,[1,2],10,3,'cat');

%% Transiogram: Spatial structure measure for categorical fields
% A transition probability of class occurrence is not a new concept, 
% but it is until recently that it has been regarded as a spatial structure 
% measure and its relationship with the indicator (cross)variogram 
% or indicator (cross)covariance was discussed [Carle and Fogg 1996, 1997].
% Transiograms can be obtained by direct computation (no parametric model 
% involved) from exhaustive sample data or from a probabilistic model
% of an underlying random field, e.g., a truncated multivariate Gaussian
% field. In analogy to mathematical forms of variograms in Gaussian 
% random fields, a transiogram can also be defined in a similar way. 
% The exponential form of transiograms was also discussed in Li (2006), 
% in the absence, however, of a clear underlying theory. The parameters
% of a set of auto- and cross- transiogram models can be obtained by 
% fitting fitting a Linear Model of Coregionalization (LMC) based on
% sample data.
%
%% Model transiograms from a truncated (clipped) bivariate Gaussian field
% Computing theoretical (model) transiograms from a bivariate Gaussian 
% field using function |corrbigauss2indic|.
%
% GeostatsToolbox function |corrbigauss2indic|
%
% Inputs:
%
% * ModPars   = (nStruct x 4,6,9) array with correlogram parameters for
%               a standard Gaussian random function
% * lagDist   = (nDir x 1) cell array, with each cell containing the
%               (nLag(idir)) lag distances for each direction.
% * azims     = (nDir x 1) array with lag azimuth angles for each
%               direction
% * dips      = (nDir x 1) array with lag dip angles for each
%               direction
% * globProps = (1 x nClass) array with global (stationary) class
%               proportions
% * MeasSpecs = (nMeas x 3) array with tail & head class (order i and j 
%               in array globProps) & indicator structural measure type 
%               to compute for each direction
%
% For more details, please type |help corrbigauss2indic|.
classprops=[0.35,0.4,0.25];
modpars=[3,0.999,0,50,50,0;0,0.001,0,50,50,0];
meas_tp= [1 1 5;1 2 5;1 3 5; 2 1 5;2 2 5;2 3 5;3 1 5; 3 2 5; 3 3 5];
lagdist100=[1:100];
% Calculate 
tp=corrbigauss2indic(modpars,{lagdist100},0,[],classprops,meas_tp);
% Plot the resulting transition transiogram out.
figure;
for i=1:3
    for j = 1:3
        subplot(3,3,(i-1)*3+j); 
        plot(tp.struct{1}(:,1,(i-1)*3+j),tp.struct{1}(:,2,(i-1)*3+j));
        xlabel('distance');
        ylabel('transition probability');
        txt=sprintf('Class %d-%d',i,j);
        title(txt);
    end
end
%% Computing transiograms from sample data
% Example of using GeostatToolbox function |raster2structdir| to compute
% transiograms from exhaustive sample data on a regular raster.
%
% GeostatsToolbox function |raster2structdir|
%
% Inputs:
%
% * RastIn    = (nPix x nCol) array with input attribute rasters in GeoEAS format
% * Gridspecs = (nDim x 3) array with grid specifications for RastIn
% * DirOffs   = (nDir x nDim) array [ixoff iyoff izoff] with # of pixels
%               along each direction specifying the magnitude and 
%               direction of each lag vector
% * nLag = (1 x nDir) vector with number of lags per direction,
% * MeasSpecs = (nMeas x 3) array with column #s for tail & head 
%               attribute rasters & structural measure specs to compute
%               for each direction
% * iOmni     = scalar with flag to average directional measures into
%               a single pseudo omni-directional measure (=1) or not
%
% For more details, please type |help raster2structdir|.

 simcat_indic=class2indic(simcat,1);
 nsize =100;
 sim_tp=raster2structdir(simcat_indic,gridspecs,[1 0;0 1],[nsize,nsize],meas_tp,0);
 sim_tp_data=sim_tp.struct{1};
 figure;
 ip=0;
 nplot=3;
 for j=1:nplot
     for k=1:nplot
               ip = ip + 1;
               subplot(nplot,nplot,ip);
               plot(sim_tp_data(1:100,1,(j-1)*nplot+k),sim_tp_data(1:100,2,(j-1)*nplot+k),'-');               
               xlabel('distance');
               ylabel('transition probability');
               txt=sprintf('Class %d-%d',j,k);
               title(txt);
      end
 end   

%%
% Example of using GeostatToolbox function |scatter2struct| to compute
% transiograms from sample data available at scattered locations.
%
% GeostatsToolbox function |scatter2struct|
%
% Inputs:
%
%  * DataIn   = (nDat x nCol) array with sample coordinates and data values
%              NOTE: Missing values, flagged as NaNs, are not allowed in coordinate columns (see below)
%  * icolC    = (1 x nDim) array with column #s (in DataIn) of sample
%              coordinates. Dimensionality is determined by this array
%  * DirSpecs = (nDir x 1,3,6) array with direction specifications
%              NOTE: This array can be EMPTY, in which case an
%                    omni-directional mode is selected
%  * LagMidTol = (nDir x 1) cell array with each cell containing a
%              (nLag(idir) x 2) array [lagMid lagTol] with:
%              lagMid = midpoint of distance class
%              lagTol = tolerance for each class
%              NOTE: Distance classses are thus defined as:
%                    [ lagMid-lagTol, lagMid+lagTol ]
%              NOTE: This could also be a regular array NOT a cell array
%                    In this case, the # of directions nDir is determined
%                    by array DirSpecs and the same distances & tolerances 
%                    are used for ALL directions
%  * MeasSpecs = (nMeas x 3) array with column #s (in DataIn) for attributes
%              to consider as tail & head variables, and structural 
%              measure type to compute FOR EACH direction
%
% For more details, please type |help scatter2struct|.

coords_indic=class2indic(coords,3,[1,2]);
meas2= [3 3 5;3 4 5;3 5 5; 4 3 5;4 4 5;4 5 5;5 3 5;5 4 5; 5 5 5];
coords_tp=scatter2struct(coords_indic,[1,2],[0,100,100],[1:50';ones(50,1)']',meas2);
coords_tp_data=coords_tp.struct{1};
figure
 ip=0;
 nplot=3;
 for j=1:nplot
     for k=1:nplot
               ip = ip + 1;
               subplot(nplot,nplot,ip);
               plot(coords_tp_data(1:50,1,(j-1)*nplot+k),coords_tp_data(1:50,2,(j-1)*nplot+k),'-');               
               xlabel('distance');
               ylabel('transition probability')
               txt=sprintf('Class %d-%d',j,k);
               title(txt);
      end
 end
 
%% Fitting transiogram models to sample data
% Example of using the GeostatsToolbox function |fitcrosstrans| for 
% fitting transiogram models (via iterative weighted least squares) 
% to empirical transiograms.
%
% GeostatsToolbox function |fitcrosstrans|
%
% Inputs:
%
% *  ModPars contains the common structure specs
% * CoeffMat0 contains initial partial sill estimates, assumed to sum to
%    the respective class proportions
% * distin (Lagsize*1) contains lag distance
% * Gammahat (nVar*nvar*nLagsize)contains sample transiogram values corresponding to distin.
% * iWeight weights when computing square of difference
% * optimPars structure to specify the stop condition, [Max Interations, Converge Tolerance], e.g.[1000,0.1]
%  means run 1000 interations or the difference of each interations is
%  smaller than 0.1, whichever happen first.
%
% Outputs:
%
% * coeff: the output improved versions of CoeffMat0
% * ModParsMulti containts all auto- and cross-transiogram models
%   with auto-sills summing to 1-p, and cross-sill summing to p
temp=reshape(sim_tp_data(1:100,2,:),100,9)';
gammahat=reshape(temp,3,3,100);
modpars0=[1,0.9,0,5,5,0];
coeffmat=repmat(classprops,3,1)*eye(3);
[coeff,modmulti]=fitcrosstrans(lagdist100',gammahat,modpars0,coeffmat,1,[10000,0.01]);

%% The need for fusing spatial transition probabilities
% As indicator cross-covariances, transition probabilities are also
% two-point statistics; that is, they quantify spatial structure based
% only on pairs of class labels at different locations or pixels. An 
% efficient method that can integrate such two-point statistics into 
% a multi-point measure of spatial structure, needs to be developed 
% for capturing the complex spatial patterns found in categorical fields.

%% Transition probability fusion: Spatial Markov Chain model
% Based on the concept of transiogram and under the assumption of 
% conditional independence, Li (2007) proposed fusing spatial transition
% probabilities using a single Markov Chain (spatial Markov Chain model)
% moving randomly within a stationary random field, the so called a 
% Markov Chain Random Field (MCRF), where the unknown state at an 
% arbitrary location depends entirely only on the states at its
% nearest neighbours (Markovian property).

%% Transition probability fusion: Permanence of odds ratios
% The conditional indepence assumption provides great flexibility at
% the cost of an over-simplified fusion model. This may lead to incorrect
% conclusions due to the overwhelming dependencies in a spatial context.
% Hence, the consequences of this simplification should be investigated
% carefully [Journel, 2002]. The assumption of permanence of odds ratios
% is another way to fuse transition probabilities into a conditional 
% (multi-point) probability of class occurrence. The idea behind this
% new assumption is that ratios of information increments are typically
% more stable than the increments themselves. Compared to the assumption
% of conditional independence, the permanence of odds ratios assumption
% avoids the calculation of the marginal probability in Bayesian expansion.
% It can be easily demonstrated that the model of permanence of odds 
% ratios is a more general, hence more flexible and less unrealistic, 
% form of conditional independence.

%% Transition probability fusion: Modified permanence of odds ratios
% As stated above, the permanence of ratios approximation also assumes 
% a certain form of independence between class labels at different 
% locations. To relax this assumption, an exponent factor , $\tau_n$, 
% is introduced to account for information redundancy between location 
% pairs. The most challenging aspect of this model is the computation
% of fusion weights, $\tau_n$.

%% Indicator sequential simulation 
% Stochastic simulation is a broadly accepted tool for studying the 
% properties of a statistical model through the generation of alternative
% realizations (simulated attribute images in a spatial context) from
% that model. Indicator sequential simulation is a frequently used 
% simulation method for generating realizations of categorical fields. 
% In what follows, unconditional indicator sequential simulation is 
% performed first, so that the spatial patterns implied by the 
% statistical model are revealed without the influence or possible 
% noise due to the lack of sample data. To illustrate the ability to
% constrain realizations to observed class labels, conditional 
% simulation is then conducted.
%
% In the following examples, the truncated multivariate Gaussian model 
% is taken as reference model, against which fusion, hence approximation, 
% models (Spatial Markov Chain, Permanence of odds ratios and modified 
% permanence of odds ratios) are compared. This reference is selected 
% because the posterior multi-point probability in this case can be 
% computed numerically with great precision and without any independence
% approximation. 

%% GeostatsToolbox function |sisimtrans|
% Function to perform sequential indicator simulation (SISIM) at scattered
% locations or grid nodes of a raster, using various types of transition
% probability fusion.
%
% Inputs:
%
% * CoordsSim = [nx xmin xsiz; ny xmin ysiz]; 
%                nx, xmin, xsiz: # of nodes, minimum x-coordinate,
%                and fixed cell size in x-direction
% * DataIn    = (nDat x nCol) array with sample data (conditioning data)
%               Missing values, flagged as NaNs in any column, are ignored
% * icolC     = (1 x nDim) array with column # (in DataIn) of sample
%               coordinates
% * icolV     = (1 x 1,1+nPrior) array with column #s (in DataIn) for data
% * classLabs = (1 x nClass) array with class labels, e.g., [0 1 2]
% * ModPars   = (nStruct x 6) array with covariogram parameters for
%               a standard Gaussian random function.
% * globProps = (1 x nClass) array with global (regional) class proportions
% * srchSpecs = (1 x 1,3,6) array with search ellipsoid parameters:
% * srchData  = (1 x 3,4) array [ndMin ndMax nsMax (ndOctMax)] with 
%               data parameters for local search:
% * simPars   = (1 x 2,3) array [nSim seed (nMulti)] with simulation 
%               parameters:
% * iweights = different transtion probablity combination methods.
%   
%
% Outputs:
%
% * Out       = output SINGLE precision array with simulated class labels; 
%               each column contains simulated labels for one realization.
%               Out has size (nx*ny*nSim)
%
modpars=[0,0.001,0,5,5,0;1,0.999,0,5,5,0];
classprops=[0.3,0.45,0.25];

%% Unconditional simulation of a truncated multivariate Gaussian field

% set the size of each realization = 100;
nsize=100;
% configuration of modpars
modpars = [0,0.001,0,5,5,0;1,0.999,0,5,5,0];
% class proportion for each class 
classprops=[0.3,0.45,0.25];
% Set iweight=12 to conduct simulation in truncated multivariate Gaussian
% field. Please note that this operation will take a while to finish 
% due to the high dimensional intergral involved.
trunc_gauss=sisimtrans2([nsize,1,1;nsize,1,1],[],1,1,[1 2 3],modpars,classprops,[0,30,30],[50,15],[10,234,0],12);
figure;
rastermap(gridspecs100,trunc_gauss,10,'cat');
 
%% Unconditional simulation using the Spatial Markov Chain model

% Set iweight=0 to conduct unconditional simulation using spatial Markov
% chain model (under the conditional independence assumption)
smc=sisimtrans2([nsize,1,1;nsize,1,1],[],1,1,[1 2 3],modpars,classprops,[0,30,30],[10,15],[50,234,0],0);
figure;
rastermap(gridspecs100,smc,10,'cat');



%% Unconditional simulation using the permanence of odds ratios model

% Set iweight=3 to conduct unconditional simulation using permanence of
% odds ratios model
tau1=sisimtrans2([nsize,1,1;nsize,1,1],[],1,1,[1 2 3],modpars,classprops,[0,30,30],[10,15],[50,234,0],3);
figure;
rastermap(gridspecs100,tau1,10,'cat');


%% Unconditional simulation using the modified permanence of odds ratios model

% Set iweight=6 to conduct unconditional simulation using modified 
% permanence of odds ratios model
tau2=sisimtrans2([nsize,1,1;nsize,1,1],[],1,1,[1 2 3],modpars,classprops,[0,30,30],[10,15],[50,234,0],6);
figure;
rastermap(gridspecs100,tau2,10,'cat');


%% Conditioning data
% Conditioning data that contain 1000 points.
figure;
scattermap(coords,[1,2],10,3,'cat');

%% Conditional simulation using the spatial Markov Chain model

% set iweight=0 to conduct conditional simulation of spatial Markov chain
% model
smc_cond=sisimtrans2([nsize,1,1;nsize,1,1],cond,1,1,[1 2 3],modpars,classprops,[0,30,30],[10,15],[10,234,0],0);
figure;
rastermap(gridspecs100,smc_cond,10,'cat');



%% Conditional simulation using the permanence of odds ratios model

% set iweight=3 to perform conditional simulation of permanence of odds
% ratios model
tau1_cond=sisimtrans2([nsize,1,1;nsize,1,1],cond,1,1,[1 2 3],modpars,classprops,[0,30,30],[10,15],[10,234,0],3);
figure;
rastermap(gridspecs100,tau1_cond,10,'cat');


%% Conditional simulation using the modified permanence of odds ratios model

% set iweight=6 to conduct conditional simulation of modified permanence of
% odds ratios model
tau2_cond=sisimtrans2([nsize,1,1;nsize,1,1],cond,1,1,[1 2 3],modpars,classprops,[0,30,30],[10,15],[10,234,0],6);
figure;
rastermap(gridspecs100,tau2_cond,10,'cat');


%% Transiogram Reproduction
% As mentioned earlier, auto- and cross-class transiogram are used as 
% measures of spatial pattern. To check how well different fusion 
% models reproduce the spatial patterns implied in the reference model 
% (truncated multivariate Gaussian model), we need to check how well 
% the auto- and cross- transiogramw of that reference model are reproduced
% in each model. 
%
% From the results presented hereafter, one can appreciate that the 
% modified permanence of odds ratios model generates the best results
% in terms of transiogram reproduction. The reproduced transiograms 
% (red solid lines) almost overlap with the transiograms of the 
% reference model (green solid lines). The differences between the 
% blue lines and green lines indicate that the spatial Markov chain 
% model does not reproduce the pattern implied in the reference model 
% very well. 
%
% We apply the same procedure for checking transiogram reproduction in
% realizations of conditional simulation and arrive at the same conclusion.

nClass=3
fh=100;
methodnum=4;
clear tpmat;
nsim=50;
out = nan(nsize*nsize,nsim);
out=cell(methodnum,1);
linetype={'g-','b-','c-','r-','b--','c--','r--','m-','y-','k-','r:','g:','b:','c:','m:','y:','k:'};
names={'TGS','SMC','Tau w t=1','Tau w OK'};
out{1}=trunc_gauss;
out{2}=smc;
out{3}=tau1;
out{4}=tau2;
gridspecs=gridspecs100;
meas= [1 1 5;1 2 5;1 3 5; 2 1 5;2 2 5;2 3 5;3 1 5;3 2 5; 3 3 5];
for m = 1:methodnum
   tptotew=zeros(nsize+1,10,9);
   tptotsn=zeros(nsize+1,10,9);
   result =out{m};
   for i = 1:nsim
      progress=sprintf('converting %d class data to indicator data',i);
      disp(progress);
       indic=class2indic(result,i);
       ucwtp=raster2structdir(indic,gridspecs,[1 0;0 1],[nsize,nsize],meas,0);
       tptotew=tptotew+ucwtp.struct{1};
       tptotsn=tptotsn+ucwtp.struct{2};
   end
   tptotew=tptotew/nsim;
   tptotsn=tptotsn/nsim;
   figure(fh);
   ip = 0;
   for i = 1:nClass
       for j = 1:nClass
           ip = ip + 1;
           subplot(nClass,nClass,ip);
           hold on;
           plot(tptotew(:,1,(i-1)*nClass+j),tptotew(:,2,(i-1)*nClass+j),linetype{m},'LineWidth',1);
           xlim([0 nsize]);   
           xlabel('distance');
           ylabel('transition probability');
           txt=sprintf('Class %d-%d',i,j);
           title(txt);
       end
   end
end
hold off;


%% References
% Carle, S. F. and Fogg, G. E. (1996), ‘Transition probability-based indicator geostatistics’,
% Mathematical Geology 28(4), 453–476.
%
% Carle, S. F. and Fogg, G. E. (1997), ‘Modeling spatial variability with one and multidi-
% mensional continuous-lag Markov chains’, Mathematical Geology 29(7),
% 891–918.
%
% Li, W. (2007), ‘Markov Chain random fields for estimation of categorical variables’,
% Mathematical Geology 39, 321–335.
%
% Journel, A. (2002), ‘Combining knowledge from diverse sources: an alternative to traditional data independence hypothesis’, Mathematical Geology 34(5),
% 573–596.