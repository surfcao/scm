function Out=dist2msemivar(D,sizeD,ModParsIso,nStruct)
% Isotropic model semivariogram values for given distances
%
%% DESCRIPTION: dist2msemivar.m
% Function to compute semivariogram values at specified distances (in D),
% given a theoretical model specification (in ModParsIso).
% NO ERROR CHECKING IS PERFORMED, since this function is called by others.
% Similar to dist2modelstruct, but with no error checking and no imeas
%
%% SYNTAX:
%   Out = dist2msemivar(D,sizeD,ModParsIso,nStruct);
%
%% INPUTS:
%   D            = input array of distance values
%   sizeD        = (1 x 2) array with size of distance matrix D
%   ModParsIso   = array with isotropic model parameters
%                  e.g., for 2 nested structures
%                  [mtype(1) msill(1) mrange(1) extraPar(1); 
%                   mtype(2) msill(2) mrange(2) extraPar(2)];
%                   mtype = functional type:
%                   (0=nugget, 1=spherical, 2=exponential,
%                    3=Gaussian, 4=power, 5=cosine, 6=dampened cosine,
%                    7=dampened sine, 8=logarithmic, 9=Whittle, 10=cubic,
%                    11=generalized Cauchy, 12=stable)
%                   extraPar is a non-negative free parameter.
%                   For a nugget structure, mrange is a free parameter
%                   For power models, msill = free, mrange = slope, 
%                                     extraPar = power
%                   For logarithmic models extraPar = microstructure
%                   For cosine, dampened sine cosine models,
%                   mrange = wavelength ~ 2*size of underlying objects
%                   For dampened cosine model: extraPar = length at which
%                   95% of hole-effect is dampened out
%   nStruct      = scalar with # of nested structures
%
%% OUTPUTS:
%   Out          = output array, of same size as D, with model-derived
%                  semivariogram values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%   Out = dist2msemivar(D,sizeD,ModParsIso,nStruct);

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
TINY     = 1/10^10;

%% Get some input parameters
mtype    = ModParsIso(:,1);
msill    = ModParsIso(:,2);
mrange   = ModParsIso(:,3);
extraPar = ModParsIso(:,4);

%% Declare output
Out = zeros(sizeD);

%% Loop over all structures, and compute covariance
for ist=1:nStruct
   switch mtype(ist);
      case 0 %Nugget model
         Out = Out + msill(ist)*(D ~= 0);
      case 1 %Spherical model
         Dr = min(D./mrange(ist),1);
         Out = Out + msill(ist)*( 1.5*Dr-0.5*Dr.^3 );
      case 2 %Exponential model
         Out = Out + msill(ist)*( 1- exp(-3*D/mrange(ist)) );
      case 3 %Gaussian model
         Out = Out + msill(ist)*( 1 - exp ( -3*(D/mrange(ist)).^2 ) );
      case 4 %Power model
         Out = Out + mrange(ist)*(D.^extraPar(ist));
      case 5 %Sine model
         Out = Out + msill(ist)*( 1 - cos(2*pi*D/mrange(ist)) );
      case 6 %Dampened cosine model
         Out = Out + msill(ist)*...
             ( 1 - exp(-3*D./extraPar(ist)).*cos(2*pi*D/mrange(ist)) );
      case 7 %Cardinal sine model
         Dr = 2*pi*D/mrange(ist);
         OutTmp = msill(ist)*( 1 - ( sin(Dr)./max(TINY,Dr) ));
         OutTmp( Dr <= TINY ) = 0;
         Out    = Out + OutTmp;
      case 8 %Logarithmic model
          X = sqrt(D.^2 + extraPar(ist)^2)./mrange(ist);
          Out = Out + msill(ist)*(log(X));
      case 9 %Whittle's correlation
         Dr = D/mrange(ist);
         OutTmp = msill(ist).*( 1 - Dr.*besselk(1,Dr) );
         OutTmp(isnan(OutTmp)) = 0;
         Out = Out + OutTmp;
      case 10 %Cubic model
         Dr = min(D./mrange(ist),1);
         Out = Out + msill(ist)*( 1 - (1 - 7*Dr.^2 + 35/4*Dr.^3 ...
                                    - 7/2*Dr.^5 + 3/4*Dr.^7 ));
      case 11 %Generalized Cauchy model
         Out = Out + msill(ist)*...
             ( 1 - (1 + (D/mrange(ist)).^2) .^(-extraPar(ist)) );
      case 12 %Stable model
         Out = Out + msill(ist)*...
             ( 1 - exp ( -3*(D/mrange(ist)).^extraPar(ist) ) );
      otherwise
         error('Invalid structure type');
   end
end