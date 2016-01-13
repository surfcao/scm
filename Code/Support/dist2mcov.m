function C = dist2mcov(D,sizeD,CovParsIso,nStruct)
% Isotropic model covariance values for given distances
%
%% DESCRIPTION: dist2mcov.m
% Function to compute covariance values at specified distances (in D),
% given a theoretical covariance model specification (in CovParsIso).
% Addititive nested structures are handled via multiple rows in CovParsIso
% NOTE: No error checking is performed.
%
%% SYNTAX:
%   C = dist2mcov(D,sizeD,CovParsIso,nStruct);
%
%% INPUTS:
%   D          = input array of distance values
%   sizeD      = (1 x 2) array with size of distance matrix D
%   CovParsIso = array with isotropic covariance model parameters
%                e.g., for 2 nested structures
%                [ctype(1) csill(1) crange(1) extraPar(1); 
%                 ctype(2) csill(2) crange(2) extraPar(2)];
%                ctype = functional type:
%                (0=nugget, 1=spherical, 2=exponential,
%                 3=Gaussian, 4=power, 5=cosine, 6=dampened cosine,
%                 7=dampened sine, 8=logarithmic, 9=Whittle, 10=cubic,
%                 11=generalized Cauchy, 12=stable)
%                extraPar is a non-negative free parameter.
%                For a nugget structure, crange is a free parameter
%                For power models, csill = free parameter, crange = slope,
%                                  extraPar = power
%                For logarithmic models extraPar = microstructure
%                For cosine, dampened sine cosine models,
%                crange = wavelength ~ 2*size of underlying objects
%                For dampened cosine model: extraPar = length at which
%                95% of hole-effect is dampened out
%   nStruct    = scalar with # of nested structures
%
%% OUTPUTS:
%   C          = output covariance array, of same size as D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%   C = dist2mcov(D,sizeD,CovParsIso,nStruct);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2005                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get some parameters
TINY     = 1/10^10;
ctype    = CovParsIso(:,1);
csill    = CovParsIso(:,2);
crange   = CovParsIso(:,3);

%% Declare output
C = zeros(sizeD);

%% Loop over all structures, and compute covariance
for ist=1:nStruct
   switch ctype(ist);
      case 0 %nugget model
         %C = C + csill(ist)*(D<=TINY);
         C = C + csill(ist)*(D==0);
      case 1 %Spherical model
         r = max(crange(ist),TINY);
         C = C + csill(ist)*(1.0 - (1.5*min(D/r,1)-0.5*min(D/r,1).^3) );
      case 2 %Exponential model
         r = max(crange(ist),TINY);
         C = C + csill(ist)*( exp(-3*D/r) );
      case 3 %Gaussian model
         r = max(crange(ist),TINY);
         C = C + csill(ist)*( exp ( -3*(D/r).^2 ) );
      case 4 %Power model
         pow = CovParsIso(ist,end);
         %if pow <=0.0 || pow >= 2.0
         %   error('For power model, 0 < extraPar (=power) < 2');
         %end
         %BIGNUM   = 10^10;
         %C = C + BIGNUM - csill(ist)*(D.^r);
         C = C + 10^10 - crange(ist)*(D.^pow);
         % Allows anisotropic slopes (ranges) for given power (extraPar)
      case 5 %Cosine model
         r = max(crange(ist),TINY);
         C = C + csill(ist)*( cos(2*pi*D/r) );
      case 6 %Dampened cosine model
         extraPar = CovParsIso(ist,4);
         %if extraPar < TINY
         %   error('Invalid dampening parameter'); 
         %end
         r = max(crange(ist),TINY);
         C = C + csill(ist)*( exp(-3.0*D./extraPar) .* cos(2*pi*D/r) );
      case 7 %Cardinal sine model
         r = max(crange(ist),TINY);
         Dtmp = 2*pi*D/r;
         Ctmp = csill(ist)*( sin(Dtmp)./max(TINY,Dtmp) );
         %Ctmp = csill(ist)*( sin(pi*D/r) );
         Ctmp(Dtmp <= TINY) = csill(ist);
         C    = C + Ctmp;
      case 8 %Logarithmic model
         extraPar = CovParsIso(ist,4);
         r = max(crange(ist),TINY);
         X = sqrt(D.^2 + extraPar^2)./r;
         C = C + 10^10 - csill(ist)*(log(X));
      case 9 %Whittle's correlation
         r = max(crange(ist),TINY);
         Ctmp = csill(ist).*(D/r).*besselk(1,D/r);
         inan = find(isnan(Ctmp));
         if isempty(inan) == 0
            Ctmp(inan) = csill(ist);
         end
         C = C + Ctmp;
      case 10 %Cubic model
         r = max(crange(ist),TINY);
         Dt = min(D/r,1);
         C = C + csill(ist)*(1.0 - 7*Dt.^2 + 35/4*Dt.^3 ...
                                 - 7/2*Dt.^5 + 3/4*Dt.^7  );
      case 11 %Generalized Cauchy model
         extraPar = CovParsIso(ist,4);
         %if extraPar <= 0
         %   error('For Cauchy model (11), extraPar must be > 0');
         %end
         r = max(crange(ist),TINY);
         C = C + csill(ist)*( 1 + (D/r).^2 ).^(-extraPar);
      case 12 %stable model
         extraPar = CovParsIso(ist,4);
         %if extraPar > 2
         %   error('For stable model (12), extraPar must be <= 2');
         %end
         r = max(crange(ist),TINY);
         C = C + csill(ist)*( exp ( -3*(D/r).^extraPar ) );
   end
end