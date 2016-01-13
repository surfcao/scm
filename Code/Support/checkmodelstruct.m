function checkmodelstruct(ModPars,nDim)
% Check consistency of parameters in models of spatial structure
%
%% DESCRIPTION: checkmodelstruct.m
% Function to check the consistency of a parameter array ModPars
% specifying a theoretical model of spatial structure.
% NOTE: Only consistency between type and extraPar with respect to
%       dimensionality is checked.
%
%% SYNTAX:
%   checkmodelstruct(ModPars,nDim);
%
%% INPUTS:
%   ModPars   = array with isotropic model parameters
%               e.g., for 2 nested structures
%                  [mType(1) mSill(1) mRange(1) extraPar(1)]; 
%                   mType(2) mSill(2) mRange(2) extraPar(2)];
%                   mType = functional type:
%                   (0=nugget, 1=spherical, 2=exponential,
%                    3=Gaussian, 4=power, 5=cosine, 6=dampened cosine,
%                    7=dampened sine, 8=logarithmic, 9=Whittle, 10=cubic,
%                    11=generalized Cauchy, 12=stable)
%   nDim      = scalar with dimensionality
%               NOTE: Dimensionality is determined by nDim, not by the
%                     # of columns in ModPars
%
%% OUTPUTS:
%   None. If an inconsistency is found, an error message is reported,
%         and the program is terminated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%   checkmodelstruct(ModPars,nDim);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2005                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error checking
if ~isscalar(nDim)
    error('nDim must be a scalar');
end
if ~(nDim == 1 || nDim == 2 || nDim == 3)
    error('Dimensionality must be 1, 2 or 3');
end
%% Get some input parameters
[nStruct,nCols] = size(ModPars);
switch nDim
    case 1
        if nCols ~= 4
            error('In 1D, you need 4 columns in array ModPars');
        end
    case 2
        if nCols ~= 6
            error('In 2D, you need 6 columns in array ModPars');
        end
    case 3
        if nCols ~= 9
            error('In 3D, you need 9 columns in array ModPars');
        end
end
mType = ModPars(:,1);
%mSill = ModPars(:,2);
extraPar = ModPars(:,end);

%% Check sills and ranges
TINY = 1/10^10;
%if any( mSill < TINY )
%    error('Sill parameters must be positive');
%end
switch nDim
    case 1
        mRange = ModPars(:,3);
    case 2
        mRange = ModPars(:,[4 5]);
    case 3
        mRange = ModPars(:,[6 7 8]);
end
if any( mRange(mType~=0) < TINY ) 
    error('Range parameter(s) must be positive');
end

%% Loop over all structures, and check extra parameter
for ist=1:nStruct
   switch mType(ist);
      case 0 %Nugget model
      case 1 %Spherical model
      case 2 %Exponential model
      case 3 %Gaussian model
      case 4 %Power model
         if extraPar(ist) <=0 || extraPar(ist) >= 2
            error('For power model (4), extraPar must lie in (0,2)');
         end
      case 5 %Cosine model
      case 6 %Dampened cosine model
         if extraPar(ist) <= 0
         error('For dampened cosine (6), dampening parameter must be > 0'); 
         end
         switch nDim
             case 1
                 if extraPar(ist) < 0
               error('For 1D dampened cosine (6), extraPar must be >= 0');
                 end
             case 2
                 if extraPar(ist) < 1
               error('For 2D dampened cosine (6), extraPar must be >= 1');
                 end
             case 3
                 if extraPar(ist) < sqrt(3)
          error('For 3D dampened cosine (6), extraPar must be >= sqrt(3)');
                 end
         end
      case 7 %Dampened sine model
      case 8 %Logarithmic model
      case 9 %Whittle's model
      case 10 %Cubic model
      case 11 %Generalized Cauchy model
         if extraPar(ist) <= 0
            error('For Cauchy model (11), extraPar must be > 0');
         end
      case 12 %Stable model
         if extraPar(ist) > 2 || extraPar(ist) < 0
            error('For stable model (12), extraPar must lie in [0,2]');
         end
      otherwise
         error('Invalid structure type');
   end
end

%% Finished
%disp(' ');
%disp('Passed structural model specification check');