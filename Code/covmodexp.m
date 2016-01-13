function gam = covmodexp(modPars,lagDist)
% Exponential model semivariogram values for given distances (isotropic)
%
%% DESCRIPTION: vmodexp.m
% Function to compute semivariogram values at given distances in lagDist,
% from an isotropic EXPONENTIAL model with parameters in modPars.
% NOTE: The 1st entry in modPars pertains to the nugget effect.
%
%% SYNTAX:
%   gam = vmodexp(modPars,lagDist);
%
%% INPUTS:
%   modPars      = (1 x 3) array [nugg sill range] with isotropic 
%                  semivariogram model parameters:
%                       nugg  = partial sill for nugget effect
%                       sill  = partial sill for exponential model
%                       range = practical range for exponential model  
%   lagDist      = (nLag x 1) array of distance values
%
%% OUTPUTS:
%   gam          = (nLag x 1) array with model semivariogram values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%   gam = vmodexp(modPars,lagDist);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2007                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nugg = modPars(1); sill = modPars(2); range = modPars(3);
%% Error checking
if numel(modPars) ~= 3
    error('modPars must have 3 entries');
end

%% Exponential covariogram model
gam = exp(-3.*lagDist./modPars(3));
gam = modPars(1) + modPars(2).*gam;
% Nugget contribution applies only to non-zero distances
gam(lagDist<eps) = 0;