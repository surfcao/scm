function Out = geoeas2matlab(dataIn,nxnynz)
% Transform GeoEAS array into Matlab array
%
%% DESCRIPTION: geoeas2matlab.m
% Function to transform a SINGLE GeoEAS-formatted raster (dataIn)
% i.e., a single column, to a MATLAB array that can be viewed using 
% imagesc (in 2D) or slice (in 3D).
% NOTE: No error checking is performed
% NOTE: Avoid using unecessary 1s in nxnynz
%
%% SYNTAX:
%   Out = geoeas2matlab(dataIn,nxnynz);
% 
%% INPUTS:
%   dataIn    = input GeoEAS-formatted raster of dimensions:
%               in 1D: (nx x 1)
%               in 2D: (nx*ny x 1)
%               in 3D: (nx*ny*nz x 1)
%   nxnynz    = (1 x nDim) vector with raster size
%               in 1D: [nx]
%               in 2D: [nx ny] 
%               in 3D: [nx ny nz]
%
%% OUTPUTS:
%   Out       = output MATLAB array of dimensions:
%               in 1D: (nx x 1)
%               in 2D: (ny x nx) 
%               in 3D: (ny x nx x nz)
%
%% NOTES: 
%   In 3D, z increases downwards in Matlab's image coordinate system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%   Out = geoeas2matlab(dataIn,nxnynz);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2005                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some error checking
%nDim = length(nxnynz);
%if size(dataIn,2) ~= 1
%   error('dataIn must be a column vector');
%end
%if length(dataIn) ~= prod(nxnynz(:));
%   error('Length of dataIn incompatible with nxnynz');
%end
%if any(nxnynz == 1)
%    error('Found a 1 in nxnynz');
%end

%% Proceed according to dimensionality
nDim = sum(nxnynz~=1);
switch nDim
case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D case

Out = dataIn;

case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D case

Out = reshape(dataIn,nxnynz(1),nxnynz(2));
Out = permute(Out, [2 1]);
Out = Out(nxnynz(2):-1:1,:);

case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D case

Out = reshape(dataIn,nxnynz(1),nxnynz(2),nxnynz(3));
Out = permute(Out, [2 1 3]);
Out = Out(nxnynz(2):-1:1,:,:);

end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch nDim