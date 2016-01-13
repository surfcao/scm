function Out = indgeoeas2ixiyiz(indexIn,ngIn,nDim,nNodes)
% Cartesian offsets corresponding to indices from GeoEas array
%
%% DESCRIPTION: indgeoeas2ixiyiz.m
% Function to compute the Cartesian offsets (ix/iy/iz) of a set of
% grid node indices from a GeoEas-formatted grid (up to 3D).
% No NaNs allowed. No error checking is performed.
%
%% SYNTAX:
%   Out = indgeoeas2ixiyiz(indexIn,ngIn);
%   Out = indgeoeas2ixiyiz(indexIn,ngIn,nDim,nNodes);
%
%% INPUTS:
%   indexIn   = (nNodes x 1) array with indices of grid nodes from a 
%               GSLIB-formatted Cartesian raster
%   ngIn      = (1 x nDim) array with # of grid nodes of GeoEAS raster
%               along each direction, e.g., in 2D [nx ny]
%   nDim      = Optional scalar with dimensionality
%   nNodes    = Optional scalar with # of elements in indexIn
%
%% OUTPUTS:
%   Out       = (nNodes x nDim) array with Cartesian (not image) offsets 
%               (ix/iy/iz) for the indices in input vector indexIn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%    Out = indgeoeas2ixiyiz(indexIn,ngIn);
%    Out = indgeoeas2ixiyiz(indexIn,ngIn,nDim,nNodes);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2005                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get some input parameters, and declare output array
if nargin == 2
    nDim   = sum(ngIn~=1);
    nNodes = numel(indexIn);
end
Out = zeros(nNodes,nDim,'single');

%% Proceed according to dimensionality
switch nDim
    case 1
        % Compute actual Cartesian coordinates
        Out = indexIn;
    case 2
        % Compute Cartesian ix/iy offsets
        Out(:,2) = fix( (indexIn - 1)./ngIn(1) ) + 1;
        Out(:,1) = indexIn - (Out(:,2)-1).*ngIn(1);
    case 3
        nxy = ngIn(1)*ngIn(2);
        % Compute Cartesian ix/iy/iz offsets
        Out(:,3) = fix( (indexIn-1)./nxy ) + 1;
        Out(:,2) = fix( ( indexIn - (Out(:,3)-1).*nxy-1 )./ngIn(1) ) + 1;
        Out(:,1) = indexIn-(Out(:,3)-1).*nxy-(Out(:,2)-1).*ngIn(1);
end