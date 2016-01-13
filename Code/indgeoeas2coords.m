function Out = indgeoeas2coords(indexIn,GridSpecs,nDim,nNodes)
% Coordinates of indices from GeoEas array
%
%% DESCRIPTION: indgeoeas2coords.m
% Function to compute the Cartesian coordinates of a set of
% grid node indices from a GeoEas-formatted grid (up to 3D).
% NOTE: No NaNs allowed. 
%       No error checking is performed.
%
%% SYNTAX:
%   Out = indgeoeas2coords(indexIn,GridSpecs);
%   Out = indgeoeas2coords(indexIn,GridSpecs,nDim,nNodes);
%
%% INPUTS:
%   indexIn   = (nNodes x 1) array with indices of grid nodes from a 
%               GeoEAS-formatted Cartesian raster
%   GridSpecs = (nDim x 3) matrix with grid specifications
%               e.g., in 2D: [nx xmin xsiz ;  ny ymin ysiz]
%               nx   = # of gridnodes in x-direction
%               xmin = x-coordinate of grid origin (lower left corner)
%               xsiz = cell size along x-direction
%   nDim      = Optional scalar with dimensionality
%   nNodes    = Optional scalar with # of elements in indexIn
%
%% OUTPUTS:
%   Out       = (nNodes x nDim) array with Cartesian (not image) 
%               coordinates for the indices in input vector indexIn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%   Out = indgeoeas2coords(indexIn,GridSpecs);
%   Out = indgeoeas2coords(indexIn,GridSpecs,nDim,nNodes);

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
    nDim = size(GridSpecs,1);
    nNodes = numel(indexIn);
end
Out = zeros(nNodes,nDim);

%% Proceed according to dimensionality
switch nDim
    case 1
        % Compute actual Cartesian coordinates
        Out = ( indexIn - 1 )*GridSpecs(1,3) + GridSpecs(1,2);
    case 2
        % Compute Cartesian ix/iy offsets
        Out(:,2) = fix( (indexIn-1)./GridSpecs(1,1) ) + 1;
        Out(:,1) = indexIn - (Out(:,2)-1).*GridSpecs(1,1);
        % Compute actual Cartesian coordinates
        Out(:,1) = ( Out(:,1) - 1 )*GridSpecs(1,3) + GridSpecs(1,2);
        Out(:,2) = ( Out(:,2) - 1 )*GridSpecs(2,3) + GridSpecs(2,2);
    case 3
        nxy = GridSpecs(1,1)*GridSpecs(2,1);
        % Compute Cartesian ix/iy/iz offsets
        Out(:,3) = fix( (indexIn-1)./nxy ) + 1;
        Out(:,2) = fix( ( indexIn - (Out(:,3)-1).*nxy-1 ) ...
                            ./GridSpecs(1,1) ) + 1;
        Out(:,1) = indexIn-(Out(:,3)-1).*nxy-(Out(:,2)-1).*GridSpecs(1,1);
        % Compute actual Cartesian coordinates
        Out(:,1) = ( Out(:,1) - 1 )*GridSpecs(1,3) + GridSpecs(1,2);
        Out(:,2) = ( Out(:,2) - 1 )*GridSpecs(2,3) + GridSpecs(2,2);
        Out(:,3) = ( Out(:,3) - 1 )*GridSpecs(3,3) + GridSpecs(3,2);
end