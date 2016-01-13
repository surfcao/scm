function Out = coords2sqdist(Coords1,Coords2)
% Squared Euclidean distance array from 1 or 2 sets of coordinates
%
%% DESCRIPTION: coords2sqdist.m
% Function to compute squared Euclidean distance between all pairs of 
% points in 1 or 2 sets of locations with given coordinates.
% No NaNs are allowed, since this function is internal to many other ones.
% NO ERROR CHECKING IS PERFORMED.
%
%% SYNTAX:
%    Out = coords2sqdist(Coords1);
%    Out = coords2sqdist(Coords1,Coords2);
%
%% INPUTS:
%    Coords1    = (nD1 x nDim) array with coordinates set #1
%    Coords2    = Optional (nD2 x nDim) array with coordinates set #2
%
%% OUTPUTS: 
%    Out        = (nD1 x nD1) or (nD1 x nD2) array of Euclidean distances
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%    Out = coords2sqdist(Coords1);
%    Out = coords2sqdist(Coords1,Coords2);

%% CREDITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2005                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get some parameters
[nD1,nDim] = size(Coords1);

%% Single set of locations
if nargin == 1
    Out = zeros(nD1,nD1);
    % Loop over # of locations in set #1
    for ii=1:nD1
        % Loop over # of dimensions & sum up squared vector components
        distSq = zeros(nD1,1);
        for id = 1:nDim
            distSq = distSq + (Coords1(ii,id) - Coords1(:,id)).^2;
        end
        % Fill-in current portion of output array
        Out(:,ii) = distSq;
    end
    return
end

%% Loop over set of coordinates with fewer entries
nD2 = size(Coords2,1);
Out = zeros(nD1,nD2);
if nD1 <= nD2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop over all locations of set #1 (row-wise loop)
    for ii=1:nD1
        % Loop over # of dimensions & sum up squared vector components
        distSq = zeros(nD2,1);
        for idim=1:nDim
            distSq = distSq + (Coords1(ii,idim) - Coords2(:,idim)).^2;
        end
        % Fill-in current portion of output array
        Out(ii,:) = distSq;
    end
    
else  %%% nD1 > nD2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop over all locations of set #2 (column-wise loop)
    for jj=1:nD2
        % Loop over # of dimensions & sum up squared vector components
        distSq = zeros(nD1,1);
        for idim=1:nDim
            % Compute vector component for this dimension
            distSq = distSq + (Coords1(:,idim) - Coords2(jj,idim)).^2;
        end
        % Fill-in current portion of output array
        Out(:,jj) = distSq;
    end
    
end          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end switch nD1 <= nD2