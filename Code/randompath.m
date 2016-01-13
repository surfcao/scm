function out = randompath(iGrid,ngIn,seedIn,nMulti)
%  Indices in a random sequence (or within nested grids)
%
%% DESCRIPTION: randompath.m
% Function to compute a random path for visiting nodes in sequential 
% simulation. No seed is provided, because it is set outside this function,
% in the main sequential simulation function.
%
%% SYNTAX:
%   out = randompath(iGrid,ngIn,seedIn,nMulti);
%
%% INPUTS:
%   iGrid     = scalar with flag for grid mode
%                   iGrid  = 0 -> scattered locations
%                   iGrid ~= 0 -> locations at grid nodes
%   ngIn      = if iGrid  = 0 -> ngIn = nLoc = scalar with # of locations
%               if iGrid ~= 0 -> ngIn = ng = (1 x nDim) array with # of 
%                                grid nodes per direction, e.g., [ngx ngy]
%               NOTE: Dimensionality is determined by this array
%   seedIn    = scalar with random # seed
%   nMulti    = scalar with # of multiple nested grid to consider
%               If nMulti =  0, output indices are completely random
%               If nMulti ~= 0, the output indices are random within 
%                               these nMulti nested grids
%
%% OUTPUTS:
%   out       = (ngIn x 1) or (prod(ngIn) x 1) array (single precision) 
%               with indices of grid nodes defining a visiting sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% SYNTAX:
%   out = randompath(iGrid,ngIn,seedIn,nMulti);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                              May 2004                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error checking
if nargin ~= 4
    error('You need 4 input arguments for this function');
end
if ~isscalar(seedIn)
    error('seed must be a scalar');
end
if iGrid == 0 && nMulti > 0
    error('For scattered locations, iGrid = 0, nMulti must be 0'); 
end
if ~isscalar(nMulti)
    error('nMulti must be a scalar');
end
if nMulti < 0
    error('nMulti must be >= 0');
end
if nMulti > 0 && rem(nMulti,fix(nMulti)) ~= 0
    error('nMulti must not have decimals');
end

%% Get some parameters
rand('state',seedIn);
if iGrid ~= 0
    ng = ngIn(:)';
    nLoc = prod(ng);
    nDim = numel(ng);
else
    nLoc = ngIn;
end

%% Initialize indices with uniform random numbers
out = rand(single(nLoc),single(1));

%% Stratified random path or multiple grid search
if iGrid ~= 0 && nMulti > 0
    s1 = single(1);
    switch nDim
    case 1
        for im=1:nMulti
            im4 = im*4;
            nnx = max( 1 , ng(1)/im4 );
            ix = (s1:nnx)';
            jx = ones(size(ix),'single');
            jx(ix>1) = ix(ix>1)*im4;
            index = jx;
            out(index) = out(index) - im;
        end
    case 2
        for im=1:nMulti
            im4 = im*4;
            nnx = max( 1 , ng(1)/im4 );
            nny = max( 1 , ng(2)/im4 );
            [ix,iy] = ndgrid(s1:nnx,s1:nny);
            jx = ones(size(ix),'single');
            jy = ones(size(iy),'single');
            jx(ix>1) = ix(ix>1)*im4;
            jy(iy>1) = iy(iy>1)*im4;
            index = jx(:) + (jy(:)-1)*ng(1);
            out(index) = out(index) - im;
        end
    case 3
        ngXY = ng(1)*ng(2);
        for im=s1:nMulti
            im4 = im*4;
            nnx = max( 1 , ng(1)/im4 );
            nny = max( 1 , ng(2)/im4 );
            nnz = max( 1 , ng(3)/im4 );
            [ix,iy,iz] = ndgrid(s1:nnx,s1:nny,s1:nnz);
            jx = ones(nnx,nny,nnz,'single');
            jy = ones(nnx,nny,nnz,'single');
            jz = ones(nnx,nny,nnz,'single');
            jx(ix>1) = ix(ix>1)*im4;
            jy(iy>1) = iy(iy>1)*im4;
            jz(iz>1) = iz(iz>1)*im4;
            index = jx(:) + (jy(:)-1)*ng(1) + (jz(:)-1)*ngXY;
            out(index) = out(index) - im;
        end
    end
end


%% For-loop version of above
%{
if iGrid ~= 0 && nMulti > 0
    ng3 = [1 1 1];
    ng3(1:nDim) = ng;
    for im=1:nMulti
        im4 = im*4;
        nnz = max( 1 , ng3(3)/im4 );
        nny = max( 1 , ng3(2)/im4 );
        nnx = max( 1 , ng3(1)/im4 );
        jz = 1; jy = 1; jx = 1;
        for iz=1:nnz
            if nnz > 1; jz = iz*im4; end;
            for iy=1:nny
                if nny > 1; jy = iy*im4; end;
                for ix=1:nnx
                    if nnx > 1; jx = ix*im4; end;
                    index = jx + (jy-1)*ng3(1) + (jz-1)*ngXY;
                    out(index) = out(index) - im;
                end
            end
        end
    end
end
%}

%% Sort uniform random deviates = randperm
[dump,out] = sort(out);
clear dump;
out = single(out);

%% Finished
disp(' ');
if iGrid ~= 0 && nMulti > 0
    tit = ['Computed random path over ',num2str(nMulti),' nested grids'];
    disp(tit);
else
    disp('Computed random path');
end