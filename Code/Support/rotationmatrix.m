function [RotMat,iAnis] = rotationmatrix(AnisSpecs)
% Anisotropy rotation matrix for vector component rotation
%
%% DESCRIPTION: rotationmatrix.m
% Function to compute matrices of rotation and shrinkage (up to 3D), 
% given anisotropy specifications in array AnisSpecs.
% The dimensionality of the problem (nDim) is determined by
% the # of columns of AnisSpecs: =1->1D, =3->2D, =6->3D.
% The number of structures is determined by the # of rows of AnisSpecs.
%
%% SYNTAX:
%    [RotMat,iAnis] = rotationmatrix(AnisSpecs);
%
%% INPUTS:
%    AnisSpecs  = (nStruct x 1,3,6) array with anisotropy specifications
%                 for each nested structure.
%                 2D example with 2 nested structures: 
%                 [ ang1(1)  rng1(1) rng2(1); ang1(2) rng1(2) rng2(2) ] 
%                 3D example with 1 nested structure: 
%                 [ ang1  ang2 ang3 rng1 rng2 rng3 ]
%                 ALL INPUT ANGLES ARE DEFINED AS IN GSLIB:
%                 ang1: positive degrees clockwise from North.
%                 ang2: negative degrees down from horizontal.
%                 ang3: third rotation angle (plunge).
%                 rng1/rng2/rng3 = radii of anisotropy ellipsoid
%                 along directions specified by ang1/ang2/ang3.
%
%% OUTPUTS:
%    RotMat     = (nDim x nDim x nStruct) array with rotation/shrinkage 
%                 matrix for each nested structure.
%    iAnis      = (1 x nStruct) binary array with flags for 
%                 anisotropy (=1) or isotropy (=0)
%                 NOTE: iAnis(ist) = 1 if all ranges are the same
%                       for structure ist
%
%% NOTES: 
%    All rotation input angles are in degrees (GSLIB style)
%    but they are converted to angles increasing anti-clockwise
%    when observer stands at positive side of rotation axis
%    and is looking towards the origin (e.g., ang1 is positive 
%    anti-clockwise when looking from positive Z-axis towards origin)
%    ang1->  rotation angle of XY-plane around Z-axis (rads)
%    ang2->  rotation angle of XZ-plane around Y-axis (rads)
%    ang3->  rotation angle of YZ-plane around X-axis (rads)
%    anis1-> 1st anisotropy ratio: rng2/rng1
%    anis2-> 2nd anisotropy ratio: rng3/rng1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%    [RotMat,iAnis] = rotationmatrix(AnisSpecs);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                               May 2005                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fixed parameters
TINY = 1/10^20;

%% Determine dimensionality
[nStruct,nCols] = size(AnisSpecs);
%if ~(nCols == 1 | nCols == 3 | nCols == 6)
%    error('Array AnisSpecs must have 1, 3, or 6 columns');
%end

%% Proceed according to dimensionality
switch nCols
 case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D case
 RotMat = ones(1,1,nStruct);          % nothing to rotate in 1D
 iAnis = zeros(1,nStruct);

 case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D case
 RotMat = zeros(2,2,nStruct);
 iAnis = ones(1,nStruct);
 for ist = 1:nStruct
    % convert input angles to mathematical angles
    ang1  = (90-AnisSpecs(ist,1)).*(pi/180);
    if AnisSpecs(ist,2) == AnisSpecs(ist,3)
        iAnis(ist) = 0;
    end
    %rng1  = max(AnisSpecs(ist,2),TINY);
    %if rng1 < TINY; error('Found entry in AnisSpecs(:,2)=rng1 < 0'); end;
    %rng2  = AnisSpecs(ist,3);
    anis1 = AnisSpecs(ist,3)./max(AnisSpecs(ist,2),TINY);
    % Establish rotation matrix
    RotMatTmp = [ cos(ang1) sin(ang1) ; ...
                 -sin(ang1) cos(ang1) ];
    % Compute affine transformation matrix
    f = 1 / max(anis1,TINY);
    F = diag([1 f]);
    % Apply rotation, and then affine transform
    RotMat(:,:,ist) = F*RotMatTmp;
    % Match GSLIB
    RotMat(:,:,ist) = RotMat(:,:,ist)';
 end

 case 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3D case
 RotMat = zeros(3,3,nStruct);
 iAnis = ones(1,nStruct);
 for ist = 1:nStruct
    % convert input angles to mathematical angles
    ang1  = (90-AnisSpecs(ist,1)).*(pi/180);
    ang2  = -AnisSpecs(ist,2).*(pi/180);
    ang3  = AnisSpecs(ist,3).*(pi/180);
    if AnisSpecs(ist,4) == AnisSpecs(ist,5) ...
            && AnisSpecs(ist,4) == AnisSpecs(ist,6)
        iAnis(ist) = 0;
    end
    rng1  = max(AnisSpecs(ist,4),TINY);
    %if rng1 < TINY; error('Found entry in AnisSpecs(:,4)=rng1 < 0'); end;
    %rng2  = AnisSpecs(ist,5); rng3  = AnisSpecs(ist,6);
    anis1 = AnisSpecs(ist,5)./rng1;
    anis2 = AnisSpecs(ist,6)./rng1;
    % Establish rotation matrices
    RotMat1 = [ cos(ang1) sin(ang1)  0 ; ...
               -sin(ang1) cos(ang1)  0 ; ...
                  0         0        1];
    RotMat2 = [ cos(ang2) 0 -sin(ang2) ; ...
                  0    1    0       ; ...
                sin(ang2) 0  cos(ang2)];
    RotMat3 = [ 1      0         0     ; ...
                0  cos(ang3) sin(ang3) ; ...
                0 -sin(ang3) cos(ang3) ];
    % Compute affine transformation matrix
    f1 = 1 / max(anis1,TINY);
    f2 = 1 / max(anis2,TINY);
    F  = diag([1 f1 f2]);
    % Apply rotation, and then affine transform
    RotMat(:,:,ist) = F*RotMat3*RotMat2*RotMat1;
    % Match GSLIB
    RotMat(:,:,ist) = RotMat(:,:,ist)';
 end
end               %%%%%%%%%%%%%%%%%%%%% end switch nCols