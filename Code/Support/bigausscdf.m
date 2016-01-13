function out = bigausscdf(yL,yU,rho)
% Integral over rectangle of standard bivariate Gaussian density
%
%% DESCRIPTION: bigausscdf.m
% Function to integrate STANDARD Gaussian bi-variate density function
% with correlation coefficient rho, over a rectangle specified by its
% 2 lower limits (yL) and 2 upper limits (yU).
% NOTE: Based on function mvncdf of the Statistics Toolbox
% NOTE: Uses adaptive quadrature on a transformation of the t density, 
%       based on methods developed by Drezner and Wesolowsky, and by Genz, 
%       as described in the references of function mvncdf.  
% NOTE: The default absolute error tolerance for these cases is 1e-8. 
%
%% SYNTAX:
%   out = bigausscdf(yL,yU,rho);
%
%% INPUTS:
%   yL       = (1 x 2) array [yLow1 yLow2] with lower limits of integration
%   yU       = (1 x 2) array [yHigh1 yHigh2] with upper limits for integr
%   rho      = scalar with correlation coefficient
%                              
%% OUTPUTS:
%   out      = scalar with corresponding CDF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYNTAX:
%   out = bigausscdf(yL,yU,rho);

%% CREDITS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          Phaedon Kyriakidis                             %
%                       Department of Geography                           %
%               University of California Santa Barbara                    %
%                            September 2007                               %
%                                                                         %
%        Based on function mvncdf from Matlab's Statistics Toolbox        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Fixed parameters
NARGMIN   = 3;
tol = 1e-8;

%% Error checking
if ~( nargin == NARGMIN )
    tit = ['You need ',num2str(NARGMIN),...
           ' input arguments for this function'];
    disp(' ');
    disp(tit);
    error('Check # of input arguments');
end
if numel(yL) ~= 2 || numel(yU) ~= 2
    error('You need 2 entries in yL and yU');
end
if any(any(yL > yU))
    error('Entries in yL must be <= to entries in yU');
end
S  = [1 rho; rho 1];
% Test for positive definiteness
chol(S);

%% Proceed with computation
nDim = 2; nDat = 1;
yL0 = yL; yU0 = yU;
% Compute the probability over the rectangle as sums and differences
% of integrals over semi-infinite half-rectangles.  For degenerate
% rectangles, force an exact zero by making each piece exactly zero.
equalLims = (yL==yU);
yL0(equalLims) = -Inf;
yU0(equalLims) = -Inf;
out = zeros(nDat,1);
for ii = 0:nDim
    k = nchoosek(1:nDim,ii);
    for jj = 1:size(k,1)
        y = yU0; y(:,k(jj,:)) = yL0(:,k(jj,:));
        out = out + (-1)^ii * bvncdf(y, rho, tol/4);
    end
end

end
%----------------------------------------------------------------------
function p = bvncdf(b,rho,tol)
% CDF for the bivariate normal.
%
% Implements the unnumbered equation between (3) and (4) in Section 2.2 of
% Genz (2004), integrating in terms of theta between asin(rho) and +/- pi/2,
% using adaptive quadrature.

n = size(b,1);
sqrt2 = sqrt(2);
if rho == 0
    gausscdfb = 0.5.*erfc( -b./sqrt2 );
    p = cast(prod(gausscdfb,2), superiorfloat(b,rho));
else
    if rho > 0
        tmp = min(b,[],2);
        p1 =  0.5.*erfc( -tmp./sqrt2 );
        p1(any(isnan(b),2)) = NaN;
    else
        g1 = 0.5.*erfc( -b(:,1)./sqrt2 );
        g2 = 0.5.*erfc( -b(:,2)./sqrt2 );
        p1 = g1 - g2;
        p1(p1<0) = 0; % max would drop NaNs
    end
    if abs(rho) < 1
        loLimit = asin(rho);
        hiLimit = sign(rho).*pi./2;
        p2 = zeros(size(p1),class(p1));
        for i = 1:n
            b1 = b(i,1); b2 = b(i,2);
            if isfinite(b1) && isfinite(b2)
                p2(i) = quadl(@bvnIntegrand,loLimit,hiLimit,tol);
            else
                % This piece is zero if either limit is +/- infinity.  If
                % either is NaN, p1 will already be NaN.
            end
        end
    else
        p2 = zeros(class(p1));
    end
    p = cast(p1 - p2./(2.*pi), superiorfloat(b,rho));
end

    function integrand = bvnIntegrand(theta)
        % Integrand is exp( -(b1.^2 + b2.^2 - 2*b1*b2*sin(theta))/(2*cos(theta).^2) )
        sintheta = sin(theta);
        cossqtheta = cos(theta).^2; % always positive
        integrand = exp(-((b1*sintheta - b2).^2 ./ cossqtheta + b1.^2)/2);
    end
end