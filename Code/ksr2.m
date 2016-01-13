function r=ksr2(x,y,x_test,h)
% function to fit empirical transiogram values.
% Guofeng Cao (cao@geog.ucsb.edu)
% Based on Yi Cao's (Cranfield University) ksr method

% x: input lag distances
% y: input empirical transiogram values, same size as x
% x_test: lag distance to be estimated
% h: band width

% Check input and output
error(nargchk(2,6,nargin));
error(nargoutchk(0,1,nargout));
if numel(x)~=numel(y)
    error('x and y are in different sizes.');
end

x=x(:);
y=y(:);

r.n=length(x);
if nargin<4
    % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
    hx=median(abs(x-median(x)))/0.6745*(4/3/r.n)^0.2;
    hy=median(abs(y-median(y)))/0.6745*(4/3/r.n)^0.2;
    h=sqrt(hy*hx);
    N = size(x,1)
    if h<sqrt(eps)*N
         error('There is no enough variation in the data. Regression is meaningless.')
     end
elseif ~isscalar(h)
    error('h must be a scalar.')
end
% Gaussian kernel function
kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);


z=kerf((x_test-x)/h);
sumz=sum(z);
if sumz==0
    r =0;
else
    r=sum(z.*y)/sumz;
end
return;
