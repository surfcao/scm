function [nll,g,H] = PenalizedKernelMatrix(w,K,nCols,gradFunc,lambda,varargin)

if nargout <= 1
    [nll] = gradFunc(w,varargin{:});
elseif nargout == 2
    [nll,g] = gradFunc(w,varargin{:});
else
    [nll,g,H] = gradFunc(w,varargin{:});
end

nInstances = size(K,1);
% temporily
w = reshape(w,[nInstances nCols]);

% for non square matrix
%[r,c]=size(K);
%s = max(r,c);
%KK = zeros(s);
%KK(1:r,1:c)=K;
%K=KK;


for i = 1:nCols
    nll = nll+lambda*sum(w(:,i)'*K*w(:,i));
end

if nargout > 1
    g = reshape(g,[nInstances nCols]);
    for i = 1:nCols
    g(:,i) = g(:,i) + 2*lambda*K*w(:,i);
%       g(:,i) = g(:,i) + 0.5*lambda*K*w(:,i);
    end
    g = g(:);
end
return;