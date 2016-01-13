% only work in 2d case
function out = dist2trans(globprops,lagdist, modpars,bplot)
%% Description
%  Get transition probability given lag distance, global proportion and
%  modepars.
%% Compute rotation matrix for all structures
modpars(:,2) = modpars(:,2)./sum(modpars(:,2));
nStruct  = size(modpars,1);
nClass = numel(globprops);

%% Compute Gaussian thresholds
globprops = double(globprops(:));
globprops(globprops > (1-eps)) = 1-eps;
globprops(globprops < eps) = eps;
globprops=globprops./sum(globprops);
thresh = cumsum(globprops(1:nClass-1));
thresh = -sqrt(2).*erfcinv(2.*thresh);
thresh = [-inf; thresh; inf];

nD1= size(lagdist,1);
tpmat = cell(nD1,1);
for ii=1:nD1
    covTot = 0;
    for ist =1:nStruct
        p = modpars(ist,[1 2 4 end]);
        covTot = covTot + dist2mcov(lagdist(ii),[1 1],p,1);
    end
    for iT=1:nClass
        for iH=1:nClass
            yLowT = thresh(iT);
            yLowH = thresh(iH);
            yLow = [yLowT yLowH];
            yHighT = thresh(iT+1);
            yHighH = thresh(iH+1);
            yHigh = [yHighT yHighH];
            jointProb = bigausscdf(yLow,yHigh,covTot);
            pT = globprops(iT);
            tpmat{ii}(iT,iH)= jointProb./pT;
        end
    end
end
out=tpmat;
if bplot ==1
    tp=zeros(nD1,nClass,nClass);
    figure;
    index=0;
    for c=1:nClass
        for k=1:nClass
            index = index +1;
            subplot(nClass,nClass,index);
            for l = 1:nD1
                tp(l,c,k)=tpmat{l}(c,k);
            end
            plot(lagdist,tp(lagdist,c,k),'-');
            xlim([min(lagdist),max(lagdist)]);
            hold on;
            
        end
    end
    hold off;
    clear tp;
else
    return;
end
    



       








