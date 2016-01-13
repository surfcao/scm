% This file is to demonstrate how to use Multinomial Logistic Mixed Model
% for prediction of categorical spatial data (please see the enclosed mk.pdf for more about this method). 
% Syntheic dataset is used for demonstrating advantages of multinomial logistic
% mixed model over the current BME (Bayesian Maximum Entropy) method.

clear all;
load mk;

% %%synthetic data generation.
% nclass =3;
% classprops=[0.3,0.45,0.25];
% classlabels = [1,2,3];
% nsize = 64;
% gridspecs = [nsize,0.5, 1;nsize,0.5, 1];
% modpars = [0,0.001,0,10,10,0;2,0.999,0,10,10,0];
% nsim = 10;
% 
% % sim = cholsim(1,gridspecs,modpars, 0, [nsim,234,1]);
% sim = mafftsim(gridspecs,[],modpars,0,[nsim,234,1]);
% % figure;
% % rastermap(gridspecs, sim2,2,'cont');
% 
% %% converting simulated continuous data to class labels
% simcat = simgauss2class(sim,classlabels,classprops);
% %% Sample the reference data
% r=fix(rand(8000,1)*4095) + 1;
% X = simcat(:,1);
% unknownIndex = unique(r);
% knownIndex = setdiff(ind,unknownIndex);
% X(unknownIndex)=NaN;

%% empirical transiogram 
meas2=[1 1 5;1 2 5; 1 3 5; 2 1 5; 2 2 5; 2 3 5; 3 1 5; 3 2 5; 3 3 5];
nsize = 64;
nclass = 3;
tptot=zeros(nsize+1,10,9);
indic=class2indic(X,1);
ucwtp=raster2structdir(indic,gridspecs,[1 0;0 1],[nsize,nsize],meas2,1);
tptot=ucwtp.struct{1};
fh = 100;
figure(fh);
ip = 0;
for i = 1:nclass
    for j = 1:nclass
        ip = ip + 1;
        subplot(nclass,nclass,ip);
        hold on;
        plot(tptot(1:50,1,(i-1)*nclass+j),tptot(1:50,2,(i-1)*nclass+j),'g-','LineWidth',1);
        xlim([0 50]);
        txt=sprintf('Class %d-%d',i,j);
        title(txt);
    end
end
%% transiogram prediction
X_trans=X;
srchSpecs= [0 64 64];
ndMax = 20;
nClasses = 3;

props = zeros(size(X,1),3);    
for i = 1:size(X,1)
    fprintf('computing %d of %d points... \n', i,size(X,1)); 
    IndexClose = searchnodes2(X,gridspecs,i,srchSpecs,ndMax);
    ValClose = X(IndexClose,:);    
    neighCoords = indgeoeas2coords(IndexClose,gridspecs);
    centerCoords = indgeoeas2coords(i,gridspecs);
    nCount = size(IndexClose,1);
    dist= coords2sqdist(centerCoords,neighCoords);
    Cpmf = ones(nClasses,1);
    for c = 1:nClasses
        for j = 1:nCount
            Cpmf(c) = ksr2(tptot(1:63,1,(c-1)*nClasses+ValClose(j)),tptot(1:63,2,(c-1)*nClasses+ValClose(j)),sqrt(dist(j)),0.1)*Cpmf(c);
        end
    end
    Cpmf = Cpmf.*classprops';
    Cpmf = Cpmf./sum(Cpmf);
    [junk yhat] = max(Cpmf,[],1);
    props(i,:)=Cpmf(:);
    X_trans(i)=yhat;
end
correct_rate_bme = size(find(simcat(:,1)-X_trans==0),1)/4096;

%% Multinomial Logistic Mixed Model: Covariance Kernel Desgin
%% 
meas2=[1 1 3;1 2 3; 1 3 3; 2 1 3; 2 2 3; 2 3 3; 3 1 3; 3 2 3; 3 3 3];
nsize = 64;
nclass = 3;
tptot2=zeros(nsize+1,10,9);
indic=class2indic(X,1);
ucwtp=raster2structdir(indic,gridspecs,[1 0;0 1],[nsize,nsize],meas2,1);
tptot2=ucwtp.struct{3};
fh = 100;
figure(fh);
ip = 0;
for i = 1:nclass
    for j = 1:nclass
        ip = ip + 1;
        subplot(nclass,nclass,ip);
        hold on;
        plot(tptot2(1:50,1,(i-1)*nclass+j),tptot2(1:50,2,(i-1)*nclass+j),'g-','LineWidth',1);
        xlim([0 50]);
        txt=sprintf('Class %d-%d',i,j);
        title(txt);
    end
end
beta=[0,0.25,10];
covObj=@covmodgauss;
betahat=nlinfit(tptot2(1:63,1,1),tptot2(1:63,2,1),covObj,beta);
modparsarray{1}=[3,betahat(2),0,betahat(3),betahat(3),0];
betahat=nlinfit(tptot2(1:63,1,5),tptot2(1:63,2,5),covObj,beta);
modparsarray{2}=[3,betahat(2),0,betahat(3),betahat(3),0];
betahat=nlinfit(tptot2(1:63,1,9),tptot2(1:63,2,9),covObj,beta);
modparsarray{3}=[3,betahat(2),0,betahat(3),betahat(3),0];
% betahat=nlinfit(tptot2(1:63,1,4),tptot2(1:63,2,4),covObj,beta);
% modparsarray{4}=[3,betahat(2),0,betahat(3),betahat(3),0];
%%
X_MM=X;
labels =X(knownIndex);
nClasses = 3;
srchSpecs= [0 64 64];
knownSize = size(knownIndex,2);
options.Display = 'full';
props_kernel=zeros(size(X,1),3);
kernel=zeros(knownSize,knownSize);
neighCoords = indgeoeas2coords(knownIndex', gridspecs);
%modparsarray={[2,0.1574,0,5.4397,5.4397,0];[2,0.1368,0,4.3119,4.3119,0];[2,0.1138,0,10.6141,10.6141,0]}
for i = 1:knownSize
    fprintf('computing %d of %d points... \n', i,knownSize);
    centerCoords = indgeoeas2coords(knownIndex(i),gridspecs);
    for k =1:nClasses
        kernel(i,:)= kernel(i,:)+ classprops(k)*coords2mcov2(centerCoords,modparsarray{k},neighCoords);
    end
end

lambda = 0.01;
funObj = @(w)MLogisticLoss2(w,kernel,labels,nClasses);
u= minFunc(@PenalizedKernelMatrix,randn((knownSize)*(nClasses-1),1),options,kernel, nClasses-1,funObj,lambda);
% u= fminunc(@PenalizedKernelMatrix,randn((knownSize)*(nClasses-1),1),options,kernel, nClasses-1,funObj,lambda);
u= reshape(u,[knownSize nClasses-1]);
u = [u zeros(knownSize,1)];

for i = 1:size(X)
    uxx=zeros(1,knownSize);
    fprintf('computing %d of %d points... \n', i,size(simcat,1));
    centerCoords = indgeoeas2coords(i,gridspecs);
    for k =1:nClasses
        uxx(1,1:knownSize)= uxx(1,1:knownSize)+ classprops(k)*coords2mcov2(centerCoords,modparsarray{k},neighCoords);
    end
    expon = uxx*u;
    props_kernel(i,:)=exp(expon)/sum(exp(expon),2);
    [junk yhat] = max(expon,[],2);
    X_MM(i)=yhat;
end
correct_rate_mixed= size(find(simcat(:,1)-X_MM==0),1)/4096;

%% Reproduced Transiogram Check.
meas2=[1 1 5;1 2 5; 1 3 5; 2 1 5; 2 2 5; 2 3 5; 3 1 5; 3 2 5; 3 3 5];
nsize = 64;
nclass = 3;
repro_tptot_ref=zeros(nsize+1,10,9);
repro_tptot_trans=zeros(nsize+1,10,9);
repro_tptot_kernel=zeros(nsize+1,10,9);

indic=class2indic(simcat,1);
ucwtp=raster2structdir(indic,gridspecs,[1 0;0 1],[nsize,nsize],meas2,1);
repro_tptot_ref=ucwtp.struct{3};
indic=class2indic(X_trans,1);
ucwtp=raster2structdir(indic,gridspecs,[1 0;0 1],[nsize,nsize],meas2,1);
repro_tptot_trans=ucwtp.struct{3};
indic=class2indic(X_MM,1);
ucwtp=raster2structdir(indic,gridspecs,[1 0;0 1],[nsize,nsize],meas2,1);
repro_tptot_kernel=ucwtp.struct{3};

fh = 100;
figure(fh);
ip = 0;
for i = 1:nclass
    for j = 1:nclass
        ip = ip + 1;
        subplot(nclass,nclass,ip);
        hold on;
        plot(repro_tptot_ref(1:50,1,(i-1)*nclass+j),repro_tptot_ref(1:50,2,(i-1)*nclass+j),'r-','LineWidth',1);
        plot(repro_tptot_trans(1:50,1,(i-1)*nclass+j),repro_tptot_trans(1:50,2,(i-1)*nclass+j),'b-','LineWidth',1);
        plot(repro_tptot_kernel(1:50,1,(i-1)*nclass+j),repro_tptot_kernel(1:50,2,(i-1)*nclass+j),'g-','LineWidth',1);
        hold off;
        xlim([0 50]);
        txt=sprintf('Class %d-%d',i,j);
        title(txt);
    end
end
return;