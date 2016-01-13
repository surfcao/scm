% This file is to demonstrate how to use Multinomial Logistic Mixed Model
% for prediction of categorical spatial data (please see the enclosed mk.pdf for more about this method). 
% Jura dataset is used for demonstrating advantages of multinomial logistic
% mixed model over the current BME (Bayesian Maximum Entropy) method.

clear all;
load mk;

%% Jura data preprocessing
minx = min(rock(:,1));
maxx = max(rock(:,1));
miny = min(rock(:,2));
maxy = max(rock(:,2));
xpixel=(maxx-minx)/63;
ypixel=(maxy-miny)/63;
xorig = minx-0.50*xpixel;
yorig = miny-0.50*ypixel;
jura_gridspecs = [64,xorig,xpixel;64,yorig,ypixel];
jura_data = NaN(64,64);
for i = 1:size(rock,1)
    x = rock(i,1);
    y = rock(i,2);
    ix = ceil((x-xorig)/xpixel);
    iy = ceil((y-yorig)/ypixel);
    jura_data(ix,iy) = rock(i,3);
end
jura_data = reshape(jura_data, 4096,1);

%% Empirocal Transiogram 
meas3=[4 4 5;4 5 5; 4 6 5; 4 7 5; 5 4 5; 5 5 5; 5 6 5; 5 7 5; 6 4 5; 6 5 5; 6 6 5; 6 7 5;7 4 5; 7 5 5;7 6 5; 7 7 5];
ucwtp=scatter2struct2(rock,[1,2],[],[[0:0.1:3]',0.1*ones(size([0:0.1:3]'))],meas3);
tptot3=ucwtp.struct{1};
nclass=4;
fh = 100;
figure(fh);
hold on;
ip = 0;
for i = 1:nclass
    for j = 1:nclass
        ip = ip + 1;
        subplot(nclass,nclass,ip);
        hold on;
        plot(tptot3(1:30,1,(i-1)*nclass+j),tptot3(1:30,2,(i-1)*nclass+j),'g-','LineWidth',2);
        xlim([0 3]);
        txt=sprintf('Class %d-%d',i,j);
        title(txt);
    end
end
%% Apply BME on Jura dataset
kk = convhull(rock(:,1), rock(:,2)); 
jura_data_trans = jura_data;
jura_classprops = [0.2046,0.3282,0.2432,0.2239];
nClasses_jura = 4;
coords = rock(:,[1 2]);
nClasses = 4;
nCount = 5;
for i = 1:size(jura_data_trans,1)
    fprintf('computing %d of %d points... \n', i,size(jura_data_trans,1));
    centerCoords = indgeoeas2coords(i,jura_gridspecs);
    if inpolygon(centerCoords(1,1),centerCoords(1,2),rock(kk,1),rock(kk,2))==1
        dist= sum((repmat(centerCoords,size(coords,1),1)-coords).^2,2);
        [sortval, sortpos]=sort(dist,'ascend');
        Cpmf = ones(nClasses,1);
        for c = 1:nClasses
            for j = 1:nCount
                Cpmf(c) = ksr2(tptot3(1:30,1,(c-1)*nClasses+rock(sortpos(j),3)),tptot3(1:30,2,(c-1)*nClasses+rock(sortpos(j),3)),sqrt(sortval(j)),0.05)*Cpmf(c);
            end
        end
        Cpmf = Cpmf.*jura_classprops';
        Cpmf = Cpmf./sum(Cpmf);
        [junk yhat] = max(Cpmf,[],1);
        jura_data_trans(i)=yhat;
    end    
end
%% transiogram jura validation    
kk = convhull(rock_validation2(:,1), rock_validation2(:,2)); 
jura_classprops = [0.2046,0.3282,0.2432,0.2239];
nClasses_jura = 4;
srchSpecs= [0 10,10];
coords = rock(:,[1 2]);
nClasses = 4;
nCount = 10;
nEq = 0;
nTotal = 0;
for i = 1:size(rock_validation2,1)
    fprintf('computing %d of %d points... \n', i,size(rock_validation2,1));
    centerCoords = rock_validation2(i,[1,2]);
%   centerCoords = indgeoeas2coords(i,jura_gridspecs);
    if inpolygon(centerCoords(1,1),centerCoords(1,2),rock_validation2(kk,1),rock_validation2(kk,2))==1
        nTotal = nTotal + 1;
        dist= sum((repmat(centerCoords,size(coords,1),1)-coords).^2,2);
        [sortval, sortpos]=sort(dist,'ascend');
        Cpmf = ones(nClasses,1);
        for c = 1:nClasses
            for j = 1:nCount
                Cpmf(c) = ksr2(tptot3(1:30,1,(c-1)*nClasses+rock(sortpos(j),3)),tptot3(1:30,2,(c-1)*nClasses+rock(sortpos(j),3)),sqrt(sortval(j)),0.05)*Cpmf(c);
            end
        end
        Cpmf = Cpmf.*jura_classprops';
        Cpmf = Cpmf./sum(Cpmf);
        [junk yhat] = max(Cpmf,[],1);
        if rock_validation2(i,3)==yhat
            nEq = nEq + 1;
        end
    end
end

%% Kernel design for Jura dataset
meas4=[4 4 3;5 5 3; 6 6 3; 7 7 3];
ucwtp2=scatter2struct2(rock,[1,2],[],[[0:0.1:3]',0.1*ones(size([0:0.1:3]'))],meas4);
tptot4=ucwtp2.struct{1};
nclass=4;
fh = 100;
figure(fh);
ip = 0;
for i = 1:nclass
        ip = ip + 1;
        subplot(nclass,1,ip);
        hold on;
        plot(tptot4(1:30,1,i),tptot4(1:30,2,i),'g-','LineWidth',2);
        xlim([0 3]);
        txt=sprintf('Class %d-%d',i,i);
        title(txt);
end
beta=[0,0.25,1];
covObj=@covmodgauss;
betahat=nlinfit(tptot4(1:30,1,1),tptot4(1:30,2,1),covObj,beta);
modparsarray{1}=[3,betahat(2),0,betahat(3),betahat(3),0];
betahat=nlinfit(tptot4(1:30,1,2),tptot4(1:30,2,2),covObj,beta);
modparsarray{2}=[3,betahat(2),0,betahat(3),betahat(3),0];
betahat=nlinfit(tptot4(1:30,1,3),tptot4(1:30,2,3),covObj,beta);
modparsarray{3}=[3,betahat(2),0,betahat(3),betahat(3),0];
betahat=nlinfit(tptot4(1:30,1,4),tptot4(1:30,2,4),covObj,beta);
modparsarray{4}=[3,betahat(2),0,betahat(3),betahat(3),0];
%% Mixed Model
kk = convhull(rock(:,1), rock(:,2)); 
jura_data_MM= jura_data;
jura_classprops = [0.2046,0.3282,0.2432,0.2239];
options.Display = 'full';
options.MaxIter=5000;
options.MaxFunEvals=5000;
nClasses = 4;
rocksize=size(rock,1);
neighCoords = rock(:,[1,2]);
kernel2 = zeros(rocksize,rocksize);
knownSize=size(neighCoords,1);
for i = 1:knownSize
    fprintf('computing %d of %d points... \n', i,knownSize);
    centerCoords = neighCoords(i,:);
    for k =1:nClasses
        kernel2(i,:)= kernel2(i,:)+ jura_classprops(k)*coords2mcov2(centerCoords,modparsarray{k},neighCoords);
    end
end
lambda = 0.01;
funObj = @(w)MLogisticLoss2(w,kernel2,rock(:,3),nClasses);
u= minFunc(@PenalizedKernelMatrix,randn((rocksize)*(nClasses-1),1),options,kernel2, nClasses-1,funObj,lambda);
% u= fminunc(@PenalizedKernelMatrix,randn((rocksize)*(nClasses-1),1),options,kernel2, nClasses-1,funObj,lambda);
u = reshape(u, rocksize , nClasses-1);
u = [u zeros(rocksize,1)];

for i = 1:size(jura_data_MM,1)
    fprintf('computing %d of %d points... \n', i,size(jura_data_MM,1));
    centerCoords = indgeoeas2coords(i,jura_gridspecs);
    ukernel = zeros(1,rocksize);
    if inpolygon(centerCoords(1,1),centerCoords(1,2),rock(kk,1),rock(kk,2))==1
        for k =1:nClasses
             ukernel(1,1:rocksize)= ukernel(1,1:rocksize)+ jura_classprops(k)*coords2mcov2(centerCoords,modparsarray{k},neighCoords);
         end
        [junk yhat] = max(ukernel*u,[],2);
        jura_data_MM(i)=yhat;
    end 
end
%% Jura Mixed Model validation
kk = convhull(rock_validation2(:,1), rock_validation2(:,2)); 
jura_classprops = [0.2046,0.3282,0.2432,0.2239];
% jura_classprops = jura_classprops.^2;
nVars = nClasses;
options.Display = 'full';
coords = rock(:,[1 2]);
nClasses = 4;
nCount = 20;
rocksize=size(rock,1);
neighCoords = rock(:,[1,2]);
kernel3 = zeros(rocksize,rocksize);
for i = 1:size(neighCoords,1)
    fprintf('computing %d of %d points... \n', i,knownSize);
    centerCoords = neighCoords(i,:);
    for k =1:nClasses
        kernel3(i,:)= kernel3(i,:)+ jura_classprops(k)*coords2mcov2(centerCoords,modparsarray{k},neighCoords);
    end
end
lambda = 0.01;
funObj = @(u)MLogisticLoss2(u,kernel3,rock(:,3),nClasses);
u= minFunc(@PenalizedKernelMatrix,randn((rocksize)*(nClasses-1),1),options,kernel3, nClasses-1,funObj,lambda);
%u= fminunc(@PenalizedKernelMatrix,randn((rocksize)*(nClasses-1),1),options,kernel3, nClasses-1,funObj,lambda);
u = reshape(u, rocksize, nClasses-1);
u = [u zeros(rocksize,1)];
nEq = 0;
nTotal = 0;
for i = 1:size(rock_validation2,1)
    fprintf('computing %d of %d points... \n', i,size(rock_validation2,1));
    ukernel = zeros(1,rocksize);
    centerCoords =  rock_validation2(i,[1,2]);
    if inpolygon(centerCoords(1,1),centerCoords(1,2),rock_validation2(kk,1),rock_validation2(kk,2))==1
        nTotal = nTotal + 1;
        for k =1:nClasses
             ukernel(1,1:rocksize)= ukernel(1,1:rocksize)+ jura_classprops(k)*coords2mcov2(centerCoords,modparsarray{k},neighCoords);
        end
        [junk yhat] = max(ukernel*u,[],2);
        if rock_validation2(i,3)==yhat
            nEq = nEq + 1;
        end
    end    
end

