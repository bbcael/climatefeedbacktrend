clear all; close all; clc;
mu = [(-.01:.0001:.03)]'; % mu relative to FINAL year -- faster computationally
load fht.mat; % time-series ensemble
r = f-h; 
fd = f;
td = t;
r = repmat(r,1,1,length(mu));
t = repmat(t,1,1,length(mu));
t = cumsum(permute(repmat(1-mu.*(49:-1:0),1,1,2237),[3 2 1]).*t,2);
t0 = cumsum(permute(repmat(1+0.*(0:49),1,1,2237),[3 2 1]).*t,2);

%%
for i = 1:size(t,1);
    i
    x = squeeze(t(i,:,:))'; y = squeeze(r(i,:,:))';
    x0 = squeeze(t0(i,:,:))';
    [~,l,b] = regression(x,y);
    [~,l0,b0] = regression(x0,y);
    res = sum((y-(l.*x+b)).^2,2);
    res0 = sum((y-(l0.*x0+b0)).^2,2);
    [R(i),indx] = min(res);
    [R0(i),indx0] = min(res0);
    dBIC(i) = (50.*log(R(i)./50)+2.*log(50))-(50.*log(R0(i)./50)+1.*log(50));
    dAIC(i) = (50.*log(R(i)./50)+2.*2)-(50.*log(R0(i)./50)+1.*2);
    L(i) = l(indx);
    M(i) = mu(indx);
    L0(i) = l(indx0);
end
clear i x y l b x0 l0 b0 res res0 indx indx0 R R0 t0 mu;
clear r t;
save results_take7_backwards.mat;

%%

% transform to find lambda 1970 and mu relative to it
L1970 = L.*(1-49.*M);
M1970 = M.*L./L1970;

