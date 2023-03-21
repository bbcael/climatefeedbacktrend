% figure S3
clear all; close all; clc; load fht.mat;

for i = 1:2237; % 30-year sliding window method
    for j = 15:34;
        T = t(i,1:end-1)';
        R = diff(f(i,:)'-h(i,:)');
       [~,l(i,j),~] = regression(T((j-14):(j+15))',R((j-14):(j+15))');
    end
end
L = 1./16.0886.*prctile(l(:,15:34),[17 50 83]); L = L'; % median & 66% CI

%%

% figure S2a, S1, post-hoc justification of ansatz
% for S1, use same code to fit medians without taking residuals
% for S2a, use same code with mu = 0

clear all; close all; clc; load fht.mat;
mu = [(-.01:.001:.03)]';
r = f-h;
r = repmat(r,1,1,length(mu));
t = repmat(t,1,1,length(mu));
t = cumsum(permute(repmat(1-mu.*(49:-1:0),1,1,2237),[3 2 1]).*t,2);
for i = 1:size(t,1); % quadratic fit to residuals
    i
    x = squeeze(t(i,:,:))'; y = squeeze(r(i,:,:))';
    [~,l,b] = regression(x,y);
    res = sum((y-(l.*x+b)).^2,2);
    [R(i),indx] = min(res);
    L(i) = l(indx);
    M(i) = mu(indx);
    resid = y(indx,:)-L(i).*x(indx,:)+b(indx);
    X = x(indx,:);
    fit2 = fit(X',resid','poly2');
    ci = confint(fit2);
    U(i) = abs(((ci(2,1)-ci(1,1))./3.96));
    P(i) = fit2.p1;
end

%%

% figure 3a
% for figure 3b, do the same forward-stepping, with historical forcing, same kappa, and either lambda(1970) or lambda(t), & plot temperature difference
clear all; close all; clc; % load results from main analysis
F = 4+.3.*randn(1,2237);
K = .58+.08.*randn(1,2237);
E = 9.67+.8.*randn(1,2237);
L = L1970./16.09;

for i = 1:2237;
    T1(1,i) = 0;
    T2(1,i) = 0;
    for j = 2:140;
        T1(j,i) = T1(j-1,i) + 1./E(i).*(F(i).*j./70 - (K(i)+L(i)).*T1(j-1,i));
        T2(j,i) = T2(j-1,i) + 1./E(i).*(F(i).*j./70 - (K(i)+L(i).*(1+49.*M1970(i))).*T2(j-1,i));
    end
    [~,y_2_1(i)] = min(abs(T1(:,i)-2));
    [~,y_2_2(i)] = min(abs(T2(:,i)-2));
    [~,y_1p5_2(i)] = min(abs(T2(:,i)-1.5));
    [~,y_1p5_1(i)] = min(abs(T1(:,i)-1.5));
end
d1p5 = y_1p5_2-y_1p5_1;
d2 = y_2_2-y_2_1;

d1p5e0 = 1.5*70*(K+L)./F-1.5*70*(K+(1+49*M1970).*L)./F;
d2e0 = 2*70*(K+L)./F-2*70*(K+(1+49*M1970).*L)./F;
