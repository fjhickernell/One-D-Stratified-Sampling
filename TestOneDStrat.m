%% Test stratified sampling convergence

gail.InitializeWorkspaceDisplay
tic
d = 5;
m=20;
nvec = 2.^(0:m)';
nlen = length(nvec);
nmax = nvec(end);
rng = 474747;

nrep = 500;
xpts(nmax,d) = 0;
cub(nmax,nrep) = 0;
h(4,1) = 0;

a = (1:d)'.^-4;
testfun = @(t) exp(t*a); trueint = exp(sum(a.^2)/2);
%testfun = @(t) t; trueint = 0;
%testfun = @(t) t.^2; trueint = 1;

vartransform = @(x) norminv(x); transpdf = @(t) prod(normpdf(t),2); label1 = 'Gaussian';
unitfunt = @(t) testfun(t) .* prod(normpdf(t),2)./transpdf(t);
unitfun = @(x) unitfunt(vartransform(x));
for r = 1:nrep
   xpts = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
   cub(:,r) = cumsum(unitfun(xpts))./(1:nmax)';
end
err = abs(trueint - cub(nvec,:));
rmse = sqrt(sum(err.*err,2));
figure
h(1) = loglog(nvec,rmse,'.','Color',MATLABBlue);
hold on
h(2) = loglog([1; nmax],rmse(nlen).*[nmax; 1],'-','Color',MATLABBlue);

df = 15; vartransform = @(x) tinv(x, df); transpdf = @(t) tpdf(t,df); label2 = ['t-' int2str(df) 'df'];
% loc = 0; sc = 1; vartransform = @(x) loc + sc*log(x./(1-x));  
%     exppart = @(t) exp(-(t-loc)/sc); transpc = @(part) part./(sc*(1 + part).^2);
%     transpdf = @(t) prod(transpc(exppart(t)),2); label2 = 'logistic';
unitfunt = @(t) testfun(t) .* prod(normpdf(t),2)./transpdf(t);
unitfun = @(x) unitfunt(vartransform(x));
for r = 1:nrep
   xpts = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
   cub(:,r) = cumsum(unitfun(xpts))./(1:nmax)';
end
err = abs(trueint - cub(nvec,:));
rmse = sqrt(sum(err.*err,2));
h(3) = loglog(nvec,rmse,'.','Color',MATLABOrange);
h(4) = loglog([1; nmax],rmse(nlen).*[nmax^(3/2); 1],'-','Color',MATLABOrange);
legend(h,{label1,'$O(n^{-1}$)',label2,'$O(n^{-3/2})$'},'box','off')
set(gca,'XTick',10.^(0:2:6),'Ytick',10.^(-4:2:4))
print('-depsc',['ConvergeRateStrat_' label2 '_d' int2str(d) '.eps'])
toc
