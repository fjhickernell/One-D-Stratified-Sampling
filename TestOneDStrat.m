%% Test stratified sampling convergence


gail.InitializeWorkspaceDisplay
d = 1;
a = ones(d,1);
%testfun = @(x) exp(x*a); trueint = exp(a - a^2/2);
testfun = @(x) x; trueint = 0;
unitfun = @(x) testfun(norminv(x));
nvec = 2.^(0:12)';
nlen = length(nvec);
nmax = nvec(end);
rng = 474747;

nrep = 50;
xpts(nmax,d) = 0;
cub(nmax,nrep) = 0;
for r = 1:nrep
   xpts = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
   cub(:,r) = cumsum(unitfun(xpts))./(1:nmax)';
end
err = abs(trueint - cub(nvec,:));
rmse = sqrt(sum(err.*err,2));
figure
loglog(nvec,rmse,'.')
hold on
plot([1; nmax],rmse(nlen).*[nmax^(3/2) nmax; 1 1],'-')

