%% Test stratified sampling convergence

gail.InitializeWorkspaceDisplay
tic
d = 3;
m=10;
nvec = 2.^(0:m)';
nlen = length(nvec);
nmax = nvec(end);
rng = 474747;
nrep = 500;
cub(nmax,nrep) = 0;
h(4,1) = 0;
a = (1:d)'.^-4;
dmax = 5;
xpts(nmax,dmax) = 0;

%% Make test function cases
[tf(1:3).d] = deal(3);
[tf(4:6).d] = deal(5);
[tf([1 4]).testfun] = deal(@(t) exp(t*a)); 
[tf([1 4]).trueint] = deal(exp(sum(a.^2)/2)); 
[tf([1 4]).testfunname] = deal('exp');
[tf([2 5]).testfun] = deal(@(t) t*a); 
[tf([2 5]).trueint] = deal(0);  
[tf([2 5]).testfunname] = deal('t');
[tf([3 6]).testfun] = deal(@(t) sum((t .* a').^2,2)); 
[tf([3 6]).trueint] = deal(sum(a.^2)); 
[tf([3 6]).testfunname] = deal('tsq');
ntf = 6;
return

%% Make variable transformations
df = 10; 
    vt(1).vartransform = @(x) tinv(x, df);
    vt(1).transpdf = @(t) prod(tpdf(t,df),2); 
    vt(1).label2 = ['t ' int2str(df) ' df'];
loc = 0; sc = 1; 
    vt(2).vartransform = @(x) loc + sc*log(x./(1-x));
    exppart = @(t) exp(-(t-loc)/sc); transpc = @(part) part./(sc*(1 + part).^2);
    vt(2).transpdf = @(t) prod(transpc(exppart(t)),2); 
    vt(2).label2 = 'logistic';
ntv = 2;

for ii = 1:ntf
    vartransform = @(x) norminv(x); transpdf = @(t) prod(normpdf(t),2); label1 = 'Gaussian';
    unitfunt = @(t) tf(ii).testfun(t) .* prod(normpdf(t),2)./transpdf(t);
    unitfun = @(x) unitfunt(vartransform(x));
    for r = 1:nrep
       xpts(:,tf(ii).d) = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
       cub(:,r) = cumsum(unitfun(xpts(:,tf(ii).d)))./(1:nmax)';
    end
    err = abs(trueint - cub(nvec,:));
    rmsegauss = sqrt(sum(err.*err,2));
    figure(ii)
    h(1) = loglog(nvec,rmse,'.','Color',MATLABBlue);
    hold on
    h(2) = loglog([1; nmax],rmse(nlen).*[nmax; 1],'-','Color',MATLABBlue);

    for jj = 1:ntv
        unitfunt = @(t) tf(ii).testfun(t) .* prod(normpdf(t),2)./vt(jj).transpdf(t);
        unitfun = @(x) unitfunt(vt(jj).vartransform(x));
        for r = 1:nrep
           xpts(:,tf(ii).d) = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
           cub(:,r) = cumsum(unitfun(xpts(:,tf(ii).d)))./(1:nmax)';
        end
        err = abs(trueint - cub(nvec,:));
        rmse = sqrt(sum(err.*err,2));
        h(3) = loglog(nvec,rmse,'.','Color',MATLABOrange);
        h(4) = loglog([1; nmax],rmse(nlen).*[nmax^(3/2); 1],'-','Color',MATLABOrange);
        legend(h,{label1,'$O(n^{-1}$)',label2,'$O(n^{-3/2})$'},'box','off')
        set(gca,'XTick',10.^(0:2:6),'Ytick',10.^(-4:2:4))
        print('-depsc',['ConvergeRateStrat_' tf(ii).testfunname '_' vt(jj).label2 '_d' int2str(tf(ii).d) '.eps'])
    end
end
toc

