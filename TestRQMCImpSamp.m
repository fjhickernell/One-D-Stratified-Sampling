%% Test stratified sampling convergence

gail.InitializeWorkspaceDisplay
tic
d = 3;
m=8;
nvec = 2.^(0:m)';
nlen = length(nvec);
nmax = nvec(end);
rng = 474747;
nrep = 500;
cub(nmax,nrep) = 0;
cubiid(nmax,nrep) = 0;
h(6,1) = 0;


%% Make test function cases
[tf(1:3).d] = deal(3);
[tf(4:6).d] = deal(5);
dmax = max([tf(:).d]);
xpts(nmax,dmax) = 0;
xptsiid(nmax,dmax) = 0;
a = (1:dmax)'.^-4;
[tf([1 4]).testfun] = deal(@(t,d) exp(t*a(1:d))); 
[tf([1 4]).trueint] = deal(@(d) exp(sum(a(1:d).^2)/2)); 
[tf([1 4]).testfunname] = deal('exp');
[tf([1 4]).testfuntitle] = ...
    deal(@(d)['$\exp(' coeftostr(a(1)) '\textit{t}_1 + \cdots + ' ...
    coeftostr(a(d)) '\textit{t}_{' int2str(d) '})$'] );
[tf([2 5]).testfun] = deal(@(t,d) 1 + t*a(1:d)); 
[tf([2 5]).trueint] = deal(@(d) 1);  
[tf([2 5]).testfunname] = deal('t');
[tf([2 5]).testfuntitle] = ...
    deal(@(d)['$1 + ' coeftostr(a(1)) '\textit{t}_1 + \cdots + ' ...
    coeftostr(a(d)) '\textit{t}_{' int2str(d) '}$']);
[tf([3 6]).testfun] = deal(@(t,d) sum((t .* a(1:d)').^2,2)); 
[tf([3 6]).trueint] = deal(@(d)sum(a(1:d).^2)); 
[tf([3 6]).testfunname] = deal('tsq');
[tf([3 6]).testfuntitle] = deal(@(d)['$' coeftostr(a(1)^2) '\textit{t}_1^2 + \cdots + ' ...
    coeftostr(a(d)^2) '\textit{t}_{' int2str(d) '}^2$']);
ntf = length(tf);


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
nvt = length(vt);
rmse(m+1,nvt)=0;



for ii = 1:ntf
    vartransform = @(x) norminv(x); transpdf = @(t) prod(normpdf(t),2); label1 = 'Gaussian';
    unitfunt = @(t) tf(ii).testfun(t,tf(ii).d) .* prod(normpdf(t),2)./transpdf(t);
    unitfun = @(x) unitfunt(vartransform(x));
    for r = 1:nrep
        xpts(:,1:tf(ii).d) = net(scramble(sobolset(tf(ii).d),'MatousekAffineOwen'),nmax);
        cub(:,r) = cumsum(unitfun(xpts(:,1:tf(ii).d)))./(1:nmax)';
        xpts(:,1:tf(ii).d) = rand(nmax,tf(ii).d);
        cubiid(:,r) = cumsum(unitfun(xpts(:,1:tf(ii).d)))./(1:nmax)';
    end
    err = abs(tf(ii).trueint(tf(ii).d) - cub(nvec,:))/tf(ii).trueint(tf(ii).d);
    rmsegauss = sqrt(mean(err.*err,2));
    err = abs(tf(ii).trueint(tf(ii).d) - cubiid(nvec,:))/tf(ii).trueint(tf(ii).d);
    rmseiid = sqrt(mean(err.*err,2));

    for jj = 1:nvt
        unitfunt = @(t) tf(ii).testfun(t,tf(ii).d) .* prod(normpdf(t),2)./vt(jj).transpdf(t);
        unitfun = @(x) unitfunt(vt(jj).vartransform(x));
        for r = 1:nrep
            xpts(:,1:tf(ii).d) = net(scramble(sobolset(tf(ii).d),'MatousekAffineOwen'),nmax);
            cub(:,r) = cumsum(unitfun(xpts(:,1:tf(ii).d)))./(1:nmax)';
        end
        err = abs(tf(ii).trueint(tf(ii).d) - cub(nvec,:))/tf(ii).trueint(tf(ii).d);
        rmse(:,jj) = sqrt(mean(err.*err,2));
        figii = (ii-1)*nvt+jj;
        figure(figii)
        h(1) = loglog(nvec,rmsegauss,'.','Color',MATLABBlue);
        hold on
        h(2) = loglog([1; nmax],rmsegauss(nlen).*[nmax; 1],'-','Color',MATLABBlue);
        h(3) = loglog(nvec,rmseiid,'.','Color',MATLABGreen);
        h(4) = loglog([1; nmax],rmseiid(nlen).*[sqrt(nmax); 1],'-','Color',MATLABGreen);
        h(5) = loglog(nvec,rmse(:,jj),'.','Color',MATLABOrange);
        h(6) = loglog([1; nmax],rmse(nlen,jj).*[nmax^(3/2); 1],'-','Color',MATLABOrange);
        legend(h(1:6),{label1,'$O(n^{-1}$)',['IID ' label1],'$O(n^{-1/2}$)', ...
            vt(jj).label2,'$O(n^{-3/2})$'},'box','off')
        axis([1 3e6 1e-6 1e3])
        set(gca,'XTick',10.^(0:2:6),'Ytick',10.^(-6:2:2))
        title(tf(ii).testfuntitle(tf(ii).d))
        ylabel('Relative Error')
        print('-depsc',['ConvergeRateStrat_' tf(ii).testfunname '_' ...
            vt(jj).label2 '_d' int2str(tf(ii).d) '_m' int2str(m) '.eps'])
        %if figii == 3, keyboard, end
    end
end
toc



function out = coeftostr(in,dig)
if nargin < 2
    dig = 3;
end
if in == 1
    out = '';
else
    out = num2str(in,dig);
end

end
