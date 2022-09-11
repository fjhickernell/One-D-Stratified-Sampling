%% Test stratified sampling convergence

%% Initialize
gail.InitializeWorkspaceDisplay
tstart = tic;
tlast = tstart;
m = 20;
nvec = 2.^(0:m)';
nlen = length(nvec);
nmax = nvec(end);
rng = 474747;
nrep = 500;
cub(nmax,nrep) = 0;
cubiid(nmax,nrep) = 0;
h(6,1) = 0;


%% Make test function cases
[tf([1:3 7:9]).d] = deal(3);
[tf([4:6 10:12]).d] = deal(5);
dmax = max([tf(:).d]);
xpts(nmax,dmax) = 0;
xptsiid(nmax,dmax) = 0;
a = (1:dmax)'.^-4;
[tf([1 4 7 10]).testfun] = deal(@(t,d) exp(t*a(1:d))); 
[tf([1 4]).trueint] = deal(@(d) prod((exp(a(1:d))- 1)./a(1:d))); 
[tf([7 10]).trueint] = deal(@(d) exp(sum(a(1:d).^2)/2)); 
[tf([1 4 7 10]).testfunname] = deal('exp');
[tf([1 7]).testfuntitle] = ...
    deal(@(d)['$\exp(' coeftostr(a(1)) '\textit{t}_1 + ' ...
    coeftostr(a(2)) '\textit{t}_2 + ' ...
    coeftostr(a(3)) '\textit{t}_3)$'] );
[tf([4 10]).testfuntitle] = ...
    deal(@(d)['$\exp(' coeftostr(a(1)) '\textit{t}_1 + \cdots + ' ...
    coeftostr(a(d)) '\textit{t}_{' int2str(d) '})$'] );
[tf([2 5 8 11]).testfun] = deal(@(t,d) 1 + t*a(1:d)); 
[tf([2 5]).trueint] = deal(@(d) 1 + sum(a(1:d))/2);  
[tf([8 11]).trueint] = deal(@(d) 1);  
[tf([2 5 8 11]).testfunname] = deal('t');
[tf([2 8]).testfuntitle] = ...
    deal(@(d)['$1 + ' coeftostr(a(1)) '\textit{t}_1 + ' ...
    coeftostr(a(2)) '\textit{t}_2 + ' ...
    coeftostr(a(3)) '\textit{t}_3$']);
[tf([5 11]).testfuntitle] = ...
    deal(@(d)['$1 + ' coeftostr(a(1)) '\textit{t}_1 + \cdots + ' ...
    coeftostr(a(d)) '\textit{t}_{' int2str(d) '}$']);
[tf([3 6 9 12]).testfun] = deal(@(t,d) sum((t .* a(1:d)').^2,2)); 
[tf([3 6]).trueint] = deal(@(d)sum(a(1:d).^2)/3); 
[tf([9 12]).trueint] = deal(@(d)sum(a(1:d).^2)); 
[tf([3 6 9 12]).testfunname] = deal('tsq');
[tf([3 9]).testfuntitle] = deal(@(d)['$' coeftostr(a(1)^2) '\textit{t}_1^2 + ' ...
    coeftostr(a(2)^2) '\textit{t}_2^2 + ' ...
    coeftostr(a(3)^2) '\textit{t}_3^2$']);
[tf([6 12]).testfuntitle] = deal(@(d)['$' coeftostr(a(1)^2) '\textit{t}_1^2 + \cdots + ' ...
    coeftostr(a(d)^2) '\textit{t}_{' int2str(d) '}^2$']);
[tf(1:6).weight] = deal(@(t) 1);
[tf(1:6).weightname] = deal('uniform');
[tf(7:12).weight] = deal(@(t) prod(normpdf(t),2));
[tf(7:12).weightname] = deal('stdGauss');
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


%% Initialize
rmseSobnat(nlen,ntf) = 0;
rmseIIDnat(nlen,ntf) = 0;
rmseTrans(nlen,ntf,nvt)=0;

%whii = 1:ntf;
whii = [1 7];
%whjj = 1:nvt;
whjj = 2;
%%Perform calculations
for ii = whii
    disp([ii 0])
    tlast = tic;
    yesGauss = strcmp(tf(ii).weightname,'stdGauss');
    if yesGauss
       vartransform = @(x) norminv(x); transpdf = @(t) prod(normpdf(t),2); 
       unitfunt = @(t) tf(ii).testfun(t,tf(ii).d) .* tf(ii).weight(t)./transpdf(t);
       unitfun = @(x) unitfunt(vartransform(x));
    else %uniform
       unitfun = @(x) tf(ii).testfun(x,tf(ii).d);
    end
    %natural variable transformation with Sobol' and IID
    rmseSobnat(:,ii) = compSobolIID('Sobol',tf(ii).d,nvec,nrep,unitfun,tf(ii).trueint(tf(ii).d));   
    rmseIIDnat(:,ii) = compSobolIID('IID',tf(ii).d,nvec,nrep,unitfun,tf(ii).trueint(tf(ii).d));   
   
    if yesGauss %try alternative transforms
        for jj = whjj
           if toc(tlast) > 2
              disp([ii jj])
              tlast = tic;
           end
           unitfunt = @(t) tf(ii).testfun(t,tf(ii).d) .* tf(ii).weight(t) ./ vt(jj).transpdf(t);
           unitfun = @(x) unitfunt(vt(jj).vartransform(x));
           rmseTrans(:,ii,jj) = ...
              compSobolIID('Sobol',tf(ii).d,nvec,nrep,unitfun,tf(ii).trueint(tf(ii).d));   
        end
    end
end

datafile = ['TestRQMCImpSampResults-' datestr(now,'yy-mm-dd-HH-MM-SS')];
save(datafile,'m','ntf','nvt','nvec','nlen','nmax','tf','vt','rmseSobnat', ...
    'rmseIIDnat','rmseTrans','whii','whjj');
PlotTestRQMCImpSamp(datafile)
toc(tstart)



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


function rmse = compSobolIID(whpts,d,nvec,nrep,unitfun,trueint)
nmax = nvec(end);
nnvec = length(nvec);
cub(nnvec,nrep) = 0;
for r = 1:nrep
   if strcmp(whpts,'Sobol')
     xpts = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
   else
     xpts = rand(nmax,d);
   end
  temp = cumsum(unitfun(xpts));
  cub(:,r) = temp(nvec)./nvec;
end
err = abs((trueint - cub)./trueint);
rmse = sqrt(mean(err.*err,2));
end

