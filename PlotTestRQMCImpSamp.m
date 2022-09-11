function PlotTestRQMCImpSamp(datafile)
close all
gail.InitializeDisplay
load(datafile)
h(6,1) = 0;
label0 = 'Sobol''';
label1 = 'natural';

for ii = whii
   yesGauss = strcmp(tf(ii).weightname,'stdGauss');
   for jj = whjj
      %if (jj == 1) || yesGauss %try alternative transforms
      figii = (ii-1)*nvt+jj;
      figure(figii)
      h(1) = loglog(nvec,rmseSobnat(:,ii),'.','Color',MATLABBlue);
      hold on
      h(2) = loglog([1; nmax],rmseSobnat(nlen,ii).*[nmax; 1],'-','Color',MATLABBlue);
      h(3) = loglog(nvec,rmseIIDnat(:,ii),'.','Color',MATLABGreen);
      h(4) = loglog([1; nmax],rmseIIDnat(nlen,ii).*[sqrt(nmax); 1],'-','Color',MATLABGreen);
      if yesGauss
         h(5) = loglog(nvec,rmseTrans(:,ii,jj),'.','Color',MATLABOrange);
         h(6) = loglog([1; nmax],rmseTrans(nlen,ii,jj).*[nmax^(3/2); 1],'-','Color',MATLABOrange);
         legend(h(1:6),{[label0 ' ' label1],'$O(n^{-1}$)',['IID ' label1], ...
             '$O(n^{-1/2}$)', [label0 ' ' vt(jj).label2],'$O(n^{-3/2})$'}, ...
            'box','off','location','southwest')
      else
         h(5) = loglog([1; nmax],rmseSobnat(nlen,ii).*[nmax^(3/2); 1],'-','Color',MATLABOrange);
         legend(h(1:5),{[label0 ' ' label1],'$O(n^{-1})$',['IID ' label1],'$O(n^{-1/2})$', ...
            '$O(n^{-3/2})$'},'box','off','location','southwest')
      end
      axis([1 3e6 1e-16 1e3])
      set(gca,'XTick',10.^(0:3:6),'Ytick',10.^(-16:3:2))
        %title(tf(ii).testfuntitle(tf(ii).d))
        ylabel('Relative RMS Error')
        if yesGauss
           print('-depsc',['ConvergeRateStrat_' tf(ii).testfunname '_' ...
               tf(ii).weightname '_' ... 
               vt(jj).label2 '_d' int2str(tf(ii).d) '_m' int2str(m) '.eps'])
        else
           print('-depsc',['ConvergeRateStrat_' tf(ii).testfunname '_' ...
               tf(ii).weightname '_natural_d' int2str(tf(ii).d) ...
               '_m' int2str(m) '.eps'])
        end
        %if figii == 3, keyboard, end
      %end
   end
end
