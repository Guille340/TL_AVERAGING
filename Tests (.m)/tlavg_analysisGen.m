%  DESCRIPTION: plots the standard error of transmission loss curves 
%  calculated with the Harrison & Harrison (1995) method against the BPO,
%  NRA and FC. To do so, it first loads a 'tlavgAnalysis_...' structure
%  generated with tlavg_analysisLog.m.
%
%  See also tlavg_ex1.m, tlavg_ex2.m, tlavg_analysisLog.m

%  REVISION 1.1
%  - Added function help
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  13 May 2020

fc = [15 50 150 500 1500];

mkdir(pwd,'Figures')
outDir = fullfile(pwd,'Figures');
nFreqs = length(fc);
bpoMin = zeros(1,nFreqs);
nralphaMax = zeros(1,nFreqs);
for m = 1:nFreqs
    fname = sprintf('tlavgAnalysis_r20kLog_f%d.mat',fc(m));
    data = load(fname,'bpo','nralpha','rmse','r');
    r0 = data.r{1}(1);
    bpo = data.bpo;
    nralpha = data.nralpha;
    rmse = data.rmse;
    [~,bpoMin(m)] = min(rmse(:,end));
    nralphaMax(m) = max(nralpha);
    colors = hsv(23);   
    
    figure(1)
    hold(gca,'on')
    plot(bpo,rmse(:,end),'Color',colors((m-1)*5+1,:),'LineWidth',1)
    title('Fitting error with BPO')
    xlabel('Bands Per Octave')
    ylabel('\sigma [dB]')
    set(gca,'XScale','log','YScale','log')
    box on
    ylim([0.1 6])
    
    figure(2)
    hold(gca,'on')
    plot(nralpha,rmse(end,:)','LineWidth',1)
    title('Fitting error with NRA')
    xlabel('No. Ranges/\alpha')
    ylabel('\sigma [dB]')
    set(gca,'XScale','log','YScale','log')
    box on
    ylim([0.1 6])
    
    figure(2*(m-1) + 3)
    inralpha = [1:6 11:5:81];
    plot(bpo,rmse(:,inralpha))
    legStr = split(sprintf('nr/\\alpha = %d,',nralpha(inralpha)),',');
    legend(legStr(1:end-1),'Location','SouthWest')
    title(sprintf('Fitting error with BPO and NRA \\rm(fc = %d Hz)',fc(m)))
    xlabel('Bands Per Octave')
    ylabel('\sigma [dB]')
    set(gca,'XScale','log','YScale','log')
    box on
    hplo = get(gca,'Children');
    for n = 1:length(hplo), hplo(n).Color = colors(n,:); end
    ylim([0.1 6])
    set(gcf,'Units','Normalized','OuterPosition',[0.2 0.2 0.6 0.6])
    set(get(gca,'Legend'),'NumColumns',2)
    figname = sprintf('Sigma vs BPO (NRA = all, fc = %d Hz)',fc(m));
    figpath = fullfile(outDir,figname);
    print(figpath,'-dpng','-r250')
    savefig(figpath)
    
    figure(2*(m-1) + 4)
    plot(nralpha,rmse')
    legStr = split(sprintf('bpo = %d,',bpo),',');
    legend(legStr(1:end-1),'Location','SouthWest')
    title(sprintf(['Fitting error with NRALPHA and BPO \\rm(fc = %d Hz, '...
        '\\sigma_{min} for bpo = %d)'],fc(m),bpoMin(m)))
    xlabel('No. Ranges/\alpha')
    ylabel('\sigma [dB]')
    set(gca,'XScale','log','YScale','log')
    box on
    hplo = get(gca,'Children');
    for n = 1:length(hplo), hplo(n).Color = colors(n,:); end
    hplo(end-bpoMin(m)+1).LineWidth = 1.5;
    ylim([0.1 6])
    set(gcf,'Units','Normalized','OuterPosition',[0.2 0.2 0.6 0.6])
    set(get(gca,'Legend'),'NumColumns',2)
    figname = sprintf('Sigma vs NRA (BPO = all, fc = %d Hz)',fc(m));
    figpath = fullfile(outDir,figname);
    print(figpath,'-dpng','-r250')
    savefig(figpath)
    
end
figure(1)
legStr1 = split(sprintf('f_c = %d Hz (NRA = %d),',...
    reshape([fc' nralphaMax']',1,10)),',');
legend(legStr1(1:end-1),'Location','SouthWest')
set(gcf,'Units','Normalized','OuterPosition',[0.2 0.2 0.6 0.6])
figname = 'Sigma vs BPO (NRA = 177, fc = all)';
figpath = fullfile(outDir,figname);
print(figpath,'-dpng','-r250')
savefig(figpath)

figure(2)
legStr2 = split(sprintf('f_c = %d Hz (BPO = %d),',...
    reshape([fc' 18*ones(1,nFreqs)']',1,10)),',');
legend(legStr2(1:end-1),'Location','SouthWest')
set(gcf,'Units','Normalized','OuterPosition',[0.2 0.2 0.6 0.6])
figname = 'Sigma vs NRA (BPO = 18, fc = all)';
figpath = fullfile(outDir,figname);
print(figpath,'-dpng','-r250')
savefig(figpath)

