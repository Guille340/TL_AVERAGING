%  DESCRIPTION: processes the band-averaged transmission loss curves using
%  the multi-frequency integration and the Harrison & Harrison (1995)
%  methods for multiple bandwidths (i.e. BPO values) and multiple range
%  points per averaging window. It loads single-frequency transmission loss 
%  curves withing one octave previously processed with function tlavg_ex1.m
%  or tlavg_ex2.m scripts. Note that tlavg_analysisLog only accepts results 
%  that were produced with tlavg_ex1.m as exponentially spaced range points 
%  (i.e.from a file of the form tlavgData_r20kLog_f<FREQ>.mat, where <FREQ> 
%  is the central frequency of the band).
%
%  See also tlavg_ex1.m, tlavg_ex2.m, tlavg_analysisLin.m

%  REVISION 1.1
%  - Added function help
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  10 May 2020


% Input Parameters
bpo = 1:18;
numRangesPerAlphaRes = 1:2:177;
outDir = pwd;
fname = 'tlavgData_r20kLog_f15.mat';

% Load Data
tic
tlavgData = load(fname);
fc = tlavgData.config.fc; % central frequency [Hz]
nFreqs = tlavgData.config.nFreqsPerOct; % number of frequencies in 3 octaves
nRanges = tlavgData.config.nRanges; % number of ranges
zs = tlavgData.config.zs; % source depth [m]
zr = tlavgData.config.zr; % receiver depth [m]
hw = tlavgData.config.hw; % water depth [m]
Rb = tlavgData.config.Rb; % bottom reflection factor (Rb = 0, Lloyd Mirror)
Rs = tlavgData.config.Rs; % sea surface reflection factor
r = tlavgData.r';
tl = abs(tlavgData.tl);
tlf = abs(tlavgData.tlf);

% Frequency Vector
f1 = fc*2^(-3/2);
f2 = fc*2^(3/2);
f = f1:(f2-f1)/(nFreqs-1):f2;

% Calculate Maximum Common Number of Frequencies for All Bands
bpoMax = max(bpo);
alpha = (2^(1/(2*bpoMax)) - 2^(-1/(2*bpoMax)))/1.665;
numFreqsPerAlphaRes = length(find(f >= fc*(1-alpha/2) & f <= fc*(1+alpha/2)));

% Create Folder
fol1 = 'Figures';
mkdir(outDir,fol1)

% Initialise Variables
nBands = length(bpo);
nRangesPerAlphaRes = length(numRangesPerAlphaRes); 
rd = cell(nBands,nRangesPerAlphaRes);
fd = cell(nBands,1);
tlfd = cell(nBands,1);
tld_mfavg = cell(nBands,nRangesPerAlphaRes);
tld_hhavg = cell(nBands,nRangesPerAlphaRes);
rmse = zeros(nBands,nRangesPerAlphaRes);
for m = 1:nBands
    % Create Folder
    fol2 = sprintf('bpo = %d',bpo(m));
    mkdir(fullfile(outDir,fol1),fol2)
    
    % Bandwidth
    alpha = (2^(1/(2*bpo(m))) - 2^(-1/(2*bpo(m))))/1.665;
    bw = 1.665*alpha*fc;
    
    % Resample Vector of Frequencies
    numFreqsPerAlpha = length(find(f >= fc*(1-alpha/2) & f <= fc*(1+alpha/2)));
    nFreqsRes = round(nFreqs * numFreqsPerAlphaRes/numFreqsPerAlpha);
    fStepRes = (f(end) - f(1))/(nFreqsRes-1);
    fd0 = f(1):fStepRes:f(end);
    ifd = interp1(f,1:length(f),fd0,'nearest','extrap');
    fd{m} = f(ifd);
    
    for n = 1:nRangesPerAlphaRes
        % Display Progress
        fprintf('Band %d/%d, Range %d/%d (%s)\n',m,nBands,n,...
            nRangesPerAlphaRes,datestr(toc/86400,'HH:MM:SS'))
        
        % Resample Vector of Ranges
        numRangesPerAlpha = length(find(r < r(1) + alpha*r(1)));
        nRangesRes = round(nRanges * numRangesPerAlphaRes(n)/numRangesPerAlpha);
        rStepRes = (r(end)/r(1))^(1/(nRangesRes-1));
        rd0 = [r(1) r(1)*rStepRes.^(1:nRangesRes-1)]';
        ird = interp1(r,1:length(r),rd0,'nearest','extrap');
        rd{m,n} = r(ird);
        
        % Resampled TL
        tlfd{m} = abs(tlf(ird)); % TL (single frequency)
        tld = tl(ird,ifd); % TL (multiple frequencies)
            
        % Multi-Frequency TL Average
        J1 = abs(1./tld).^2;
        gaussKernel = exp(-((fd{m}-fc)./(alpha*fc)).^2);
        Jb1 = sum(J1 .* gaussKernel,2)/sum(gaussKernel);
        tld_mfavg{m,n} = sqrt(1./Jb1);

        % Harrison & Harrison (1995) TL Average
        tld_hhavg{m,n} = tlavg(rd{m,n},tlfd{m},fc,bw);

        % Error
        rmse(m,n) = sqrt(sum(((20*log10(tld_mfavg{m,n} ./ ...
            tld_hhavg{m,n})).^2))/length(rd{m,n}));

%         fprintf('%d/%d (%d)\n',numRangesPerAlpha,numRangesPerAlphaRes(n),numRangesPerAlphaRes(n)<=numRangesPerAlpha)
%         fprintf('%d/%d (%d)\n',nRanges,nRangesRes,nRangesRes<=nRanges)
%         numRangesPer3AlphaRes_temp = length(find(rRes < rRes(1) + 3*alpha*rRes(1)));
%         fprintf('m = %d, n = %d (%d)\n',m,n,isequal(numRangesPer3AlphaRes(n), numRangesPer3AlphaRes_temp))
%         
        % Plot
        plotFig = any(numRangesPerAlphaRes(n) == [1:2:11 21:10:171]);
        if plotFig
            hfig = figure;
            hold on
            plot(rd{m,n},-20*log10(abs(tlfd{m})),'g')
            plot(rd{m,n},-20*log10(tld_mfavg{m,n}),'b','LineWidth',1.5)
            plot(rd{m,n},-20*log10(tld_hhavg{m,n}),'m','LineWidth',1.5)
            xlabel('Range [m]')
            ylabel('-TL [dB]')
            set(gca,'XScale','log')
            axis([r(1) r(end) -140 -40])
            legend(sprintf('%0.0f Hz',fc),'Multi-Frequency','Harrison & Harrison (1995)')
            title({sprintf(['Transmission Loss In Frequency Band \\rm'...
                '(fc = %0.0f, bpo = %d, nr/\\alpha = %d, nf/\\alpha = %d, '...
                '\\sigma = %0.1f dB)'],fc,bpo(m),numRangesPerAlphaRes(n),...
                numFreqsPerAlphaRes,rmse(m,n)); sprintf(['hw = %0.0f, zs = %0.0f, '...
                'zr = %0.0f, Rb = %0.1f, Rs = %0.1f'],hw,zs,zr,Rb,Rs)}); 
            box on
            set(hfig,'units','normalize','outerposition',[0.15 0.1 0.7 0.8])

            % Save Figures
            figName = sprintf(['TL Average (fc = %0.0f, nralpha = %d, '...
                'nfalpha = %d)'],fc,numRangesPerAlphaRes(n),numFreqsPerAlphaRes);
            figPath = fullfile(outDir,fol1,fol2,figName);
            savefig(hfig,figPath,'compact')
            print(hfig,figPath,'-dpng','-r250')
            close(hfig) 
        end
    end
end

% Save Data
tlavgAnalysis.config.fname = fname;
tlavgAnalysis.config.fc = fc;
tlavgAnalysis.config.zs = zs;
tlavgAnalysis.config.zr = zr;
tlavgAnalysis.config.hw = hw;
tlavgAnalysis.config.Rb = Rb;
tlavgAnalysis.config.Rs = Rs;
tlavgAnalysis.bpo = bpo;
tlavgAnalysis.nralpha = numRangesPerAlphaRes;
tlavgAnalysis.nfalpha = numFreqsPerAlphaRes;
tlavgAnalysis.r = rd;
tlavgAnalysis.f = fd;
tlavgAnalysis.tlf = tlfd;
tlavgAnalysis.tlavg_mf = tld_mfavg;
tlavgAnalysis.tlavg_hh = tld_hhavg;
tlavgAnalysis.rmse = rmse;

save('tlavgAnalysis_r20kLog_f15','-struct','tlavgAnalysis')


