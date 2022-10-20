%  DESCRIPTION: example showing a comparison between multi-frequency averaging
%  and Harrison & Harrison (1995) range averaging for the calculation of the 
%  transmission loss in a frequency band. This example plots the results from 
%  the two methods and the single frequency simulation in a semi-infinite water
%  medium (Lloyd Mirror scenario).
%
%  The script also creates and saves a 'tlAvgData_...' structure containing 
%  the transmission loss curves for all simulated frequencies. This structure
%  can be loaded with functions tlavg_analysisLog.m or tlavg_analysisLin.m
%  to plot and save the comparison of the original single-frequency 
%  transmission loss curve with those generated with the multi-frequency 
%  and the Harrison & Harrison (1995) averaging methods.
%
%  See also lloydMirror.m, tlavg.m, tlavg_ex1.m, tlavg_analysisLog.m

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  12 May 2020

% 1) Input Parameters
fc = 5000; % central frequency [Hz]
nFreqsPerOct = 500; % number of frequencies in octave
bwFact = 0.019254; % bandwidth factor (bwFact = bw/fc, 0.231 for third-octave)
bwLevel = -3; % relative level at which bandwidth is calculated [dB]
nRanges = 10000; % number of ranges
rangeScale = 'lin'; % scale of range vector ('lin' or 'log')
rLim = [1.5e1 30e3]; % start and end range [m]
zs = 140; % source depth [m]
zr = 150; % receiver depth [m]

% 2) Range and Frequency Vectors
if strcmp(rangeScale,'lin')
    r = rLim(1):(rLim(2)-rLim(1))/(nRanges-1):rLim(2); % range vector (linear)
else % if RANGESCALE = 'log' or other
    rStep = (rLim(2)/rLim(1)).^(1/(nRanges-1));
    r = [rLim(1) rLim(1)*rStep.^(1:nRanges-1)];
end
f1 = fc*2^(-3/(2*36));
f2 = fc*2^(3/(2*36));
f = f1:(f2-f1)/(nFreqsPerOct-1):f2;

% 3) In-Band Transmission Loss from Multi-Frequency Simulation
bw = bwFact*fc;
bwEnergy = 10^(bwLevel/10);
alpha = bw/(2*sqrt(-log(bwEnergy))*fc);
tl = zeros(nRanges,nFreqsPerOct);
tic
for m = 1:nFreqsPerOct
    fprintf('Frequency %d/%d (%s)\n',m,nFreqsPerOct,datestr(toc/86400,'HH:MM:SS'))
    tl(:,m) = lloydMirror(r,zr,zs,f(m));
end
J1 = abs(1./tl).^2;
gaussKernel = exp(-((f-fc)./(0.88*alpha*fc)).^2);
Jb1 = sum(J1 .* gaussKernel,2)/sum(gaussKernel);
tlAvg_multiFreq = sqrt(1./Jb1);

% 4) In-Band Transmission Loss from Harrison & Harrison (1995) Average
tlf = lloydMirror(r,zr,zs,fc);
tlAvg_harrison = tlavg(r,tlf,fc,bw);

% 5) Compare Results
figure
hold on
plot(r,-20*log10(abs(tlf)),'g','LineWidth',0.5)
plot(r,-20*log10(abs(tlAvg_multiFreq)),'b','LineWidth',1.5)
plot(r,-20*log10(abs(tlAvg_harrison)),'m','LineWidth',1.5)
xlim(rLim)
xlabel('Range [m]')
ylabel('-TL [dB]')
set(gca,'XScale','log')
legend(sprintf('%0.0f Hz',fc),'Multi-Frequency','Harrison & Harrison (1995)')
title({sprintf(['Transmission Loss In Frequency Band \\rm(BwFact = %0.3f, '...
    'nFreqs = %d)'],bwFact,nFreqsPerOct);sprintf('zs = %0.0f, zr = %0.0f',...
    zs,zr)}); 
box on

% % 6) Save
% tlavgData.config.fc = fc;
% tlavgData.config.nFreqsPerOct = nFreqsPerOct;
% tlavgData.config.nRanges = nRanges;
% tlavgData.config.rangeScale = rangeScale;
% tlavgData.config.rLim = rLim;
% tlavgData.config.zs = zs;
% tlavgData.config.zr = zr;
% tlavgData.config.hw = hw;
% tlavgData.config.Rb = Rb;
% tlavgData.config.Rs = Rs;
% tlavgData.f = f;
% tlavgData.r = r;
% tlavgData.tlf = tlf;
% tlavgData.tl = tl;
% 
% save('tlavgData_r20kLog_f15','-struct','tlavgData')