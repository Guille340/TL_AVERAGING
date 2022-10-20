%  tlb = TLAVG(r,tlf,f,bw)
%
%  DESCRIPTION: calculates the multi-frequency average over a bandwidth BW for
%  a single-frequency, range-dependent transmission loss curve TLF using the
%  range-averaging method from Harrison & Harrison (1995). 
%
%  In their paper, Harrison & Harrison (H&H) propose a method that further 
%  simplifies the calculation of the transmission loss within a frequency band, 
%  eliminating the need for a large number of simulations TL(r) curves at 
%  closely spaced frequencies within the band. The method takes advantage of 
%  the equivalence between spatial and frequency averaging to provide a fast 
%  calculation of the transmission loss TL(r) within a band centred at a 
%  frequency F using only transmission loss TL(r) at that single frequency 
%  (TLF). In particular, the solution applies a range-moving average to the
%  relative intensiy curve (I = 1/TLF^2) using a Gaussian window of width 
%  proportional to range.
%
%  The function uses a matrix kernel to speed up the processing. To avoid
%  an excessive use of memory, TLF is split into small sections to limit
%  the size of the kernel. With this approach, the algorithm can handle any 
%  vector length for TLF, at the time it optimises the speed.
%  
%  INPUT VARIABLES
%  - r: vector of ranges for TLF [m]
%  - tlf: pressure transmission loss vector at R ranges. The TLF is unitless.
%    If the transmission loss level tl_dB is available instead, the pressure 
%    transmission loss can be easily calculated as TLF = 10.^(tl_dB/20).
%  - f: frequency at which TLF was calculated [Hz]
%  - bw: bandwidth of the band over which to calculate the frequency-averaged 
%    approximation TLB of the single-frequency transmission loss TLF.
%
%  OUTPUT VARIABLES
%  - tlb: frequency-average approximation of the single-frequency transmission
%    loss TLF.
%
%  INTERNALLY CALLED FUNCTIONS
%  - None
%
%  CONSIDERATIONS & LIMITATIONS
%  - This is a very powerful method for simplifying transmission loss
%    calculations within a frequency band, but it has its limitations. After
%    an in-depth analysis of its performance I found out that the accuracy
%    of the approximation is closely related to the range resolution and
%    the characteristics of the averaging window. The main conclusions I 
%    arrived at for evenly-spaced range points (relevant for AcTUP and
%    most modelling scenarios) were:
%
%    1. The range resolution step dr must be dr <= rmin/(37.5*bpo) to get
%    an error at range of less than 0.3 dB, where rmin is the minimum range
%    that meets the error condition, and bpo is the number of bands per
%    octave (i.e. 3 indicates an equivalent frequency average over a
%    third-octave band).
%
%    2. An "overfit" approximation may occur when the relative averaging band
%    is too narrow. To avoid overfitting, the number of bands per octave must
%    be bpo < fc * rmin * 1e-5. A larger bpo is preferred to a small
%    value that may result in "oversmoothing" (see next point).
%
%    3. An "oversmoothed" approximation will occur when the relative averaging 
%    band is too wide. To avoid oversmoothing, the number of bands per octave 
%    must be bpo >= 3 to bring the standard error down to 0.5 dB or less.
%
%  REFERENCES
%  - Harrison, C.H, Harrison, J.A. (1995). "A simple relationship between
%    frequency and range averages for broadband sonar", J. Acoust. Soc. Am., 
%    97(2), 1314-1317.
%
%  See also tlavg_ex1.m, tlavg_ex2.m, lloydMirror.m, imageSourceModel.m
%  

%  REVISION 1.1
%  - Added function help
%
%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  12 May 2020

function tlb = tlavg(r,tlf,f,bw)
    
    % Error Control
    if ~isnumeric(r) || ~isnumeric(tlf) || ~isvector(r) || ~isvector(tlf) ...
            || ~isequal(length(r),length(tlf))
        error(['Input arguments R and TLF must be numeric vectors of '...
            'the same length'])
    end
    
    if r <= 0
        warning(['One or more values in R are negative or zero and '...
            'will be ignored'])
        ipos = r > 0;
        r = r(ipos);
        tlf = tlf(ipos);
    end
    
    if ~isnumeric(f) || f <= 0 || length(f) > 1
        error('Input argument F must be a positive single number')
    end
    
    if ~isnumeric(bw) || bw < 0 || length(bw) > 1
        error('Input argument F must be a single number >= 0')
    end
    
    tlfAbs = abs(tlf);
    if ~all(tlfAbs >= 1)
        if all(tlfAbs > 0 & tlfAbs <= 1)
            warning(['TLF must contain complex/real sound pressure '...
                'transmission loss (TL) values (ABS(TL) >= 1). Possible '...
                'inverse input (TLF = 1/TL)'])
        elseif all(tlfAbs <= 0) || all(tlfAbs >= 0)
            warning(['TLF must contain complex/real sound pressure '...
                'transmission loss (TL) values (ABS(TL) >= 1). Possible '...
                'input decibels (TLF = +/- 20log10(TL)'])
        end
    end
    
    % Dimensions
    tlfSize = size(tlf);
    r = r(:);
    tlf = tlf(:);
    
    % Remove NaN Values
    ival = ~isnan(r) & ~isnan(tlf);
    r = r(ival);
    tlf = tlf(ival);
            
    % Calculate Bytes of Averaging Kernel Matrix
    if isa(r,'float')
        nBytes = 8; % double precission initial guess
        if isa(r,'single')
            nBytes = 4;
        end
    end
    bytesPerMbyte = 1024^2;
    maxMbytes = 50; 
    nRanges = length(r);
    nCols = min([nRanges, ceil(maxMbytes*bytesPerMbyte/(nBytes*nRanges))]);
    nBlocks = ceil(nRanges/nCols);
    
    % Process Range Average in Blocks
    J = abs(1./tlf(:)).^2; % relative intensity at frequency F
    alpha = bw/(2*sqrt(log(2))*f); % fractional bandwidth (-3dB, half power)
    Jb = J; % initialise relative intensity in frequency band
    if alpha > 0 % if the bandwidth is not zero (bw = 0 -> single frequency)
        for m = 1:nBlocks-1
            i1 = nCols*(m-1) + 1; % start index
            i2 = nCols*m; % ending index
            r0 = r(i1:i2)'; % vector of central ranges for gaussian masks            
            gaussKernel = exp(-((r-r0)./(alpha*r0)).^2); % gaussian masks
            Jb(i1:i2) = sum(J .* gaussKernel)./sum(gaussKernel);
        end
        i1 = nCols*(nBlocks-1) + 1;
        i2 = nRanges;
        r0 = r(i1:i2)';
        gaussKernel = exp(-((r-r0)./(alpha*r0)).^2);
        Jb(i1:i2) = sum(J .* gaussKernel)./sum(gaussKernel);
    end
    
    tlb = nan(tlfSize); % initialise TLB with NaN values
    tlb(ival) = sqrt(1./Jb); % transmission loss in band (magnitude)
    
    
    