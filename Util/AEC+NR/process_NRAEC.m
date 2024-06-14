function res = process_NRAEC(sig,p)
% Processing using a noise reduction (NR) preceding an acoustic echo
% cancellation (AEC) (NR-AEC). The filters are computed across the entire
% data before applying the filters to the data.
%
% INPUT: 
% sig           Struct      Struct containing the following input signals:
% -m            TXM         M-microphone microphone signal of length T samples.
%                           m=s+n+es+en.
% -s            TXM         M-microphone desired speech signal of length T samples.
% -n            TXM         M-microphone near-end room noise signal of length T samples.
% -es           TXM         M-microphone far-end room speech component in the echo 
%                           signal of length T samples.
% -en           TXM         M-microphone far-end room noise component in the echo 
%                           signal of length T samples.
% -l            TXL         L-loudspeaker loudspeaker signal of length T samples. 
%                           l=ls+ln.
% -ls           TXL         L-loudspeaker far-end room speech component in the 
%                           loudspeaker signal of length T samples.
% -ln           TXL         L-loudspeaker far-end room noise component in the 
%                           loudspeaker signal of length T samples.
% p             Struct      Struct containing the following parameters:
% -ref          1X1         Reference microphone.
% -rank_s       1X1         [Optional] Rank to be used in the 'GEVD' procedure 
%                           for the desired speech correlation matrix. See
%                           compute_NR.m.
% -sensitivity  String      Sensitivity of the standard deviation in the 
%                           voice acitivity detector (VAD) formula, 
%                           see VAD.m.
% -fs           1X1         Sampling rate [Hz].
% -M            1X1         Number of microphones.
% -L            1X1         Number of loudspeakers.
% -N            1X1         Discrete Fourier transform (DFT) size. 
%                           See WOLA_analysis.m  and WOLA_synthesis.m
% -win          NX1         Window. See WOLA_analysis.m and WOLA_synthesis.m
% -shift        1X1         Frame shift. See WOLA_analysis.m and WOLA_synthesis.m
% -Lfhat        1X1         Number of coefficients in 
%                           Normalised least mean square (NLMS) estimated 
%                           AEC filter.
%                           See compute_AEC.m 
% -mu           1X1         NLMS stepsize. See compute_AEC.m 
% -alpha        1X1         NLMS regularisation factor. See compute_AEC.m 
%
% OUTPUT:         
% res           Struct      Struct containing the processed signals:
% -NR           Struct      Struct containing the processed signals after
%                           NR.
% --W           MXMXN/2+1   Short-time Fourier transform (STFT) NR filters
%                           to estime the desired speech in each of the M
%                           microphones for the positive N/2+1 frequency bins.
% --w           2N-1XMXM    Distortion equivalent NR filters of length 2N-1
%                           to estimate the desired speech in each of the M
%                           microphones.
% --m           TXM         See INPUT.
% --s           TXM         See INPUT.
% --n           TXM         See INPUT.
% --es          TXM         See INPUT.
% --en          TXM         See INPUT.
% --l           TXM         See INPUT.
% --ls          TXM         See INPUT.
% --ln          TXM         See INPUT.
% -AEC          Struct      Struct containing the processed signals after
%                           AEC.
% --fhat        LfhatXLXM   NLMS estimated AEC filter
%                           between each of the L loudspeakers and M microphones 
%                           of length Lfhat samples.
% --m           TXM         See INPUT.
% --s           TXM         See INPUT.
% --n           TXM         See INPUT.
% --es          TXM         See INPUT.
% --en          TXM         See INPUT.
% --l           TXM         See INPUT.
% --ls          TXM         See INPUT.
% --ln          TXM         See INPUT.
%
% v1.0
% LICENSE: This software is distributed under the terms of the MIT license (See LICENSE.md).
% AUTHOR:  Arnout Roebben
% CONTACT: arnout.roebben@esat.kuleuven.be
% CITE: A. Roebben, T. van Waterschoot, and M. Moonen, "Cascaded noise 
% reduction and acoustic echo cancellation based on an extended noise 
% reduction," in EUSIPCO, Lyon, France, Aug. 2024, pp. .
% and
% A. Roebben, “Github repository: Cascaded noise reduction and acoustic echo 
% cancellation based on an extended noise reduction,”
% https://github.com/Arnout-Roebben/NRAEC_vs_NRextAEC, 2024.
%
% A preprint is available at
% A. Roebben, T. van Waterschoot, and M. Moonen, "Cascaded noise 
% reduction and acoustic echo cancellation based on an extended noise 
% reduction," 2024, arXiv:2406.08974.

%% Initialisation
res = struct(); % Struct containing results

%% NR
% Compute NR filters in WOLA domain
res.NR = compute_NR(sig,p);

% Convert NR filters to time domain
res.NR.w = nan(2*p.N-1,p.M,p.M);
for m=1:p.M % Loop across microphones
     % Compute time domain filter to retrieve the speech in microphone m
     res.NR.w(:,:,m) = WOLA2distortion(permute(res.NR.W(:,m,:),[3 1 2]), ...
         p.win,p.win,p.shift);
end

% Apply NR filters
for m=1:p.M % Loop across microphones
    res.NR.m(:,m) = sum(fftfilt(res.NR.w(:,:,m),sig.m),2);
    res.NR.s(:,m) = sum(fftfilt(res.NR.w(:,:,m),sig.s),2);
    res.NR.n(:,m) = sum(fftfilt(res.NR.w(:,:,m),sig.n),2);
    res.NR.es(:,m) = sum(fftfilt(res.NR.w(:,:,m),sig.es),2);
    res.NR.en(:,m) = sum(fftfilt(res.NR.w(:,:,m),sig.en),2);
end
% Loudspeaker signals are not affected by filter
res.NR.l = sig.l;
res.NR.ls = sig.ls;
res.NR.ln = sig.ln;

%% AEC
% Compute AEC filters
res.AEC = compute_AEC(res.NR,p);

% Apply AEC filters
for m=1:p.M
    res.AEC.m(:,m) = res.NR.m(:,m) - sum(fftfilt(res.AEC.fhat(:,:,m),res.NR.l),2);
    res.AEC.es(:,m) = res.NR.es(:,m) - sum(fftfilt(res.AEC.fhat(:,:,m),res.NR.ls),2);
    res.AEC.en(:,m) = res.NR.en(:,m) - sum(fftfilt(res.AEC.fhat(:,:,m),res.NR.ln),2);
end
% Speech, noise and loudspeaker signals are not affected by filter
res.AEC.s = res.NR.s;
res.AEC.n = res.NR.n;
res.AEC.l = res.NR.l;
res.AEC.ls = res.NR.ls;
res.AEC.ln = res.NR.ln;

end