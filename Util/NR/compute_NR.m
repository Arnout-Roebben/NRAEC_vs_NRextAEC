function res = compute_NR(sig,p)
% Processing using a noise reduction (NR) using generalised eigenvalue
% decomposition (GEVD) approximation according to [Section 2.4,1] and [2]. The filters are 
% computed in the short-time Fourier transform (STFT) domain. Only one filter is 
% computed for the entire data. 
% 
% INPUT: 
% sig           Struct          Struct containing the following input signals:
% -m            TXM             M-microphone microphone signal of length T samples.
%                               m=s+n+es+en.
% -s            TXM             M-microphone desired speech signal of length T samples.
% -n            TXM             M-microphone near-end room noise signal of length T samples.
% -es           TXM             M-microphone far-end room speech component in the echo 
%                               signal of length T samples.
% -en           TXM             M-microphone far-end room noise component in the echo 
%                               signal of length T samples.
% -l            TXL             L-loudspeaker loudspeaker signal of length T samples. 
%                               l=ls+ln.
% -ls           TXL             L-loudspeaker far-end room speech component in the 
%                               loudspeaker signal of length T samples.
% -ln           TXL             L-loudspeaker far-end room noise component in the 
%                               loudspeaker signal of length T samples.
% p             Struct          Struct containing the following parameters:
% -ref          1X1             Reference microphone.
% -rank_s       1X1             [Optional] Rank to be used in the 'GEVD' procedure 
%                               for the desired speech correlation matrix. See
%                               compute_NR.m.
% -sensitivity  String          Sensitivity of the standard deviation in the 
%                               voice acitivity detector (VAD) formula, 
%                               see VAD.m.
% -fs           1X1             Sampling rate [Hz].
% -M            1X1             Number of microphones.
% -L            1X1             Number of loudspeakers.
% -N            1X1             Discrete Fourier transform (DFT) size. 
%                               See WOLA_analysis.m  and WOLA_synthesis.m
% -win          NX1             Window. See WOLA_analysis.m and WOLA_synthesis.m
% -shift        1X1             Frame shift. See WOLA_analysis.m and WOLA_synthesis.m
%
% OUTPUT:
% res           Struct          Struct containing the NR processed signals:
% -W            M+LXM+LXN/2+1   Short-time Fourier transform (STFT) NR filters
%                               to estimate the desired speech in each of the
%                               M+L microphones and loudspeakers 
%                               for the positive N/2+1 frequency bins.
% -m            TXM             See INPUT.
% -s            TXM             See INPUT.
% -n            TXM             See INPUT.
% -es           TXM             See INPUT.
% -en           TXM             See INPUT.
% -l            TXM             See INPUT.
% -ls           TXM             See INPUT.
% -ln           TXM             See INPUT.
%
% [1] A. Spriet, “Adaptive filtering techniques for noise reduction and acoustic
% feedback cancellation in hearing aids,” Ph.D. dissertation, KU Leuven,
% Leuven, Belgium, 2004.
% [2] R. Serizel, M. Moonen, B. Van Dijk, and J. Wouters, “Low-rank
% Approximation Based Multichannel Wiener Filter Algorithms for Noise
% Reduction with Application in Cochlear Implants,” IEEE/ACM Transactions
% on Audio, Speech, and Language Processing, vol. 22, no. 4, pp.
% 785–799, Apr. 2014.
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
% Parameters
m_f = WOLA_analysis(sig.m,p.win,p.N,p.shift); % STFT of microphone signal
% STFT of desired speech signal (for voice activity detector (VAD) purposes only)
s_f = WOLA_analysis(sig.s,p.win,p.N,p.shift); 
VADs = VAD(s_f,p.sensitivity,p.ref); % VAD speech
M = size(m_f,1); % Number of microphones
K = size(m_f,2); % Number of frames
N = size(m_f,3); % Number of frequency bins

% Preallocate memory
res = struct(); % Struct to hold results
res.W = nan(M,M,K,N); % NR filter

%% Processing
% Placeholder for the microphone signal whenever VADs(k,n)=1
dm = cell(N,1);
% Placeholder for the microphone signal whenever VADs(k,n)=0
dne = cell(N,1);
for k = 1:K % Loop over frames
    for n=1:N % Loop over bins
        % Store microphone signal whenever VADs(k,n)=1
        if VADs(k,n)
            dm{n}(:,(sum(VADs(1:k,n)==1))) = m_f(:,k,n);
        else
         % Store microphone signal whenever VADs(k,n)=0
            dne{n}(:,(sum(VADs(1:k,n)==0))) = m_f(:,k,n);
        end
    end
end

% Collect the correlation matrices
% Placeholder for the microphone correlation matrix whenever VADs(k,n)=1
Rmm = nan(M,M,N);
% Placeholder for the microphone correlation matrix whenever VADs(k,n)=0
Rne_ne = nan(M,M,N);
for n=1:N % Loop over bins
    % Compute the microphone correlation matrix whenever VADs(k,n)=1
    Rmm(:,:,n) = dm{n}*dm{n}'/(sum(VADs(:,n)==1));
    % Compute the microphone correlation matrix whenever VADs(k,n)=0
    Rne_ne(:,:,n) = dne{n}*dne{n}'/(sum(VADs(:,n)==0));
end

% Compute the NR filter
% Calculate the NR filter using the GEVD approximation
res.W = updateMWFGEVDMultichannel(Rmm,Rne_ne,repmat(p.rank_s,N,1)); 

end