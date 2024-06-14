function res = compute_NRext_adaptive(sig,p)
% Processing using an extended noise reduction (NRext) using generalised eigenvalue
% decomposition (GEVD) approximation. The filters are computed in the
% short-time Fourier transform (STFT) domain. The filters are computed adaptively.
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
%
% OUTPUT:
% res           Struct      Struct containing the NR processed signals:
% -W            MXMXN/2+1XK Short-time Fourier transform (STFT) NR filters
%                           to estime the desired speech in each of the M
%                           microphones for the positive N/2+1 frequency
%                           bins and for each of the K frames.
% -m            TXM         See INPUT.
% -s            TXM         See INPUT.
% -n            TXM         See INPUT.
% -es           TXM         See INPUT.
% -en           TXM         See INPUT.
% -l            TXM         See INPUT.
% -ls           TXM         See INPUT.
% -ln           TXM         See INPUT.
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
l_f = WOLA_analysis(sig.l,p.win,p.N,p.shift); % STFT of loudspeaker signal
% STFT of desired speech signal (for voice activity detector (VAD) purposes only)
s_f = WOLA_analysis(sig.s,p.win,p.N,p.shift); 
% STFT of far-end room speech component in the loudspeaker signal 
% (for voice activity detector (VAD) purposes only)
es_f = WOLA_analysis(sig.es,p.win,p.N,p.shift); 
VADs = VAD(s_f,p.sensitivity,p.ref); % VAD speech
VADes = VAD(es_f,p.sensitivity,p.ref); % VAD speech echo
M = size(m_f,1); % Number of microphones
L = size(l_f,1); % Number of loudspeakers
K = size(m_f,2); % Number of frames
N = size(m_f,3); % Number of frequency bins

% Preallocate memory
res = struct(); % Struct to hold results
res.W = nan(M+L,M+L,K,N); % NRext filter

%% Processing
% Collect the correlation matrices
% Placeholder for the extended microphone correlation matrix whenever 
% VADs(k,n)=1 and VADes(k,n)=1
Rmm = nan(M+L,M+L,N,K+1);
Rmm(:,:,:,1) = repmat(1e-6*eye(M+L),[1 1 N]);
% Placeholder for the extended microphone correlation matrix whenever 
% VADs(k,n)=0 and VADes(k,n)=0
Rnen_nen = nan(M+L,M+L,N,K+1);
Rnen_nen(:,:,:,1) = repmat(1e-6*eye(M+L),[1 1 N]);
for k = 1:K % Loop over frames
    for n=1:N % Loop over bins
        % Store extended microphone signal whenever VADs(k,n)=1 and VADes(k,n)=1
        if VADs(k,n) && VADes(k,n)
            Rmm(:,:,n,k+1) = p.lambda*Rmm(:,:,n,k) +  ...
                (1-p.lambda)*cat(1,squeeze(m_f(:,k,n)),permute(l_f(:,k,n),[1 3 2]))*...
                cat(1,squeeze(m_f(:,k,n)),permute(l_f(:,k,n),[1 3 2]))';
        else
            Rmm(:,:,n,k+1) = Rmm(:,:,n,k);
        end
        % Store extended microphone signal whenever VADs(k,n)=0 and VADes(k,n)=0
        if ~VADs(k,n) && ~VADes(k,n)
            Rnen_nen(:,:,n,k+1) = p.lambda*Rnen_nen(:,:,n,k) + ...
                (1-p.lambda)*cat(1,squeeze(m_f(:,k,n)),permute(l_f(:,k,n),[1 3 2]))*...
                cat(1,squeeze(m_f(:,k,n)),permute(l_f(:,k,n),[1 3 2]))';  
        else
            Rnen_nen(:,:,n,k+1) = Rnen_nen(:,:,n,k);
        end
    end

    % Compute the NRext filter
    res.W(:,:,k,:) = updateMWFGEVDMultichannel(Rmm(:,:,:,k+1),Rnen_nen(:,:,:,k+1),repmat(p.rank_ses,N,1));
end
end