function X = WOLA_analysis(x,win,N,shift)
% Weighted overlap add (WOLA) analysis filterbank. Only the positive
% frequencies 0-fs/2 are returned.
% 
% INPUT:
% x         TXM     Vector in time domain of length T samples.
% win       NX1     Window.
% N         1X1     Discrete Fourier transform (DFT) size.
% shift     1X1     Frame shift.
%
% OUTPUT:
% X         MXKXN   Frequency matrix with K number of frames.
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
M = size(x,2); % Number of microphones
K = floor((length(x)-shift)/(N-shift)); % Number of frames K (See doc STFT)
X = nan(M,K,N/2+1); % Placeholder for the STFT-transformed result

%% Processing
% Convert to STFT domain
for l=1:K
    X_full = fft(x((l-1)*shift+1:(l-1)*shift+N,:).*repmat(win,1,M),N,1);
    X(:,l,:) = X_full(1:N/2+1,:).';
end