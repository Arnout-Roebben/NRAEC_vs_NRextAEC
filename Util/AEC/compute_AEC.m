function res = compute_AEC(sig,p)
% Processing using an acoustic echo cancellation (AEC) using a normalised
% least mean square (NLMS) procedure according to [Section 2.1,1]. Only one 
% filter is computed for the  entire data.
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
% -Lfhat        1X1         Number of coefficients in 
%                           Normalised least mean square (NLMS) estimated 
%                           AEC filter.
% -mu           1X1         NLMS stepsize. 
% -alpha        1X1         NLMS regularisation factor. 
%
% OUTPUT:
% res           Struct      Struct containing the processed signals after
%                           AEC.
% -fhat         LfhatXLXM   NLMS estimated AEC filter
%                           between each of the L loudspeakers and M microphones 
%                           of length Lfhat samples.
% -m            TXM         See INPUT.
% -s            TXM         See INPUT.
% -n            TXM         See INPUT.
% -es           TXM         See INPUT.
% -en           TXM         See INPUT.
% -l            TXM         See INPUT.
% -ls           TXM         See INPUT.
% -ln           TXM         See INPUT.
%
% [1] G. Rombouts, "Adaptive filtering algorithms for acoustic echo and noise
% cancellation," Ph.D. dissertation, KU Leuven, Leuven, Belgium, 2003.
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
M = size(sig.m,2); % Number of microphones
L = size(sig.l,2); % Number of loudspeakers
T = size(sig.m,1); % Number of samples

%% Processing
% Preallocate memory
res = struct(); % Struct containing results
res.fhat = nan(p.Lfhat,L,M); % Echo path impulse response
fhat = zeros(p.Lfhat*L,1); % Inital echo path impulse response

% VAD
VAD = abs(squeeze(sig.s)) > std(sig.s)*p.sensitivity;

% Calculate echo impulse responses
for m=1:M % Loop across microphones
    for t=1:T % Loop across samples
        % Update correlation matrices
        if t>p.Lfhat
            if ~VAD(t,m)
                lbar = reshape(flip(sig.l(t-p.Lfhat+1:t,:)),[],1);
                fhat = fhat + p.mu/(p.alpha+lbar'*lbar)*lbar*(sig.m(t,m)-fhat'*lbar);
            end
        end
    end
    % Compute echo impulse response
    res.fhat(:,:,m) = reshape(fhat,[p.Lfhat,L]);
    res.fhat(:,:,m) = res.fhat(:,:,m) - mean(res.fhat(:,:,m),1);
end
end