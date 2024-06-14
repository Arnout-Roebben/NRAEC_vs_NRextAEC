function res = process_NRextAEC_adaptive(sig,p)
% Processing using an extended noise reduction (NRext) preceding an acoustic echo
% cancellation (AEC) (NR-AEC). The filters are computed adaptively.
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
% -rank_ses     1X1             [Optional] Rank to be used in the 'GEVD' procedure 
%                               for the sum of the extended desired speech 
%                               correlation matrix and the far-end room speech 
%                               component in the echo. See compute_NRext.m.
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
% -Lfhat        1X1             Number of coefficients in 
%                               Normalised least mean square (NLMS) estimated 
%                               AEC filter.
%                               See compute_AEC.m 
% -mu           1X1             NLMS stepsize. See compute_AEC.m 
% -alpha        1X1             NLMS regularisation factor. See compute_AEC.m 
%
% OUTPUT:         
% res           Struct          Struct containing the processed signals:
% -NR           Struct          Struct containing the processed signals after
%                               NR.
% --W           M+LXM+LXN/2+1XK Short-time Fourier transform (STFT) NR filters
%                               to estime the desired speech in each of the
%                               M+L microphones and loudspeakers 
%                               for the positive N/2+1 frequency bins for
%                               each of the K frames.
% --w           2N-1XM+LXM+LXK  Distortion equivalent NR filters of length 2N-1
%                               to estimate the desired speech in each of the M
%                               +L microphones and loudspeakers for each of
%                               the K frames.
% --m           TXM             See INPUT.
% --s           TXM             See INPUT.
% --n           TXM             See INPUT.
% --es          TXM             See INPUT.
% --en          TXM             See INPUT.
% --l           TXM             See INPUT.
% --ls          TXM             See INPUT.
% --ln          TXM             See INPUT.
% -AEC          Struct          Struct containing the processed signals after
%                               AEC.
% --fhat        LfhatXLXMXT     NLMS estimated AEC filter
%                               between each of the L loudspeakers and M microphones 
%                               of length Lfhat samples for each of the T
%                               input samples.
% --m           TXM             See INPUT.
% --s           TXM             See INPUT.
% --n           TXM             See INPUT.
% --es          TXM             See INPUT.
% --en          TXM             See INPUT.
% --l           TXM             See INPUT.
% --ls          TXM             See INPUT.
% --ln          TXM             See INPUT.
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
T = length(sig.m); % Number of samples in microphone signal 
K = floor((T-p.shift)/(p.N-p.shift)); % Number of frames K (See doc STFT)
update_frames = ((1:floor((T-p.shift)/(p.N-p.shift)))-1)*p.shift+p.N;

%% NRext
% Compute NRext filters in WOLA domain
res.NRext = compute_NRext_adaptive(sig,p);

% Convert NRext filters to time domain
res.NRext.w = nan(2*p.N-1,p.M+p.L,p.M+p.L,K);
for k=1:K % Loop across frames
    for m=1:p.M+p.L % Loop across microphones and loudspeakers
         % Compute time domain filter to retrieve the speech and speech echo
         res.NRext.w(:,:,m,k) = WOLA2distortion(permute(res.NRext.W(:,m,k,:),[4 1 2 3]), ...
             p.win,p.win,p.shift);
    end
end

% Initialise NR filter
w = zeros(2*p.N-1,p.M+p.L,p.M+p.L);
for m=1:p.M+p.L
   w(1,m,m) = 1;
end
% Apply NRext filters
for t=1:T % Loop across frames
    % Update filter
    if ismember(t,update_frames)
        w = res.NRext.w(:,:,:,t==update_frames);
    end   

    % Apply filter
    for m=1:p.M % Loop across microphones
        res.NRext.m(t,m) = sum(w(:,:,m).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.m(max(1,t-2*p.N+2):t,:) sig.l(max(1,t-2*p.N+2):t,:)]]),'all');
        res.NRext.s(t,m) = sum(w(:,:,m).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.s(max(1,t-2*p.N+2):t,:) zeros(size(sig.l(max(1,t-2*p.N+2):t,:)))]]),'all');
        res.NRext.n(t,m) = sum(w(:,:,m).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.n(max(1,t-2*p.N+2):t,:) zeros(size(sig.l(max(1,t-2*p.N+2):t,:)))]]),'all');
        res.NRext.es(t,m) = sum(w(:,:,m).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.es(max(1,t-2*p.N+2):t,:) sig.ls(max(1,t-2*p.N+2):t,:)]]),'all');
        res.NRext.en(t,m) = sum(w(:,:,m).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.en(max(1,t-2*p.N+2):t,:) sig.ln(max(1,t-2*p.N+2):t,:)]]),'all');
    end
    for l=1:p.L % Loop across loudspeakers
        res.NRext.l(t,l) = sum(w(:,:,p.M+l).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.m(max(1,t-2*p.N+2):t,:) sig.l(max(1,t-2*p.N+2):t,:)]]),'all');
        res.NRext.ls(t,l) = sum(w(:,:,p.M+l).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.es(max(1,t-2*p.N+2):t,:) sig.ls(max(1,t-2*p.N+2):t,:)]]),'all');
        res.NRext.ln(t,l) = sum(w(:,:,p.M+l).*flip([zeros(max(2*p.N-1-t,0),p.M+p.L);[sig.en(max(1,t-2*p.N+2):t,:) sig.ln(max(1,t-2*p.N+2):t,:)]]),'all');
    end
end

%% AEC
% Compute AEC filters
res.AEC = compute_AEC_adaptive(res.NRext,p);

% Apply AEC filters
for t=1:T
    for m=1:p.M
        res.AEC.m(t,m) = res.NRext.m(t,m) - sum(res.AEC.fhat(:,:,m,t).*flip([zeros(max(p.Lfhat-t,0),p.L);res.NRext.l(max(1,t-p.Lfhat+1):t,:)]),'all');
        res.AEC.es(t,m) = res.NRext.es(t,m) - sum(res.AEC.fhat(:,:,m,t).*flip([zeros(max(p.Lfhat-t,0),p.L);res.NRext.ls(max(1,t-p.Lfhat+1):t,:)]),'all');
        res.AEC.en(t,m) = res.NRext.en(t,m) - sum(res.AEC.fhat(:,:,m,t).*flip([zeros(max(p.Lfhat-t,0),p.L);res.NRext.ln(max(1,t-p.Lfhat+1):t,:)]),'all');
    end
end
% Speech, noise and loudspeaker signals are not affected by filter
res.AEC.s = res.NRext.s;
res.AEC.n = res.NRext.n;
res.AEC.l = res.NRext.l;
res.AEC.ls = res.NRext.ls;
res.AEC.ln = res.NRext.ln;

end