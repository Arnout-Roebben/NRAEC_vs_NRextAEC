% Compares the following filters for combined noise reduction (NR) and
% acoustic echo cancellation (AEC):
% *) NR-AEC: An NR preceding an AEC. The NR aims at supressing near-end
% room noise (and possibly the echo), while the AEC aims at supressing the
% echo.
% *) NRext-AEC: An extended NR (NRext) preceding an AEC. The NRext aims at
% supressing the near-end room noise and the far-end room noise component in 
% the echo, while the AEC aims at supressing the far-end room speech component
% in the echo (and the residual far-end room noise component in the echo).
%
% v1.0
% LICENSE: This software is distributed under the terms of the MIT license (See LICENSE.md).
% AUTHOR:  Arnout Roebben
% CONTACT: arnout.roebben@esat.kuleuven.be
% CITE: A. Roebben, T. van Waterschoot, and M. Moonen, "Cascaded noise 
% reduction and acoustic echo cancellation based on an extended noise 
% reduction," Accepted for publication in EUSIPCO, Lyon, France, Aug. 2024.
% and
% A. Roebben, “Github repository: Cascaded noise reduction and acoustic echo 
% cancellation based on an extended noise reduction,”
% https://github.com/Arnout-Roebben/NRAEC_vs_NRextAEC, 2024.
%
% A preprint is available at
% A. Roebben, T. van Waterschoot, and M. Moonen, "Cascaded noise 
% reduction and acoustic echo cancellation based on an extended noise 
% reduction," 2024, arXiv:2406.08974.

%% Check that current folder corresponds to NRAEC_vs_NRextAEC
[~,curr_dir] = fileparts(pwd);
if ~strcmp(curr_dir,"NRAEC_vs_NRextAEC")
    error('Current fulder must correspond to ''NRAEC_vs_NRextAEC''!')
end

%% Cleanup
clear; close all; clc;
rng(2,"twister"); % Fix random number generator
addpath(genpath('.')); % Add subfolders of directory to path 

%% Load audio
% Audio data consists of the following components
% -m    TXM     M-microphone microphone signal of length T samples.
%               m=s+n+es+en.
% -s    TXM     M-microphone desired speech signal of length T samples.
% -n    TXM     M-microphone near-end room noise signal of length T samples.
% -es   TXM     M-microphone far-end room speech component in the echo 
%               signal of length T samples.
% -en   TXM     M-microphone far-end room noise component in the echo 
%               signal of length T samples.
% -l    TXL     L-loudspeaker loudspeaker signal of length T samples. 
%               l=ls+ln.
% -ls   TXL     L-loudspeaker far-end room speech component in the 
%               loudspeaker signal of length T samples.
% -ln   TXL     L-loudspeaker far-end room noise component in the 
%               loudspeaker signal of length T samples.
% For the desired speech and far-end room speech component in the
% loudspeakers, sentences from the Voice Cloning Toolkit (VCTK) corpus were
% used [1], as also available at [2].
% [1] C. Veaux, J. Yamagishi, K. MacDonald, "CSTR VCTK Corpus:
% English Multi-speaker Corpus for CSTR Voice Cloning Toolkit," http:
% //homepages.inf.ed.ac.uk/jyamagis/page3/page58/page58.html, 2016.
% [2] Dietzen, T., Ali, R., Taseska, M., van Waterschoot, T.: "Data
% Repository for MYRiAD: A Multi-Array Room Acoustic Database,"
% https://zenodo.org/record/7389996, 2023.
load('.\Audio\sig.mat');
% LX1 cell array containing the LfXM impulse responses from loudspeaker l={1,...,L}
% to M microphones of Lf coefficients.
load('.\Audio\impulse.mat');

%% Processing parameters
p = struct(); % Struct containing processing parameters

% General 
p.fs = fs; % Sampling rate
p.ref = 1; % Reference microphone
p.M = size(sig.m,2); % Amount of microphones
p.L = size(sig.l,2); % Amount of loudspeakers
% If 1 the filters are adaptively updated. If 0 the filters are computed
% across the entire data, 
p.adaptive = 0; 

% Frequency transform (See also WOLA_analysis.m and WOLA_synthesis.m)
p.N = 512; % Discrete Fourier transform (DFT) size N
p.shift = p.N/2; % Frame shift for weighted overlap add (WOLA)
p.win = sqrt(hann(p.N,'periodic')); % Window 

% Voice activity detection (VAD) (See also VAD.m)
p.sensitivity = 1e-5; % Sensitivity of VAD

% Noise reduction (NR) parameters (See also compute_NR.m and compute_NR_adaptive.m)
p.rank_s = 1; % Requested rank of desired speech correlation matrix 
p.lambda = 0.995; % Smoothing factor when updating the matrices (only when p.adaptive==1)

% Extended noise reduction (NRext) parameters (See also compute_NRext.m and compute_NRext_adaptive.m)
% Requested rank of the sum of the extended desired speech  
% and far-end room speech component in the echo correlation matrix.
p.rank_ses = p.L + 1; 

% Acoustic echo cancellation (AEC) parameters (See also compute_AEC.m and compute_AEC_adaptive.m)
p.Lfhat = 128; % Length of the time-domain AEC filters
p.mu = 0.1; % Normalised least mean squares (NLMS) stepsize
p.alpha = 1e-6; % NLMS regularisation factor

% Metrics (See align_proc_unproc.m)
% Start index [sample], after which the data is used to compute the metrics.
p.Tstart = 1; 

%% Process
if p.adaptive
    NRAEC = process_NRAEC_adaptive(sig,p);    
    NRextAEC = process_NRextAEC_adaptive(sig,p);    
else
    NRAEC = process_NRAEC(sig,p);    
    NRextAEC = process_NRextAEC(sig,p);
end

%% Metrics
% Align processed and unprocessed
[proc_NRAEC,unproc_NRAEC] = align_proc_unproc(NRAEC.AEC,sig,p);
[proc_NRextAEC,unproc_NRextAEC] = align_proc_unproc(NRextAEC.AEC,sig,p);

% Calculate metrics
% NR-AEC
NRAEC.metrics = compute_metrics(proc_NRAEC);
NRAEC.metrics.sd = SD(unproc_NRAEC.s,proc_NRAEC.s);
% NRext-AEC
NRextAEC.metrics = compute_metrics(proc_NRextAEC);
NRextAEC.metrics.sd = SD(unproc_NRextAEC.s,proc_NRextAEC.s);
% Reference
metrics_ref = compute_metrics(unproc_NRAEC);

%% Visualisation
% Metrics
fprintf('NR-AEC:\n')
fprintf('\t SNR improvement: %f\n',NRAEC.metrics.snr - metrics_ref.snr);
fprintf('\t SER improvement: %f\n',NRAEC.metrics.ser - metrics_ref.ser);
fprintf('\t SD: %f\n\n',NRAEC.metrics.sd);

fprintf('NRext-AEC:\n')
fprintf('\t SNR improvement: %f\n',NRextAEC.metrics.snr - metrics_ref.snr);
fprintf('\t SER improvement: %f\n',NRextAEC.metrics.ser - metrics_ref.ser);
fprintf('\t SD: %f\n',NRextAEC.metrics.sd);

% Signals
figure; hold on
t = tiledlayout(2,2);
ax = nexttile; hold on; plot(sig.s(:,p.ref)); plot(NRAEC.AEC.s(:,p.ref)); 
title('Desired speech'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.n(:,p.ref)); plot(NRAEC.AEC.n(:,p.ref)); 
title('Noise'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.es(:,p.ref)+sig.en(:,p.ref)); 
plot(NRAEC.AEC.es(:,p.ref)+NRAEC.AEC.en(:,p.ref)); 
title('Echo'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(impulse{1}(:,p.ref),'color',[0.4940 0.1840 0.5560]); 
plot(NRAEC.AEC.fhat(:,1,p.ref,end),'color',[0.4660 0.6740 0.1880]); 
title('Echo path impulse response'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
legend('True','AEC filter')
title(t,'NR-AEC');
lg  = legend(ax,["Input" "Output"],'Orientation','Horizontal'); 
lg.Layout.Tile = 'South'; hold off;

figure; hold on;
t = tiledlayout(2,2);
ax = nexttile; hold on; plot(sig.s(:,p.ref)); plot(NRextAEC.AEC.s(:,p.ref)); 
title('Desired speech'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.n(:,p.ref)); plot(NRextAEC.AEC.n(:,p.ref)); 
title('Noise'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(sig.es(:,p.ref)+sig.en(:,p.ref)); 
plot(NRextAEC.AEC.es(:,p.ref)+NRextAEC.AEC.en(:,p.ref)); 
title('Echo'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
nexttile; hold on; plot(impulse{1}(:,p.ref),'color',[0.4940 0.1840 0.5560]); 
plot(NRextAEC.AEC.fhat(:,1,p.ref,end),'color',[0.4660 0.6740 0.1880]); 
legend('True','AEC filter')
title('Echo path impulse response'); xlabel('Time [Samples]'); ylabel('Amplitude [Arb. unit]');
title(t,'NRext-AEC');
lg  = legend(ax,["Input" "Output"],'Orientation','Horizontal'); 
lg.Layout.Tile = 'South'; hold off;