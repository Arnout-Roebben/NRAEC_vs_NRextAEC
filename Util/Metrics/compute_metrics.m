function res = compute_metrics(sig)
% Computes the metrics to evaluate the signals in sig.
%
% INPUT: 
% sig           Struct      Struct containing the following input signals:
% -m            TXM         M-microphone microphone signal of length T samples.
%                           m=s+n+es+en.
% -s            TXM         M-microphone desired speech signal of length T samples.
% -n            TXM         M-microphone noise signal of length T samples.
% -e            TXM         M-microphone echo signal of length T samples.
%
% OUTPUT:
% res           Struct      Struct containing the computed metrics:
% -snr          1X1         Signal to noise ratio.
% -ser          1X1         Signal to echo ratio.
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
res = struct(); % Struct to hold results

%% Calculate metrics
res.snr = SNR(sig.s, sig.n); % Signal-to-noise ratio
res.ser = SNR(sig.s, sig.e); % Signal-to-echo ratio                               

end