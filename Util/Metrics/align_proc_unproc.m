function [proc,unproc] = align_proc_unproc(proc,unproc,p)
% Align processed and unprocessed signals, select data starting from a
% specified starting index, and only keep data when the desired speech is
% active.
% 
% INPUT:
% proc          Struct      Struct containing the following processed
%                           signals:
% -m            TXM         M-microphone microphone signal of length T samples.
%                           m=s+n+es+en.
% -s            TXM         M-microphone desired speech signal of length T samples.
% -n            TXM         M-microphone near-end room noise signal of length T samples.
% -es           TXM         M-microphone far-end room speech component in the echo 
%                           signal of length T samples.
% -en           TXM         M-microphone far-end room noise component in the echo 
%                           signal of length T samples.
% unproc          Struct    Struct containing the following unprocessed
%                           signals:
% -m            TXM         M-microphone microphone signal of length T samples.
%                           m=s+n+es+en.
% -s            TXM         M-microphone desired speech signal of length T samples.
% -n            TXM         M-microphone near-end room noise signal of length T samples.
% -es           TXM         M-microphone far-end room speech component in the echo 
%                           signal of length T samples.
% -en           TXM         M-microphone far-end room noise component in the echo 
%                           signal of length T samples.
% p             Struct      Struct containing the following parameters:
% -sensitivity  String      Sensitivity of the standard deviation in the 
%                           voice acitivity detector (VAD) formula, 
%                           see VAD.m.
% -ref          1X1         Reference microphone.
% -start        1X1         Sample index to serve as starting index for the
%                           processed and unprocessed signals.
% 
% OUTPUT:
% proc          TXM         Struct containing the aligned processed
% signals:
% -m            TXM         See INPUT.
% -n            TXM         See INPUT.
% -e            TXM         e=es+en. See INPUT.
% unproc        TXM         Struct containing the aligned unprocessed
% -m            TXM         See INPUT.
% -n            TXM         See INPUT.
% -e            TXM         e=es+en. See INPUT.
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

%% Align processed and unprocessed
% Select start index
start = p.N-1; % Delay due to distortion function

% Align processed and unprocessed
proc.m = proc.m(start:end,p.ref);
proc.s = proc.s(start:end,p.ref);
proc.n = proc.n(start:end,p.ref);
proc.e = proc.es(start:end,p.ref) + proc.en(start:end,p.ref);

%% Adjust start index
proc.m = proc.m(p.Tstart:end);
proc.s = proc.s(p.Tstart:end);
proc.n = proc.n(p.Tstart:end);
proc.e = proc.e(p.Tstart:end);

unproc.m = unproc.m(p.Tstart:end,p.ref);
unproc.s = unproc.s(p.Tstart:end,p.ref);
unproc.n = unproc.n(p.Tstart:end,p.ref);
unproc.e = unproc.es(p.Tstart:end,p.ref) + unproc.en(p.Tstart:end,p.ref);

%% Select timestamps speech is active
VADs = abs(proc.s(:,p.ref)) > std(proc.s(:,p.ref))*p.sensitivity;
proc.m = proc.m(VADs);
proc.s = proc.s(VADs);
proc.n = proc.n(VADs);
proc.e = proc.e(VADs);
unproc.m = unproc.m(VADs);
unproc.s = unproc.s(VADs);
unproc.n = unproc.n(VADs);
unproc.e = unproc.e(VADs);

end