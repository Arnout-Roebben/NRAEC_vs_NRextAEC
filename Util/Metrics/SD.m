function SD = SD(s,sp)
% Calculates the fullband signal distortion (SD) using the time domain
% signals.
%
% INPUT:
% s     TX1     Speech signal of length T samples.
% sp    TX1     Processed speech signal of length T samples.
% fs    1X1     Sampling frequency [Hz]
% 
% OUTPUT:
% SD    1X1     Fullband SD level.
%
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

%% Fullband SD
SD = 10*log10(mean(abs(s.^2))/mean(abs(sp.^2)));


end

