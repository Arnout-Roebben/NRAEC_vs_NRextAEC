function w_distortion = WOLA2distortion(W,win_analysis,win_synthesis,shift)
% Converts filters in the short-time Fourier transform (STFT) domain to their
% corresponding distortion functions according to [1].
%
% INPUT:
% W             N/2+1 X M   M STFT domain filters characterised by their
%                           coefficients corresponding to the N/2+1 
%                           non-negative frequencies. 
% win_analysis  N X 1       WOLA analysis window.
% win_synthesis N X 1       WOLA synthesis window.
% shift         1X1         WOLA frame shift [samples]. 
%
% OUTPUT:
% w_distortion  2*N-1 X M   M distortion functions corresponding to the M
%                           WOLA domain filters in W.
%
% [1] P. Didier, T. Van Waterschoot, S. Doclo and M. Moonen, 
%     "Sampling Rate Offset Estimation and Compensation for Distributed 
%     Adaptive Node-Specific Signal Estimation in Wireless Acoustic Sensor 
%     Networks," in IEEE Open Journal of Signal Processing, vol. 4, 
%     pp. 71-79, 2023, doi: 10.1109/OJSP.2023.3243851.
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
% Check the relation between the size of W and the size of win_analysis and
% win_synthesis
if size(W,1) ~= size(win_analysis,1)/2+1 || ...
   size(W,1) ~= size(win_synthesis,1)/2+1
    error(['The amount of rows in ''W'' must equal half the amount of ' ...
        'rows in ''win_analysis'' and  ''win_synthesis'' plus one!'])
end

N = size(win_analysis,1); % Number of frequency bins
M = size(W,2); % Number of channels
w_distortion = nan(2*N-1,M); % Allocate memory for the distortion function

%% Conversion
% Convert WOLA filter to time domain
W = conj(W); 
w = real(ifft(cat(1,W,conj(flip(W(2:end-1,:))))));

for m=1:M % Loop over channels
    % Construct circulant matrix out of time domain filter
    Hmat = gallery('circul',flip(w(:,m))).'; 

    % Multiply analysis and synthesis windows with circulant matrix
    Amat = diag(win_synthesis)*Hmat*diag(win_analysis);

    % Calculate the distortion function
    for l=1:2*N-1 % Loop over coefficients of the distortion function
        % Sum over corresponding diagonal
        w_distortion(l,m) = sum(diag(Amat,l-N)); 
    end
    
end

% Normalise the distortion function using the frame shift
w_distortion = w_distortion./shift;

end