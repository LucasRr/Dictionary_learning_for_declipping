function [y, ClippingLevel] = clip_signal(x, SNR_target)
%  Clip a signal at a given SNR level, using bisection method.
%  This method should be precise up to +/- 0.001 dB
%  
%  Input:
%         - x: clean input signal
%         - SNR_target: target input SNR
%  Output: 
%         - y: signal clipped at SNR_target +/- 0.001 dB
%         - ClippingLevel: clipping level
% ------------------
%
% Author: Lucas Rencker
% Last update: 28/03/18


ClippingLevel1 = 0;
ClippingLevel2 = max(abs(x));
SNRtmp = inf;
it = 0;

% Search between ClippingLevel1 and ClippingLevel2:

while abs(SNRtmp-SNR_target) > 0.001 && it < 20
    
    it = it+1;
    
    ClippingLevel = (ClippingLevel1+ClippingLevel2)/2;
    y = max(min(x,ClippingLevel),-ClippingLevel); % clip signal
    SNRtmp = SNR(x,y); % check SNR
        
    % update search interval
    if SNRtmp < SNR_target
        ClippingLevel1 = ClippingLevel;
    else
        ClippingLevel2 = ClippingLevel;
    end        
    
end
    