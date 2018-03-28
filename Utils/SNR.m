function snr = SNR(x_ref,x)
% Signal-to-noise Ratio
%
% Inputs:
%          - x: signal to evaluate
%          - x_ref: reference clean signal
%
% Output:  - snr: SNR of signal x
% ------------------
%
% Author: Lucas Rencker
% Last update: 28/03/18

snr = 10*log10((sum(x_ref.^2)+eps)/sum((x-x_ref).^2)+eps);

snr = min(snr,80); % define 80dB as "perfect signal" (to avoid inf)

return
