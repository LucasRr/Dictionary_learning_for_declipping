function snr = SNR(y_ref,y)
% Signal-to-noise Ratio
%
% Inputs:
%          - y: signal to evaluate
%          - y_ref: reference clean signal
%
% Output:  - snr: SNR of signal y
% -------------------

snr = 10*log10((sum(y_ref.^2)+eps)/sum((y-y_ref).^2)+eps);

snr = min(snr,80); % define 80dB as "perfect signal" (to avoid inf)

return
