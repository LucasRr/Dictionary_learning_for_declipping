function Y = signal2frames(y, param)
%
% Divide signal into overlapping time-frames
% 
% Input: -y: input 1D signal
%        -param.N: frame length
%         param.hop: hop size
%         param.wa: Analysis window
% 
% Output: -Y: matrix containing each frame
% 
% ------------------
%
% Author: Lucas Rencker
% Last update: 28/03/18

N = param.N;
L = length(y);
hop = param.hop;
wa = param.wa(N);

istart = 1:hop:(L-N+1); % start index of each frame
iframes = bsxfun(@plus,repmat(istart,N,1),(0:N-1)'); % index map of each frame

Y = diag(wa) * y(iframes); % Overlapping time frames

end
