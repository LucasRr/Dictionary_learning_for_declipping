function binary_mask_mat = binary_vec2mat(binary_mask_vec, param)
% Divide binary vector mask into matrix mask, in the same way as
% signal2frames
% 
% Input: -binary_mask_vec: binary vector mask
%        -param.N: frame length
%         param.hop: hop size
% 
% Output: -binary_mask_mat: binary matrix mask
%
% ------------------
%
% Author: Lucas Rencker
% Last update: 28/03/18

N = param.N;
L = length(binary_mask_vec);
hop = param.hop;

istart = 1:hop:(L-N+1); % start index of each frame
iframes = bsxfun(@plus,repmat(istart,N,1),(0:N-1)'); % index map of each frame

binary_mask_mat = logical(binary_mask_vec(iframes)); 

end
