function A_ht = hard_threshold(A,K)
% 
% Hard thresholding function. Keeps the K largest coefficients in each column of matrix A,
% and sets the rest to zero
% 
% Input:  - A: input matrix
%         - K: sparsity parameter (number of non-zero coefficient)
%         
% Output: - A_ht: sparse output matrix, with K non-zero coefficient per column
% ------------------
%
% Author: Lucas Rencker
% Last update: 28/03/18


if K == size(A,1)
    A_ht = A;
else
    A_sort = sort(abs(A),'descend');
    A_ht = A.*bsxfun(@ge,abs(A),A_sort(K,:)); % fast implementation
end
