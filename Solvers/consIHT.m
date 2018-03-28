function [A,cost] = consIHT(Y,reliable_samples_mat,D,alg_param)
% Perform declipping using consistent Iterative Hard Thresholding [1]
% 
% Inputs:
%         - Y: matrix of size NxT containing T clipped signals of size N
%         - reliable_samples_mat: binary matrix describing the reliable indices of Y
%         - D: fixed dictionary
%         - alg_param.Nit: number of iterations
%         - alg_param.A_init: initial sparse matrix
%         - alg_param.loud: 1 to print the objective at each iteration, 0 otherwise
%         
% Outputs:
%         - A: sparse activation matrix
%         - cost: vector containing the value of the cost at each iteration
% 
% [1] : Consistent iterative hard thresholding for signal declipping, Kitic et al, ICASSP 2013 
%
%         
% --------------------- 
%
% Author: Lucas Rencker
%         Centre for Vision, Speech and Signal Processing (CVSSP), University of Surrey
%
% Contact: lucas.rencker@surrey.ac.uk
%         
% Reference: Consistent dictionary learning for signal declipping, L. Rencker, F. Bach, W. Wang, M. D. Plumbley,
%            Latent Variable Analysis and Signal Separation (LVA/ICA), Guildford, UK, 2018
%             
% Last update: 28/03/18
% 
% This code is distributed under the terms of the GNU Public License version 3 
% (http://www.gnu.org/licenses/gpl.txt).
%                
        
%% Initialize parameters

if ~isfield(alg_param, 'A_init')
    alg_param.A_init = zeros(size(D,2),size(Y,2));
end

if ~isfield(alg_param, 'loud')
    alg_param.loud = 0;
end

mu = 1/norm(D)^2; % gradient descent parameter

clipped_pos_mat = (~reliable_samples_mat & Y>=0);
clipped_neg_mat = (~reliable_samples_mat & Y<=0);

cost = NaN(alg_param.Nit+1,1); % save cost at each iteration 

%% Declip

% initialize sparse coefficient matrix:
A = alg_param.A_init;

% compute residual:
ResidualMat = Y-D*A;
ResidualMat(clipped_pos_mat) = max(ResidualMat(clipped_pos_mat),0);
ResidualMat(clipped_neg_mat) = min(ResidualMat(clipped_neg_mat),0);

cost(1) = sum(sum(ResidualMat.^2));

if alg_param.loud
    fprintf('initial cost: %.3f\n', cost(1))
end

it = 0;

while it < alg_param.Nit
    it = it+1;

    % gradient descent step:
    A = A + mu * D'*ResidualMat;
    
    % hard thresholding:
    A = hard_threshold(A, alg_param.K);
    
    % update residual:
    ResidualMat = Y-D*A;
    ResidualMat(clipped_pos_mat) = max(ResidualMat(clipped_pos_mat),0);
    ResidualMat(clipped_neg_mat) = min(ResidualMat(clipped_neg_mat),0);
    
    % compute cost:
    cost(it+1) = sum(sum(ResidualMat.^2));
    
    if alg_param.loud
        fprintf('it = %d, cost: %.3f\n', it, cost(it+1))
    end

end


end
