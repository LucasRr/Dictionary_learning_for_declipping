function [A,cost] = IHT_inpainting(Y,reliable_samples_mat,D,alg_param)
% Perform declipping using Iterative Hard Thresholding for inpainting
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
%         
        
%% Initialize parameters

if ~isfield(alg_param, 'A_init')
    alg_param.A_init = zeros(size(D,2),size(Y,2));
end

if ~isfield(alg_param, 'loud')
    alg_param.loud = 0;
end

mu = 1/norm(D)^2; % gradient descent parameter

cost = NaN(alg_param.Nit+1,1); % save cost at each iteration 

%% Declip

% initialize sparse coefficient matrix:
A = alg_param.A_init;

% compute residual:
ResidualMat = Y-D*A;
ResidualMat(~reliable_samples_mat) = 0; % residual is zero on the clipped samples

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
    
    % compute residual:
    ResidualMat = Y-D*A;
    ResidualMat(~reliable_samples_mat) = 0; % residual is zero on the clipped samples
    
    % compute cost:
    cost(it+1) = sum(sum(ResidualMat.^2));
    
    if alg_param.loud
        fprintf('it = %d, cost: %.3f\n', it, cost(it+1))
    end

end


end
