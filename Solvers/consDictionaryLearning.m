function [D,A,cost] = consDictionaryLearning(Y,reliable_samples_mat,paramDL)
% Perform declipping using consistent dictionary learning
% 
% Inputs:
%         - Y: matrix of size NxT containing T clipped signals of size N
%         - reliable_samples_mat: binary matrix describing the reliable indices of Y
%         - paramDL.K: number of non-zero atoms 
%         - paramDL.Nit: number of dictionary learning iterations
%         - paramDL.Nit_sparse_coding: number of iterations sparse coding step
%         - paramDL.Nit_dict_update: number of iterations dictionary update step
%         - paramDL.warm_start: 1 to perform warm start at each iteration
%         - paramDL.A_init: initial sparse coefficient matrix
%         - paramDL.D_init: initial dictionary
%         - paramDL.loud: 1 to print results
%         
% Outputs:
%         - D: estimated dictionary
%         - A: sparse activation matrix
%         - cost: vector containing the value of the cost at each iteration
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

%% Initialize parameters:

if ~isfield(paramDL, 'warm_start')
    paramDL.warm_start = 1;
end

if ~isfield(paramDL, 'loud')
    paramDL.loud = 0;
end

A = paramDL.A_init;
D = paramDL.D_init;

cost = NaN(2*paramDL.Nit+1,1); % save cost at each iteration 

clipped_pos_mat = (~reliable_samples_mat & Y>=0);
clipped_neg_mat = (~reliable_samples_mat & Y<=0);

% compute residual:
ResidualMat = Y-D*A;
ResidualMat(clipped_pos_mat) = max(ResidualMat(clipped_pos_mat),0);
ResidualMat(clipped_neg_mat) = min(ResidualMat(clipped_neg_mat),0);

cost(1) = sum(sum(ResidualMat.^2));

if paramDL.loud
    fprintf('initial cost: %.3f\n', cost(1))
end

%% DL iterations

it = 0;

while it < paramDL.Nit
    it = it+1;

    %% Sparse coding
    
    % parameters for sparse coding step:
    paramSC.K = paramDL.K;
    paramSC.Nit = paramDL.Nit_sparse_coding;
    
    if paramDL.warm_start
        paramSC.A_init = A; % warm_start
    else
        paramSC.A_init = zeros(size(A));
    end
    
    [A,cost_SC] = consIHT(Y,reliable_samples_mat,D,paramSC);
    
    cost(2*it) = cost_SC(end);
    
    if paramDL.loud
        fprintf('it %d, sparse coding step: cost: %.3f\n',it, cost(2*it))
    end
    
    %% Prune unused atoms:
    
    unused = (sum(A.^2,2) == 0);
    
    if sum(unused)>0 && paramDL.loud
        fprintf('  %d atoms pruned\n', sum(unused))
    end
    
    D(:,unused) = [];
    A(unused,:) = [];
    
    %% Dictionary Update
    
    paramDictUpdate.D_init = D; % initialize with previous estimate
    paramDictUpdate.Nit = paramDL.Nit_dict_update;
    
    [D, final_cost] = DictUpdate(Y,A,reliable_samples_mat,paramDictUpdate);
    
    cost(2*it+1) = final_cost;
    
    if paramDL.loud
        fprintf('it %d,   dict update step: cost: %.3f\n',it, cost(2*it+1))
    end
    
end


end % end DL

function [D, final_cost] = DictUpdate(Y,A,reliable_samples_mat,paramDictUpdate)


%% initialization:

clipped_pos_mat = (~reliable_samples_mat & Y>=0);
clipped_neg_mat = (~reliable_samples_mat & Y<=0);

D = paramDictUpdate.D_init;

% residual:
ResidualMat = Y-D*A;
ResidualMat(clipped_pos_mat) = max(ResidualMat(clipped_pos_mat),0);
ResidualMat(clipped_neg_mat) = min(ResidualMat(clipped_neg_mat),0);

mu = 1/norm(A)^2; % gradient descent parameter

%% Gradient descent:

for iter = 1:paramDictUpdate.Nit
    
    % gradient descent:
    D = D + mu * ResidualMat*A';
    
    % normalize atoms:
    D = bsxfun(@times,D,1./max(sqrt(sum(D.^2,1)),1)); % normalize atoms

    % update residual:
    ResidualMat = Y-D*A;
    ResidualMat(clipped_pos_mat) = max(ResidualMat(clipped_pos_mat),0);
    ResidualMat(clipped_neg_mat) = min(ResidualMat(clipped_neg_mat),0);
    
end

final_cost = sum(sum(ResidualMat.^2));

end

