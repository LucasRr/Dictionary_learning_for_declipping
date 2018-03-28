
% This code compares 4 different approaches for signal declipping:
%
% - IHT_inpainting.m discards the clipped samples and performs sparse decomposition 
%     on the unclipped samples, using IHT and a fixed DCT dictionary
% - DictionaryLearning_inpainting.m discards the clipped samples and performs a 
%     gradient descent-based dictionary learning on the unclipped samples
% - consIHT.m performs consistent IHT for declipping, using a fixed DCT
% dictionary [1]
% - consDictionaryLearning.m performs consistent dictionary learning for
% signal declipping, as proposed in [2]
%
%
% References:
% [1]: Consistent iterative hard thresholding for signal declipping, 
%     S. Kitic, L. Jacques, N. Madhu, M. P. Hopwood, A. Spriet, C. De Vleeschouwer, ICASSP, 2013
% 
% [2]: Consistent dictionary learning for signal declipping, 
%     L. Rencker, F. Bach, W. Wang, M. D. Plumbley,
%     Latent Variable Analysis and Signal Separation (LVA/ICA), Guildford, UK, 2018
% 
% --------------------- 
%
% Author: Lucas Rencker
%         Centre for Vision, Speech and Signal Processing (CVSSP), University of Surrey
%
% Contact: lucas.rencker@surrey.ac.uk
%                    
% Last update: 28/03/18
% 
% This code is distributed under the terms of the GNU Public License version 3 
% (http://www.gnu.org/licenses/gpl.txt).
%                

close all
clear all
% clc

addpath(genpath('../Solvers/'));
addpath(genpath('../Utils/'));

%% Parameters

param.N = 256; % size of frame
param.hop = 0.25*param.N; % hop size
param.redundancyFactor = 2; % redundancy of dictionary
param.M = param.N * param.redundancyFactor; % number of atoms
param.wa = @wHamm; % analysis window
param.ws = param.wa; % synthesis window

M = param.M;

%% Generate DCT dictionary:

D_DCT = DCT_Dictionary(param);

%% Read signal

filename = '../glockenspiel.wav';

[x, fs] = audioread(filename);

x = x/max(abs(x)); % normalize signal

%% Clip signal:

SNRInput = 3; % desired input SNR
[y, ~] = clip_signal(x, SNRInput);

SNRin = SNR(x,y);
fprintf('Input SNR: %.3f dB\n',SNRin)
 
%% Decompose signal into overlapping time-frames:

Y = signal2frames(y,param);
Nframes = size(Y,2);

% crop signals:
L = length(frames2signal(Y,param)); % length of signal
y = y(1:L);
x = x(1:L);

%% Detect reliable samples:

% Detect clipping level:
ClippingLevel = max(abs(y));

reliable_samples = y<ClippingLevel & y>-ClippingLevel;
reliable_samples_mat = binary_vec2mat(reliable_samples,param);

SNRin_clipped = SNR(x(~reliable_samples),y(~reliable_samples));

fprintf('%.1f percent of clipped samples\n', sum(~reliable_samples)/length(x)*100)

%% Reconstruct signal using IHT for inpainting:

fprintf('\n    IHT for inpainting:\n')

alg_param.K = 32; % number of non-zero coefficients
alg_param.Nit = 50; % max number of iterations
alg_param.loud = 0; % 1 to print the results
alg_param.A_init = zeros(M,Nframes); % initialize sparse matrix

[A,cost] = IHT_inpainting(Y,reliable_samples_mat,D_DCT,alg_param);

X_est_IHT = D_DCT*A;
x_est_IHT = frames2signal(X_est_IHT,param);

% figure, plot(log(cost))
% title('Objective')

SNRout_IHT = SNR(x,x_est_IHT);
SNRout_clipped = SNR(x(~reliable_samples),x_est_IHT(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout_IHT)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_IHT, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('IHT for inpainting: SNR = %.2f dB',SNRout_IHT))
axis tight

%% Reconstruct signal using dictionary learning for inpainting:

fprintf('\n    Dictionary learning for inpainting:\n')

% DL parameters:
paramDL.K = 32; 
paramDL.Nit = 50; % number of iterations
paramDL.Nit_sparse_coding = 20; % number of iterations sparse coding step
paramDL.Nit_dict_update = 20; % number of iterations dictionary update step
paramDL.warm_start = 1; % 1 to perform warm start at each iteration
paramDL.A_init = zeros(M,Nframes); % initialize sparse coefficient matrix
paramDL.D_init = DCT_Dictionary(param); % initialize dictionary
paramDL.loud = 1; % print results

[D_DL,A,cost] = DictionaryLearning_inpainting(Y,reliable_samples_mat,paramDL);

X_est_DL = D_DL*A;
x_est_DL = frames2signal(X_est_DL,param);

% figure, plot(log(cost))
% title('Objective')

SNRout_DL = SNR(x,x_est_DL);
SNRout_clipped = SNR(x(~reliable_samples),x_est_DL(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout_DL)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_DL, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('Dictionary learning: SNR = %.2f dB',SNRout_DL))
axis tight

%% Reconstruct signal using consIHT:

fprintf('\n    Consistent IHT:\n')

alg_param.K = 32; % number of non-zero coefficients
alg_param.Nit = 50; % max number of iterations
alg_param.loud = 0; % 1 to print the results
alg_param.A_init = zeros(M,Nframes); % initialize sparse matrix

[A,cost] = consIHT(Y,reliable_samples_mat,D_DCT,alg_param);

X_est_consIHT = D_DCT*A;
x_est_consIHT = frames2signal(X_est_consIHT,param);

% figure, plot(log(cost))
% title('Objective')

SNRout_consIHT = SNR(x,x_est_consIHT);
SNRout_clipped = SNR(x(~reliable_samples),x_est_consIHT(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout_consIHT)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_consIHT, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('Consistent IHT: SNR = %.2f dB',SNRout_consIHT))
axis tight

%% Reconstruct signal using consDL:

fprintf('\n    Consistent dictionary learning:\n')

% DL parameters:
paramDL.K = 32; 
paramDL.Nit = 50; % number of iterations
paramDL.Nit_sparse_coding = 20; % number of iterations sparse coding step
paramDL.Nit_dict_update = 20; % number of iterations dictionary update step
paramDL.warm_start = 1; % 1 to perform warm start at each iteration
paramDL.A_init = zeros(M,Nframes); % initialize sparse coefficient matrix
paramDL.D_init = DCT_Dictionary(param); % initialize dictionary
paramDL.loud = 1; % print results

[D_consDL,A,cost] = consDictionaryLearning(Y,reliable_samples_mat,paramDL);

X_est_consDL = D_consDL*A;
x_est_consDL = frames2signal(X_est_consDL,param);

% figure, plot(log(cost))
% title('Objective')

SNRout_consDL = SNR(x,x_est_consDL);
SNRout_clipped = SNR(x(~reliable_samples),x_est_consDL(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout_consDL)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_consDL, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('Consistent dictionary learning: SNR = %.2f dB',SNRout_consDL))
axis tight


%% Plots

samples = 46800:46900; 

% percentage of missing samples:
% sum(~reliable_samples(samples))/length(samples)*100;

figure, plot(samples, x(samples), samples, x_est_IHT(samples), samples, y(samples), '--'), axis tight
title(sprintf('IHT for inpainting: SNR = %.2f dB',SNRout_IHT))
legend('clean','estimate','clipped')
figure, plot(samples, x(samples), samples, x_est_DL(samples), samples, y(samples), '--'), axis tight
title(sprintf('Dictionary learning: SNR = %.2f dB',SNRout_DL))
legend('clean','estimate','clipped')
figure, plot(samples, x(samples), samples, x_est_consIHT(samples), samples, y(samples), '--'), axis tight
title(sprintf('Consistent IHT: SNR = %.2f dB',SNRout_consIHT))
legend('clean','estimate','clipped')
figure, plot(samples, x(samples), samples, x_est_consDL(samples), samples, y(samples), '--'), axis tight
title(sprintf('Consistent dictionary learning: SNR = %.2f dB',SNRout_consDL))
legend('clean','estimate','clipped')

%% Listen to the results
 
% We can re-project on the unclipped samples to avoid extra distortion due
% to the sparse approximation:

x_est_IHT(reliable_samples) = y(reliable_samples);
x_est_DL(reliable_samples) = y(reliable_samples);
x_est_consIHT(reliable_samples) = y(reliable_samples);
x_est_consDL(reliable_samples) = y(reliable_samples);

% Uncomment to listen to the result:

% sound(x,fs)
% sound(y,fs)
% sound(x_est_IHT,fs)
% sound(x_est_DL,fs)
% sound(x_est_consIHT,fs)
% sound(x_est_consDL,fs)
