% Evaluate the declipping performance of 1 signal, using consistent
% dictionary learning

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
[y, ClippingLevel] = clip_signal(x, SNRInput);

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

reliable_samples = y<ClippingLevel & y>-ClippingLevel;
reliable_samples_mat = binary_vec2mat(reliable_samples,param);

SNRin_clipped = SNR(x(~reliable_samples),y(~reliable_samples));

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

SNRout = SNR(x,x_est_IHT);
SNRout_clipped = SNR(x(~reliable_samples),x_est_IHT(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_IHT, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('IHT for inpainting: SNR = %.3f dB',SNRout))

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
paramDL.loud = 0; % print results

[D_DL,A,cost] = DictionaryLearning_inpainting(Y,reliable_samples_mat,paramDL);

X_est_DL = D_DL*A;
x_est_DL = frames2signal(X_est_DL,param);

% figure, plot(log(cost))
% title('Objective')

SNRout = SNR(x,x_est_DL);
SNRout_clipped = SNR(x(~reliable_samples),x_est_DL(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_DL, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('Dictionary learning: SNR = %.3f dB',SNRout))

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

SNRout = SNR(x,x_est_consIHT);
SNRout_clipped = SNR(x(~reliable_samples),x_est_consIHT(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_consIHT, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('Consistent IHT: SNR = %.3f dB',SNRout))

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
paramDL.loud = 0; % print results

[D_consDL,A,cost] = consDictionaryLearning(Y,reliable_samples_mat,paramDL);

X_est_consDL = D_consDL*A;
x_est_consDL = frames2signal(X_est_consDL,param);

% figure, plot(log(cost))
% title('Objective')

SNRout = SNR(x,x_est_consDL);
SNRout_clipped = SNR(x(~reliable_samples),x_est_consDL(~reliable_samples));

fprintf('SNRout: %.3f dB\n',SNRout)
fprintf('SNR clipped improvement: %.3f dB\n',SNRout_clipped-SNRin_clipped)

figure, plot(1:L, x, 1:L, x_est_consDL, 1:L, y, '--')
legend('clean','estimate','clipped')
title(sprintf('Consistent dictionary learning: SNR = %.3f dB',SNRout))


%% Plots

samples = 17400:17500; 

figure, plot(samples, x(samples), samples, x_est_IHT(samples), samples, y(samples), '--')
title('IHT for inpainting')
figure, plot(samples, x(samples), samples, x_est_DL(samples), samples, y(samples), '--')
title('Dictionary learning for inpainting')
figure, plot(samples, x(samples), samples, x_est_consIHT(samples), samples, y(samples), '--')
title('Consistent IHT')
figure, plot(samples, x(samples), samples, x_est_consDL(samples), samples, y(samples), '--')
title('Consistent dictionary learning')


%% Evaluate how much the dictionary has "learned"

fprintf('\n    Correlation between learned dictionary and initial DCT dictionary:\n\n')

% A high correlation means the dictionary has not learned much compared to
% the initial DCT dictionary

fprintf('Using dictionary learning for inpainting: %.3f\n', sum(sum((D_DCT'*D_DL).^2))/sum(sum((D_DCT'*D_DCT).^2)))
fprintf('Using consistent dictionary learning: %.3f\n', sum(sum((D_DCT'*D_consDL).^2))/sum(sum((D_DCT'*D_DCT).^2)))

%% Listen to the results
 
% We can re-project on the unclipped samples to avoid extra distortion due
% to the sparse approximation:

x_est_IHT(reliable_samples) = y(reliable_samples);
x_est_DL(reliable_samples) = y(reliable_samples);
x_est_consIHT(reliable_samples) = y(reliable_samples);
x_est_consDL(reliable_samples) = y(reliable_samples);

% sound(x,fs)
% sound(y,fs)
% sound(x_est_IHT,fs)
% sound(x_est_DL,fs)
% sound(x_est_consIHT,fs)
% sound(x_est_consDL,fs)
