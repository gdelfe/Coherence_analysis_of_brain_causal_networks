
clear all; close all; clc;

addpath('/mnt/pesaranlab/People/Gino/code');


Lfp = importdata('LFP_all_Shaoyu.txt');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% SPECTRAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Spectrum and PDS %%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%
% ONE CHANNEL ONLY %%%%%
% %%%%%%%%%%%%%%%%%%%%%%

% Compute the Lfp for one channel, all trials
% Lfp =  trialLfp(Trials,'FEF',1,1,'TargsOn',[-500,1300]); % Extract the Lfp fields
timeVals=linspace(-500,1300,size(Lfp,3));


W = 5;
% Compute the spectrum for each trial. Format: iTrial x times
[spec, f, err] = dmtspec_GINO(sq(Lfp(:,:,:)),[1800/1e3,W],1e3,200);

% Plot the spectrum for each trial
figure(1); plot(f,log(spec(:,:)'))

% Plot the mean spectrum across trials
figure(2); plot(f,log(mean(spec)'))

% Compute the Spectrogram - f vs time 
% TFSPEC(X, TAPERS, SAMPLING, DN, FK, PAD, PVAL, FLAG, CONTFLAG, ERRORBAR)

[spectrum, f2] = tfspec_GINO(Lfp,[0.5 5],1e3,0.005,60,2,0.05,1); % one single call
% [spectrum, f2] = tfspec_GINO(Lfp,[0.5,2],1e3,0.05,200,2,0.05,0); % one single call: FLAG = 0, computes trial by trial

figure(3); tvimage(log(spectrum)); colorbar;
% figure(3); tvimage(log(spectrum(1,:,:))); colorbar; % for trial 1
dlmwrite('spectrum.txt',log(spectrum),'delimiter','\t')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single projection on tapers %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take the single projection on the tapers (Trial by Trial)
% SpecK = 4D complex vector: Trial, Times, K(taper), frequency 
[SpecK, f] = tfsp_proj_GINO(Lfp,[0.5,6],1e3,0.05,200,2);
Sx = sq(SpecK);
Sx_Length = Sx.*conj(Sx);
figure(4); tvimage(log(Sx_Length(1,:,:))); colorbar; % for trial 1
Sx_Length_mean = mean(Sx_Length,1);
figure(5); tvimage(log(Sx_Length_mean)); colorbar; % spectrogram across trials




% Take the single projection on the tapers (Trial by Trial)
% SpecK = 4D complex vector: Trial, Times, K(taper), frequency 
[SpecK, f] = tfsp_proj_GINO(Lfp,[0.5,6],1e3,0.05,200,2);
Sx_Length = SpecK.*conj(SpecK);
figure(4); tvimage(log(Sx_Length(1,:,1,:))); colorbar; % for trial 1
Sx_Length_mean = sq(mean(Sx_Length,1));
figure(6); tvimage(log(Sx_Length_mean(:,5,:))); colorbar; % spectrogram across trials







