% Calculate spectrogram of time series in the manner of Yoon et al 2015
% by Kate

clear all; clc; 
% Load time series data
load('../../data/SHcut.mat');

% Adjustable parameters for spectrogram
% Comments show parameters paper used
windowLength = 50; % Samples per window for STFT (paper = 200)
overlapLength = windowLength-2; % Sample lag in windows for STFT (paper = 2)
nfft = 4095; % Number of DFT points for STFT (paper = unknown)
dec = 100; % Decimate data to this many samples per second (paper = 20)
series = SHZc; % Component to take transform of
downSamp = 32; % Downsample to this number of frequency bins (paper = 32)

% Adjustable parameters for spectral image creation
specWindowLength = 64; % Samples per spectral image (paper = 100)
specOverlapLength = specWindowLength-5; % Spectral image sample lag (paper = 10)

% spectrogram is the full spectrogram
% specWindows is the spectrogram divided into individual spectral images
[spectrogram, specWindows] = makeSpec(windowLength, overlapLength, nfft, dec, series, downSamp, specWindowLength, specOverlapLength, tc);


