% Calculate spectrogram of time series in the manner of Yoon et al 2015

% Load time series data
load('../data/SHcut.mat');

% Adjustable parameters for spectrogram
% Comments show parameters paper used
windowLength = 50; % Samples per window for STFT (paper = 200)
overlapLength = windowLength-2; % Sample lag in windows for STFT (paper = 2)
nfft = 127; % Number of DFT points for STFT (paper = unknown)
% For downsampling, make this a multiple of 32 - 1
dec = 500; % Decimate data to this many samples per second (paper = 20)
series = SHZc; % Component to take transform of
downSamp = 32; % Downsample to this number of frequency bins (paper = 32)

% Adjustable parameters for spectral image creation
specWindowLength = 128; % Samples per spectral image (paper = 100)
% For downsampling, make this a multiple of 64
specOverlapLength = specWindowLength-5; % Spectral image sample lag (paper = 10)
specDownSamp = 64;
% spectrogram is the full spectrogram
% specWindows is the spectrogram divided into individual spectral images
[spectrogram, specWindows] = makeSpec(windowLength, overlapLength, nfft, dec, series, downSamp, specWindowLength, specOverlapLength, specDownSamp, tc);
save('../data/sampleSpectralWindows.mat', 'specWindows');

function [powerDown, specWindows] = makeSpec(windowLength, overlapLength, nfft, dec, series, downSamp, specWindowLength, specOverlapLength, specDownSamp, tc)
% Decimate data
% In the paper, they decimate the data from 100 samples/second
% to 20 samples/second
[dt, series, tNew] = decimate(dec, series, tc);


% Compute the STFT, frequencies, and time
% Use Hamming window for tapering
% Default frequency content is onesided
samplePerS = int32(1/dt);
[s,f,t] = spectrogram(series,hamming(windowLength),overlapLength,nfft,samplePerS);

% Compute spectrogram from STFT by calculating the power of 
% the STFT
power = specCalc(s);

% Downsample spectrogram to n frequency bins
% In the paper, they downsample to 32 frequency bins
n = floor(size(power,1)/downSamp);
powerDown = downsample(power,n);

% Calculate time vector for spectrogram
T = t+tc(1);

% Plot original time series and spectrogram
% Spectrogram amplitude shown on log scale
plotSpec(tNew, T, f, series, powerDown)


% Dividing the spectrogram into spectral images
% To get a single spectral image, specify first index
windowNum = floor((size(powerDown,2)-specWindowLength)/specOverlapLength);
for n = 1:windowNum
    specWindows(n,:,:) = powerDown(:,(n-1)*specOverlapLength+1:(n-1)*specOverlapLength+specWindowLength);
end

% Downsample to 64 bins on the time axis for each spectral image
k = floor(size(specWindows,3)/specDownSamp);
for n = 1:windowNum
    specWindowsDown(n,:,:) = downsample(squeeze(specWindows(n,:,:))',k)';
end
clear specWindows
specWindows = specWindowsDown;
% Time associated with each spectral image
specT = dt:dt:specWindowLength*dt;

% Plot a sample spectral image
plotSpecIm(specT, f, specWindows, 6)
end


function [dt, series, tNew] = decimate(dec, series, tc)
dt = 1/dec;
dtLength = dt/(tc(2)-tc(1));
count = 1;
for n = 1:length(series)
    if rem(n,int32(dtLength)) == 0
        temp(count) = series(n);
        tNew(count) = tc(n);
        count = count + 1;
    end 
end
clear series
series = temp;
end

function power = specCalc(s)
for n = 1:size(s,2)
    for m = 1:size(s,1)
        power(m,n) = sqrt(real(s(m,n))^2+imag(s(m,n))^2);
    end
end
end

function plotSpec(tNew, T, f, series, powerDown)
figure
subplot(2,1,1)
plot(tNew, series)
title('Time Series')
xlabel('Time (s)')
ylabel('Amplitude')
subplot(2,1,2)
imagesc(T,f,log(powerDown))
set(gca,'YDir','normal')
colorbar('southoutside')
colormap('jet')
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
caxis([-2 2])
saveas(gcf,'../figs/sampleSpectrogram.png');
end

function plotSpecIm(specT, f, specWindows, ind)
figure
imagesc(specT,f,log(squeeze(specWindows(ind,:,:))))
set(gca,'YDir','normal')
colorbar()
colormap('jet')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
caxis([-2 2])
title('Sample Spectral Image')
saveas(gcf,'../figs/sampleSpectralImage.png');
end
