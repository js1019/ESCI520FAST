function [powerDown, specWindows] = makeSpec(windowLength, overlapLength, nfft, dec, series, downSamp, specWindowLength, specOverlapLength, tc)
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

% Time associated with each spectral image
specT = dt:dt:specWindowLength*dt;

% Plot a sample spectral image
plotSpecIm(specT, f, specWindows, 3)
end
