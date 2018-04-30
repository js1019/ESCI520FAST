function plotSpecIm(specT, f, specWindows, ind)
figure
imagesc(specT,f,log(squeeze(specWindows(ind,:,:))))
set(gca,'YDir','normal')
colorbar()
colormap('jet')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
caxis([-5 5])
title('Sample Spectral Image')
end
