%% FFT Analysis 
% Function to perform a FFT of a time-dependent signal.

function [freq, pwr] = FFTAnalysis(vector)

figure('Color', [1 1 1]);

fs = 100;  % sample frequency (Hz)
%t = 0:1/fs:30-1/fs; % 30 second span time vector


x = detrend(vector); % remove DC component
y = fft(x); % perform fast fourier transform

n = length(x);          % number of samples
P2 = (abs(y).^2)/(n*fs); % power, shift to show positive frequencies only
P1 = P2(1:n/2+1);  
P1(2:end-1) = 2*P1(2:end-1);
freq = fs*(0:(n/2))/n;
pwr = P1;
%plot results
plot(freq,P1, 'LineWidth', 1.5)
ax = gca;
ax.LineWidth = 1.75;
ax.FontName = 'CMU Serif';
ax.FontSize = 26;
box off
xlabel('Frequency (Hz)')
ylabel('Spectral Power (nm^2/Hz)')

