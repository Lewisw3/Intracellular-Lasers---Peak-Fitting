%% Function to low pass fourier filter the LASING_SPECTRA data 
% input arguments:
% spectra = Lasing spectra in format output from find_lasing_peaks_lw58()
% cut_off = cut off frequency (in Hz) for low pass filter

% Uses the FouFilter function develoed by Prof Thomas C. O'Haver
%https://terpconnect.umd.edu/~toh/spectrum/SignalProcessingTools.html
function L = spectraFilter(spectra, cut_off_freq)

% number of spectra
num = sum(spectra(1,:)~=0);
L = spectra;
for i = 1:num
    filtered = FouFilter(spectra(:,i)',30,0,cut_off_freq,2,0);
    L(:,i) = filtered;
end

end
