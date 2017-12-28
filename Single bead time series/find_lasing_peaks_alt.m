%% Alternative peak fitting routine using the findpeaksb function
% Find and Gaussian fit spectral peaks in time series using the findpeaksb
% function developed by Prof Thomas C. O'Haver, Uni of Maryland

% offers improved functionality than the matlab fit() function

function [L,int]=find_lasing_peaks_alt(SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype,window,PeakShape,extra,NumTrials,AUTOZERO,lower_lambda,upper_lambda)

% Read in ascii file containing spectrla data
mainpath = '/Users/lewisdw11/Dropbox/Cardiomyocyte Project (Large Files)/Data/Lasing/21-Mar Primary CDM';
file = 'FOV2_B2_20x_100Hz_1ms_OD13.asc';
f = fullfile(mainpath, file);
timeseries = dlmread(f);
max_no_peaks = 18;

lambda=timeseries(:,1); % take only column 1 for the wavelengths
SPECTRA=timeseries(:,2:end); % take all other columns for the spectra data

n_data=(size(timeseries,2)-1); 
disp('number of spectra loading...')
disp(n_data)
SPECTRA=detrend(SPECTRA); %detrend each spectra set 

%% Gaussian fit to each lasing peak

LASING_SPECTRA=zeros(n_data,max_no_peaks);
LASING_INTENSITY=zeros(n_data,max_no_peaks);

h=waitbar(0,'reading data (~ several mins):');
m=0;

for j=1:n_data %for each dataset... (j-th dataset)

    spectrum=SPECTRA(:,j); 

    % calculate the minimum lasing threshold as 5 standard deviations from
    % the mean

    min_lasing_threshold = mean(spectrum)+5*std(spectrum);   
    %if no lasing, go to next spectra
    if (max(spectrum) < min_lasing_threshold)  
        disp('no lasing detected in data set number:')
        disp(j)
        continue
    end

    m=m+1; %only interested in spectra with lasing 
    waitbar(j/n_data,h) %update waitbar

 
    P = findpeaksb(lambda,spectrum,SlopeThreshold,AmpThreshold,...
        smoothwidth,peakgroup,smoothtype,window,PeakShape,extra,...
        NumTrials,AUTOZERO);
    
    locs = P(:,2)'; % peak location
    heights = P(:,3)'; % peak height
    num_peaks = size(locs,2); % number of peaks
    
    % add results to output variable
    LASING_SPECTRA(m,1:num_peaks) = locs; 
    LASING_INTENSITY(m,1:num_peaks) = heights;
end

close(h)

int = LASING_INTENSITY;
L=LASING_SPECTRA;

end

