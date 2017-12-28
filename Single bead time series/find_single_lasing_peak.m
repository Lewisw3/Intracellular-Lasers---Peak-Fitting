
%fitting routine for SINGLE CELL, MULTIPLE DATA SETS

function L=find_single_lasing_peak(min_lasing_threshold, lower_lambda, upper_lambda)

%mask the timeseries data to contain information over the relevant lasing
%range
mainpath = '/Users/lewisdw11/Dropbox/Cardiomyocyte Project (Large Files)/Data/Lasing';
file = '23-Feb HL-1 (200k dish seeded on 17-Feb)/FOV1_bead1_30x_100Hz_1ms_OD13.asc';
f = fullfile(mainpath, file);
timeseries = dlmread(f);
timeseries(:,501:end) = [];
% fitting only a single peak this time
max_no_peaks = 1; 

condition = (timeseries(:,1)>lower_lambda)&(timeseries(:,1)<upper_lambda);
timeseries(~condition,:) = [];

%added some clipping of the datasets to only include 4 peaks
lambda=timeseries(:,1); % take only column 1 for the wavelengths
SPECTRA=timeseries(:,2:end); % take all other columns for the spectra data

n_data=(size(timeseries,2)-1); %number of data sets to load and process (number of columns minus 1 for the wavelength column)
disp('number of spectra loading...')
disp(n_data)
SPECTRA=detrend(SPECTRA); %detrend each spectra set 

%% Gaussian fit to each lasing peak

LASING_SPECTRA=zeros(n_data,max_no_peaks);
%LASING_STDDEV=zeros(n_data,max_no_peaks);

%PEAKS=zeros(3,n_data*max_no_peaks); %a guess at required size of array (Col 1: Peak Postion; Col 2: Variance, Col 3: Dataset number)

% j = particular spectra at a given time (i.e. a column in timeseries)
% m = those sets that have lasing resonances (i.e. strong peaks)
% i = for each peak in the spectra k = number of peaks too?
% l = number of peaks that lie outside the upper lasing limit 

h=waitbar(0,'reading data (~ several mins):');
m=0;
k=0;
l=0;
for j=1:n_data %for each dataset... (j-th dataset)

 spectrum=SPECTRA(:,j);   
    
 if (max(spectrum) < min_lasing_threshold)  %if no lasing, go to next spectra
     disp('no lasing detected in data set number:')
     disp(j)
     continue
 end
 
m=m+1; %only interested in spectra with lasing 

waitbar(j/n_data,h) %update waitbar

clear pks locs 

%[pks,locs]=findpeaks(spectrum, lambda,'MinPeakDistance',1.8);
[pks,locs]=findpeaks(spectrum, lambda,'MinPeakDistance', 1.0);

%setup windowing for Gaussian fits

%d=diff(locs); %find separarion between each peak (diff:  calculates differences between adjacent elements of X along the first array dimension whose size does not equal 1)
%a=zeros((size(d,1) + 2),1); %column vector with (number of peak locations + 1) rows 
%a(1)=d(1); % put first peak separation in first element of a
%a(2:(size(d,1))+1)=d(:); % put all elements of d in the next available elements
%a(end)=d(end); % repeat the last element
%a=a./6;     %locations halfway between each peak - WHY DO WE DIVIDE BY 6??

n=size(pks,1); 

% fit Gaussian across each window
 for i=1:n
     
     k=k+1; % not sure what this is used for anymore?
    % lower_limit=locs(i)-a(i); % create a window for each peak
    % upper_limit=locs(i)+a(i+1);
     
    % mask=(lambda>lower_limit)&(lambda<upper_limit);
     
    % X=lambda(mask);
    % Y=spectrum(mask);
     
     %f=fit(X,Y,'gauss1'); %perform fit
     f=fit(lambda,spectrum,'gauss1'); %perform fit
     C=coeffvalues(f);
     
       if C(2)> (upper_lambda *1.e9)
         l=l+1;
         failed(1,l)=j;
         failed(2,l)=i;
         continue
       end
       
     LASING_SPECTRA(m,i)=C(2);
     
% LW: this code used the k variable, but its been commented out now     
%     LASING_STDDEV(j,i)=(C(3))/(sqrt(2));
%        PEAKS(1,k)=C(2); %centroid
%        PEAKS(2,k)=(C(3))/(sqrt(2)); %standard deviation
%        PEAKS(3,k)=j; %spectrum number
     
   
 end
 
end
close(h)

clear mask
% mask=PEAKS(1,:)~=0;
% 
% Q(1,:)=PEAKS(1,mask); %centre of gaussian
% Q(2,:)=PEAKS(2,mask); %variance (SQRT for STD DEV)
% Q(3,:)=PEAKS(3,mask); %time (spectrum number)


L=LASING_SPECTRA;

end

