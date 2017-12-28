%% Function to find peak positions from a series of spectra

% The function takes in a spectral timeseries in ascii format and then
% finds/fits the peaks in each one. 

% Output:
% L = List of peak positions to be used for refractive index sensing
% int = List of peak intensity
% err = List of error intervals

% Adapted from code developed by Isla Barnard 
% Lewis Woolfson - June 2017 - lewis.woolfson@gmail.com

function [L,int,err]=find_lasing_peaks_lw58(min_peak_prominence, lower_lambda, upper_lambda)

%% User Setup 

mainpath = '/Users/?';
file = 'FOV3_B2_20x_100Hz_1ms_OD16_6hr.asc';
f = fullfile(mainpath, file);

% For spectral data over 10 minutes long you need to specify the range of
% data to read from the delimited file. For data shorter than 
% 10 minutes the command 'dlmread(f)' is sufficient.
% Comment and uncomment as required.

% FOR LARGE TIMESERIES (>10mins)
%startRow = 0; % start row number (first wavelength bin, FIXED)
%startCol = 0; % start column number (first spectrum
%maxRow = 1599; % end row number (last wavelength bin, FIXED)
%maxCol = 60000; % end column number (last spectrum)
%timeseries = dlmread(f, '', [startRow startCol maxRow maxCol]);

% FOR SMALL TIMESERIES (<10mins)
timeseries = dlmread(f);

%trim the data set
timeseries(:,100:end) = [];

% Error confidence interval (default 95%)
confidence = 0.95;

% maximum number of peaks
max_no_peaks = 18;

% Generate an output file of the peak positions in .txt format?
genOutputFile = false;

%% Gaussian fit to each lasing peak

% Remove wavelengths outside region of interest
condition = (timeseries(:,1)>lower_lambda)&(timeseries(:,1)<upper_lambda);
timeseries(~condition,:) = [];

lambda=timeseries(:,1); % take only column 1 for the wavelengths
SPECTRA=timeseries(:,2:end); % take all other columns for the spectra data
n_data=(size(timeseries,2)-1); %number of data sets 
fprintf('-----------------------------------------------\n')
fprintf('Gaussian Peak Fitting\n')
fprintf('Total number of spectra to process = %d\n', n_data)
fprintf('Loading...')
SPECTRA=detrend(SPECTRA); %detrend each spectra set 

LASING_SPECTRA=zeros(n_data,max_no_peaks); % peak positions output
LASING_INTENSITY=zeros(n_data,max_no_peaks); % peak intensity output
LASING_ERROR=zeros(n_data,max_no_peaks); % confidence intervals output

str = sprintf('Fitting %d spectra...\n', n_data); 
h = waitbar(0, str);
m=0;
k=0;
l=0;

for j=1:n_data %for each dataset... (j-th dataset)

    spectrum=SPECTRA(:,j); 
 
    % calculate the minimum lasing threshold as 5 standard deviations from
    % the mean (arb.)
    min_lasing_threshold = mean(spectrum)+5*std(spectrum);   

    %if no lasing, go to next spectra
    if (max(spectrum) < min_lasing_threshold) 
     disp('no lasing detected in data set number:')
     disp(j)
     continue
    end
 
    m=m+1; 
    
    waitbar(j/n_data,h) 

    clear pks locs 

    % locate peaks in spectrum
    [pks,locs,w,p]=findpeaks(spectrum, lambda,'MinPeakDistance',1.0,...
        'MinPeakProminence', min_peak_prominence);

    % number of peaks 
    n=size(pks,1); 

    if (n > 1) 
        % create a window around each peak
        d=diff(locs); 
        a=zeros((size(d,1) + 2),1); 
        a(1)=d(1); 
        a(2:(size(d,1))+1)=d(:); 
        a(end)=d(end); 
        a=a./6;   
    end 
  
    % fit Gaussian across each window
    for i=1:n
         k=k+1; 
         
        if (n > 1)  
            lower_limit=locs(i)-a(i); 
            upper_limit=locs(i)+a(i+1);
            mask=(lambda>lower_limit)&(lambda<upper_limit);         
            X=lambda(mask);
            Y=spectrum(mask);
        else 
            X = lambda; 
            Y = spectrum;
        end    
            
        try 
            fo = fitoptions('gauss1', 'Lower', [0, lower_lambda, 0], ...
                'Upper', [Inf, upper_lambda, Inf]);
             f=fit(X,Y,'gauss1',fo); %perform 1-term gaussian fit
             C = coeffvalues(f);  % extract coefficients
             Ci = confint(f,confidence); %extract confidence intervals  
             
             % redundant code, fit options won't allow it. 
             if (C(2) < lower_limit || C(2) > upper_limit)
                 % if it returns a centroid value outside the window
                 error('Gauss1 returned a bad fit for peak number %d at spectra number %d. Attempting gauss2 fit...', i,m);
             end                    
             
        catch            
            warning('Gauss1 returned a bad fit for peak number %d at spectra number: %d. Attempting gauss2 fit...', i,m);
            % if the 1-term gaussian fit didn't work then try a 2-term one
            try                         
                f=fit(X,Y,'gauss2');
                C = coeffvalues(f);
                Ci = confint(f,confidence);
                if (C(2) < lower_limit || C(2) > upper_limit)
                    error('Gauss2 failed for peak number %d at spectra number %d.', i,m);
                end     
            catch
                warning('Two term gaussian fit failed.');
                C(2) = 0;
                C(1) = 0;
                Ci(2,2) = 0; 
                Ci(1,2) = 0;
            end            
        end
      
       % if the fit result is greater than maximum wavelength (redundant)
       if C(2)> (upper_lambda *1.e9)
         l=l+1;
         failed(1,l)=j;
         failed(2,l)=i;
         continue
       end
     
       % add fit results from this spectrum to output variables
       LASING_ERROR(m,i)=(Ci(2,2)-Ci(1,2))/2;
       LASING_INTENSITY(m,i)=C(1);  
       LASING_SPECTRA(m,i)=C(2);   
    end
end

close(h)
clear mask
err = LASING_ERROR;
int = LASING_INTENSITY;
L=LASING_SPECTRA;

if (genOutputFile)
    fileName = strcat(file(1:end-4), '_peaks.txt');
    save(fileName, 'L', '-ascii');
end
fprintf(' Done!\n')
fprintf('-----------------------------------------------\n')

end

