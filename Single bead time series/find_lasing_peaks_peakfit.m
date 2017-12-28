%% Function to find peak positions from a series of spectra

% Adapted from code developed by Isla Barnard
% Lewis Woolfson - Jun 2017 - lewis.woolfson@gmail.com

function L = find_lasing_peaks_peakfit(min_peak_prominence, lower_lambda, upper_lambda)
% The function takes in a spectral timeseries in ascii format and then
% finds/fits the peaks in each one using the peakfit() function

% Output:
% L = List of peak positions to be used for refractive index sensing
% int = List of peak intensity
% err = List of error intervals

%% User Setup 

mainpath = '/Users/?';
file = 'FOV3_B2_20x_100Hz_1ms_OD16_6hr.asc';
f = fullfile(mainpath, file);
timeseries = dlmread(f);
timeseries(:,10:end) = [];

% Remove wavelengths outside region of interest
condition = (timeseries(:,1)>lower_lambda)&(timeseries(:,1)<upper_lambda);
timeseries(~condition,:) = [];

% Error confidence interval (default 95%)
confidence = 0.95;

% Minimum separation of the peaks
minPeakDistance = 1.0;

% Maximum number of peaks
max_no_peaks = 18;

% Peak Splitting Options
NumPeaks_TE = 1; % Number of peaks in an individual TE mode (default: 1)
NumPeaks_TM = 2; % Number of peaks in an individual TE mode (default: 1)

% autozero=0 (default) does not subtract baseline from data segment;
% autozero=1 interpolates a linear baseline from the edges of the data segment and subtracts it from the signal; 
% autozero=2, like mode 1 except that it computes a quadratic curved baseline; 
% autozero=3 compensates for a flat baseline without reference to the signal itself 
autozeroTE = 0;
autozeroTM = 0;

% Generate an output ascii file? 
genOutputFile = false;

%% Gaussian fit to each lasing peak

lambda=timeseries(:,1); % take only column 1 for the wavelengths
SPECTRA=timeseries(:,2:end); % take all other columns for the spectra data
n_data=(size(timeseries,2)-1); %number of data sets 
fprintf('-----------------------------------------------\n')
fprintf('Gaussian Peak Fitting\n')
fprintf('Total number of spectra to process = %d\n', n_data)
fprintf('Loading...')
SPECTRA=detrend(SPECTRA); %detrend each spectra set 

LASING_SPECTRA=zeros(n_data,max_no_peaks); % peak positions output
%LASING_INTENSITY=zeros(n_data,max_no_peaks); % peak intensity output
%LASING_ERROR=zeros(n_data,max_no_peaks); % confidence intervals output

% setup waitbar
str = sprintf('Fitting %d spectra...\n', n_data); 
h = waitbar(0, str);

m=0;
k=0;
l=0;

for j=1:n_data %for each spectrum...

    spectrum=SPECTRA(:,j); 

    % calculate the minimum lasing threshold as 5 standard deviations from
    % the mean (arb.)
    min_lasing_threshold = mean(spectrum)+1*std(spectrum);   
    
    if (max(spectrum) < min_lasing_threshold) 
       fprintf('no lasing detected in spectrum number: %d\n', j)
       continue
    end
    
    m =m+1;
    
    waitbar(j/n_data,h)

    clear pks locs 
    
    % locate peaks in spectrum
    [pks,locs,w,p]=findpeaks(spectrum, lambda,'MinPeakDistance',...
        minPeakDistance,'MinPeakProminence', min_peak_prominence);
      
    % number of peaks 
    n=size(pks,1);   
   
    if (n > 2) 
         % separate into TE and TM modes
        [te,tm] = te_or_tm(locs);
    elseif (n == 2)
         if(p(1) > p(2)) 
             te = locs(1);
             tm = locs(2);
         else
             te = locs(2);
             tm = locs(1);
         end 
    else
        te = locs; % adjust to TE or TM as necessary when fitting a single peak
    end
    
    % for spectra containing more than one peak
    if (n > 1)     
        % for making a window around each peak
        d=diff(locs); 
        a=zeros((size(d,1) + 2),1); 
        a(1)=d(1); 
        a(2:(size(d,1))+1)=d(:); 
        a(end)=d(end); 
        a=a./6;   
    end         
    
    for i = 1:n % for each peak in the spectrum...     
        
        k = k+1;
         
        if (n > 1) % create window around each peak
            lower_limit=locs(i)-a(i); 
            upper_limit=locs(i)+a(i+1);
            mask=(lambda>lower_limit)&(lambda<upper_limit);         
            X=lambda(mask);
            Y=spectrum(mask);
        else  % use the entire spectrum if it only contains one peak
            X = lambda; 
            Y = spectrum; 
        end 
        
        % if it is a TE mode then do a fit
        % otherwise if it is a TM mode then do a different fit
        % use peakfit() for both modes
        % use the find function to get the first non zero entry of an array
           
        
        try 
            if(ismember(locs(i),te)) 
                % TE Mode
               FitResults = peakfit([X Y],0,0,NumPeaks_TE,0,0,0,0,...
                 autozeroTE,0,0);
            else
                % TM Mode
               FitResults = peakfit([X Y],0,0,NumPeaks_TM,0,0,5,0,...
                 autozeroTM,0,0);
            end 
            
            rows2delete = [];
            for p = 1:size(FitResults,1)
                centroid = FitResults(p,2);
                if (centroid < lower_limit || centroid > upper_limit)
                    %error('Bad peakfit() result...')
                    rows2delete(end+1) = p;
                end
            end                  
            if (isempty(rows2delete)==0)
                FitResults(rows2delete,:) = []; % delete rows
            end            
        catch 
            warning('peakfit() failed for peak number %d at spectra number: %d. Attempting gauss1 fit...', i,m);
              try 
                fo = fitoptions('gauss1', 'Lower', [0, lower_lambda, 0],...
                'Upper', [Inf, upper_lambda, Inf]);
                f=fit(X,Y,'gauss1',fo); %perform 1-term gaussian fit
                FitResults = coeffvalues(f);
              catch     
                warning('Gauss1 fit failed for TM peak number %d at spectra number: %d. Attempting gauss2 fit...', i,m);
                try 
                  f=fit(X,Y,'gauss2');
                  FitResults = coeffvalues(f);
                catch 
                  warning('Gauss2 fit failed for TM peak number %d at spectra number: %d. Attempting gauss2 fit...', i,m);
                  FitResults = [0 0 0 0];
                end
              end   
        end
        
        for p = 1:size(FitResults,1)
            % extract peak centroid position
            centroid = FitResults(p,2);
            next = find(~LASING_SPECTRA(m,:) ,1);
            LASING_SPECTRA(m,next) = centroid;
        end 
    end
end

close(h)
clear mask
L = LASING_SPECTRA;
if (genOutputFile)
    fileName = strcat(file(1:end-4), '_peaks.txt');
    save(fileName, 'L', '-ascii');
end
fprintf(' Done!\n')
fprintf('-----------------------------------------------\n')
end







