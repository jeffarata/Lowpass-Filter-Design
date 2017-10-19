function [ output ] = distinguish_notes( x, fs, tol_db )
% This function takes an input signal x and its sampling frequency fs and
% returns distinct/separate fundamental frequencies of the recording.
% Tol_db is a decibel level tolerance in order to use the frequencies with
% the highest values in the log scale of the PSD.

[PSD, F] = pwelch(x, [], [], [], fs);
logPSD = 10*log10(PSD);                     % log scale of PSD

filt_freq = F(logPSD(F < 4000) > tol_db)';  % Filter to frequencies greater than tol_db

piano_freq = 440 * 2 .^ (([1:88] - 49)/12); % Piano frequencies for comparison

% Half the frequency distance between notes for more accurate tolerances.
adap_freq_tol = 400 * (2.^(([2:89] - 49)/12) - 2.^(([1:88] - 49)/12)) / 2;

% This loop sets any octaves of the lower frequencies to zero. Once it
% reaches a zero value, it assumes the rest of the entries are harmonics
% and so cuts them out. 

for ii = 1:length(filt_freq)
    [val, key] = min(abs(filt_freq(ii) - piano_freq));
    
    filt_freq( abs(filt_freq(ii)*2*ones(1,length(filt_freq)) - filt_freq) < adap_freq_tol(key) ) = 0;
       
    if filt_freq(ii) == 0
       filt_freq = filt_freq(1:ii-1);
       break
    end    
end

% PSD values to compare at the interested frequencies.
for kk = 1:length(filt_freq)
   PSDvals(kk) = logPSD(F == filt_freq(kk)); 
end

ii=1;   % initialization

while ii < length(filt_freq)    % This block boils it down to the frequencies 
                                % where the highest PSD values are
    [tmp, key] = min(abs(filt_freq(ii) - piano_freq));
    
    % if the difference in freq between two consecutive entries is close
    % enough, the one with the smaller PSD db value is removed.

    if abs(filt_freq(ii+1)-filt_freq(ii)) < adap_freq_tol(key)
        [val, idx] = min(PSDvals(ii:ii+1));
        idx = idx + ii - 1;
        filt_freq(idx) = 0;
        filt_freq = filt_freq( filt_freq ~= 0 );
        PSDvals(idx) = -120;
        PSDvals = PSDvals( PSDvals ~= -120 );
    else
        ii=ii+1;    % Increment
    end
end

output = filt_freq;

end

