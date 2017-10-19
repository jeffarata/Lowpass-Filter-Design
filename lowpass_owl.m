% Jeff Arata 
% 10/17/17

% This project and the associated files were provided by Joe Hoffbeck and
% are found in his paper "Enhance Your DSP Course With These Interesting
% Projects.pdf"

% In this project, we look at the recording of an own hooting that has some
% truck's back up beeper in the background. We will design/apply lowpass 
% filter in order to reduce that beeper sound.

clear;
clc;

[x, fs] = audioread('owl_beep.wav');
soundsc(x, fs);                         % Listen to the file

figure(1)
spectrogram(x, 512, [], 512, fs, 'yaxis');  % Spectrogram

[PSD, F] = pwelch(x, [], [], [], fs);
logPSD = 10*log10(PSD);                     % Get PSD for plotting

figure(2)           % Plotting is just nice to help see the fundamental 
plot(F, logPSD)     % frequencies and pick a tolerance to find those values.

tol_db = -25;

fun_freq = distinguish_notes( x, fs, tol_db );  % Gets frequencies of the distinct fundamentals

theta = fun_freq(2) / (fs/2);   % Angle in radians of the freq. we don't want
H_root = exp(i * theta);        % Root as a complex number

H_roots_prev = i*H_root;        % Initialization to make sure it's not the same as H_root


% This loop generates a filter of successively larger orders, get's the
% zeros of the transfer function and then compares them to the root we want
% in our transfer function - H_root - or as close to it as possible.
% Comparisons are of the angle, real parts, and imaginary parts of their
% roots. A new filter with a different order is stored only if the sum of
% the differences between the true root and the closest in the filter is
% smaller than that of the previous filter.

% The cutoff frequency of 730 was taken from Joe Hoffbeck's discussion of
% the project in his pdf, but I decided to see if I could find a more
% accurate filter by getting a specific order.

for ii = 1:300
    b = fir1(ii, 730/(fs/2));   
    H_roots = roots(b);
    
    dtheta_new = min(abs(theta - angle(H_roots)));
    dreal_new = min(abs(real(H_root) - real(H_roots)));
    dimag_new = min(abs(imag(H_root) - imag(H_roots)));
    sum_new = dtheta_new + dreal_new + dimag_new;
    
    dtheta_old = min(abs(theta - angle(H_roots_prev)));
    dreal_old = min(abs(real(H_root) - real(H_roots_prev)));
    dimag_old = min(abs(imag(H_root) - imag(H_roots_prev)));
    sum_old = dtheta_old + dreal_old + dimag_old;
    
    if sum_new < sum_old
        order = ii;
        H_roots_prev = H_roots;
    end
end
  
b = fir1(order, 730/(fs/2));    % Best filter found

figure(3)                       % Plot zeros and poles
zplane(b, 1)

figure(4)                       % Magnitude and Phase of Frequency Response
freqz(b, 1, 512)    

y = filter( b, 1, x );          % Filter the sound

soundsc(y, fs);                 % Play the sound
