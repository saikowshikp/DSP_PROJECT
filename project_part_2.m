Fs = 8000; % Sampling frequency
x=load('speech_noisy.txt')'; % Input signal
y=remove_noiseband(x,Fs);

audiowrite('Part_ii_original.wav',x,Fs);
audiowrite('Part_ii_fixed.wav',y,Fs);

function x_denoised=remove_noiseband(x,Fs)
% Parameters for Welch's method
windowLength = 700; % Length of each segment/window
overlap = 0.5; % Overlap percentage (50%)

% Compute PSD
[Pxx, f] = average_PSD(x, windowLength, overlap, Fs);
dB_Pxx=10*log10(Pxx);
norm_dB_Pxx=dB_Pxx./abs(mean(dB_Pxx));
% dPdf=zeros(length(f)-1);
% for i=2:length(f)
%     dPdf(i-1)=(norm_dB_Pxx(i)-norm_dB_Pxx(i-1))/(f(i)-f(i-1));
% end

x_denoised=x;
idx=1;
Idx=1;
check=0;
while 1
    start_freq=-1;
    for i=idx:length(f)-1
        if dB_Pxx(i)>mean(dB_Pxx)
        %band may start
        if start_freq==-1 && (dB_Pxx(i)-dB_Pxx(i+1))/dB_Pxx(i) < 0.05
            start_freq=f(i);
        end
        %band is not long enough
        if start_freq ~=-1 && (dB_Pxx(i)-dB_Pxx(i+1))/dB_Pxx(i) > 0.05
            start_freq = -1;
        end
        %found the band
        if start_freq~=-1 && f(i)>=200+start_freq
            Idx=Idx+1;
            idx=i;
            h=band_stop(1001,start_freq-100,start_freq+200+100,Fs);
            x_denoised=conv(h,x_denoised);
            break;
        end
        
        end
        if i==length(f)-1
            start_freq=-1;
            check=1;
            break;
        end
    end
    if check==1
        break;
    end
end
% h=band_stop(1001,start_freq(1)-100,start_freq(1)+200+100,Fs);
% x_denoised=conv(h,x);
end

% Function to compute PSD estimate using Welch's method
function [Pxx, f] = average_PSD(x, windowLength, overlap, Fs)
    hopSize = round(windowLength * (1 - overlap));
    numSegments = floor((length(x) - windowLength) / hopSize) + 1;
    Pxx = zeros(windowLength/2 + 1, 1);
    for i = 1:numSegments
        startIdx = (i-1)*hopSize + 1;
        endIdx = startIdx + windowLength - 1;
        segment = x(startIdx:endIdx);
        segment = segment .* hamming(windowLength); %Applying Hamming window
        segment_fft = fft(segment);
        Pxx = Pxx + abs(segment_fft(1:windowLength/2 + 1)).^2;
    end
    Pxx = Pxx / numSegments;
    f = (0:(windowLength/2)) * Fs / windowLength;
end

function h=band_stop(N,f1,f2,Fs)
n = -(N-1)/2 : (N-1)/2;
h = zeros(size(n));
h(n == 0) = 1 - (f2-f1)/Fs;
h(n ~= 0) = (sin(2*pi*f1*n(n ~= 0)/Fs)./(pi*n(n ~= 0)) - sin(2*pi*f2*n(n ~= 0)/Fs)./(pi*n(n ~= 0)));
w = hamming(N)'; % Hamming window
h = h .* w; % Apply window to impulse response
h=h/1.02464; % For unity gain
end