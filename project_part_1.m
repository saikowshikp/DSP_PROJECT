x = load("speech_with_beeps.txt");
y = remove_beeps(x, 8000);

% If Signal is Real:
audiowrite('Part_i_original.wav',x,8000);
audiowrite('Part_i_fixed.wav',y,8000);

function y = remove_beeps(x, fs)
    % Creating a Hamming Window to calculate STFT every windowSize samples
    windowSize = fs * 0.02 * 2;
    
    % Padding with zeros for window purpose
    x = [zeros(1,floor(windowSize/2)), x, zeros(1, floor(windowSize/2))];

    overlap = floor(windowSize/2);  
    segments = floor((length(x) - windowSize) / overlap) + 1;
    stft = zeros(windowSize, segments);
    power = zeros(windowSize, segments);
    hammingWindow = 0.54 - 0.46 * cos(2 * pi * (1:windowSize) / (windowSize - 1));
    
    sustainedSegments = 25; 
    threshold = zeros(1, segments);
    y = zeros((segments - 1) * overlap + windowSize, 1);

    % Computing the STFT 
    for k = 1:segments
        startSTFT = (k - 1) * overlap + 1;
        endSTFT = startSTFT + windowSize - 1;
        segment = x(startSTFT:endSTFT);
    
        windowedSegment = segment .* hammingWindow(1:windowSize);
        stft(:, k) = fft(windowedSegment);

        % Calculation of Power and Peak Detection Algorithm
        power(:, k) = 20*log(abs(stft(:, k)));
        threshold(k) = mean(power(:,k)) + std(power(:,k));          
        for l = 1:length(power(:,k))
            if(power(l,k) >= threshold(k))
                flag = 0;
                for m = 0:sustainedSegments
                    if(k-m >= 1 && power(l,k-m) >= threshold(k-m))
                        flag = flag + 1; 
                    end
                end
                if(flag == sustainedSegments + 1)
                    for m = 0:sustainedSegments
                        % Approximate Notch Filter
                        notchFilter = ones(1,320);
                        notchFilter(max(l-1,1)) = 0.5;
                        notchFilter(min(l+1,windowSize)) = 0.5;
                        notchFilter(l) = 0;
                        stft(:, k-m) = stft(:, k-m) .* notchFilter(:);
                    end 
                end
            end    
        end
    end
    for k = 1:segments
        % Inverse STFT calculation
        startISTFT = (k - 1) * overlap + 1;
        endISTFT = startISTFT + windowSize - 1;
        windowedSegment = ifft(stft(:, k), 'symmetric');    
        % Overlap-Add: As sum of cosines are constant
        y(startISTFT:endISTFT) = y(startISTFT:endISTFT) + windowedSegment;
    end
    y = transpose(y);
    y = y(floor(windowSize/2) + 1:length(x) - floor(windowSize/2));
end