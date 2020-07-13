% CI Project Phase 1

% Set the variable inputAudioName as the name of the audio file
inputAudioName = 'violin-C4';
samplingRate = 16000;

task3(inputAudioName, samplingRate);

function task3(fileName, fsNew)
    % 3.1 Read audio file
    origAudio = strcat(fileName,'.wav');
    [origData,fs] = audioread(origAudio);
    % data = arr of sampled data
    % fs   = sampling rate
    
    % 3.2 stereo -> mono
    [~, n] = size(origData);
    % ~ = number of audio samples read
    % n = number of audio channels

    if n == 2
        % sum 2 columns (stereo) to make it 1 channel (mono).
        y = origData(:, 1) + origData(:, 2); %sum(y, 2)
        peakAmp = max(abs(y)); 
        y = y/peakAmp;
        %  check the L/R channels for orig. peak Amplitudes
        peakL = max(abs(origData(:, 1)));
        peakR = max(abs(origData(:, 2))); 
        maxPeak = max([peakL peakR]);
        %apply original peak amplitude to the normalized mono mixdown 
        yMono = y*maxPeak;
    else
        yMono = origData; 
    end

    % 3.3. Play the sound in Matlab
    % sound(yMono, fs);

    % 3.4. Write the sound to a new file.
    monoAudio = strcat(fileName, '_mono', '.wav');
    audiowrite(monoAudio,yMono,fs);

    % variable of all plots
    f1 = figure;
    
    % 3.5. Plot  sound waveform as a func of the sample number
    %subplot(211);
    %plot(yMono);
    title(fileName, 'Interpreter', 'none');
    xlabel('Number of Audio Samples');
    ylabel('Amplitude');

    % 3.6 Downsample to 16kHz
    resampledAudio = resample(yMono, fsNew, fs);
    % resampledAudio is an array of downsampled data
    
    % PLOT FREQUENCY REPRESENTATION 
    % Get signal duration 
    [samplesRead, ~] = size(resampledAudio);
    disp(samplesRead)
    L = samplesRead;
    % FFT two sided spectrum? 
    Y = fft(resampledAudio);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = fsNew*(0:(L/2))/L;
    subplot(2,1,1)
    plot(f,P1) 
    
    title(strcat(fileName,'Single-Sided Amplitude Spectrum of X(t)'), 'Interpreter', 'none')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    savefig(strcat(fileName, '.fig'));
    saveas(f1,strcat(fileName, '.png'));
    
    filteredSignal = filter(Hd,resampledAudio);
    filteredFFT = fft(filteredSignal);
    P4 = abs(filteredFFT/L);
    P3 = P4(1:L/2+1);
    P3(2:end-1) = 2*P3(2:end-1);
    
    f = fsNew*(0:(L/2))/L;
    subplot(2,1,1)
    plot(f,P3) 
    
    % 3.5 (redo) after resampling, not necessary
%     subplot(312);
%     plot(resampledAudio);
%     title('Resampled Audio');
%     xlabel('Number of Audio Samples');
%     ylabel('Amplitude');
    
    % 3.7 plot
    % subplot: row col plot #
%     subplot(212)
%     dt = 1/fsNew; % pass in new FS
%     F = 1000;
%     T = 1/F;
%     [samplesRead, ~] = size(resampledAudio);
%     duration = samplesRead/fsNew; % length of resampled audio (sec) = samplesRead/samplingRate
%     t = (0:dt:duration); 
%     sig = cos(2*pi*F*t);
%     plot(t, sig);
%     xlim([0 2*T]); % display 2 cycles, or 2 periods
%     title('2 cycles of Waveform');
%     xlabel('Time (seconds)');
%     ylabel('Amplitude')
%     grid; grid minor;
%     
%     % 3.7 sound
%     pause(duration);
%     sound(sig, fsNew);
%     
%     % save figure as .png + .fig
%     saveas(f1,strcat(fileName, '.png'));
%     savefig(strcat(fileName, '.fig'));
end

