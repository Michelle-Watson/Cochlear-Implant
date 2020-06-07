
inputAudioName = '1980s-Casio-Trumpet-C5';
samplingRate = 16000;

task3(inputAudioName, samplingRate);

function task3(fileName, fsNew)
    % 3.1 Read audio file
    origAudio = strcat(fileName,'.wav');
    [origData,fs] = audioread(origAudio);
    % data = arr of sampled data
    % fs   = sampling rate
    
    % 3.2 stereo -> mono
    [m, n] = size(origData);
    % m = number of audio samples read
    % n = number of audio channels

    if n == 2
        % add the two columns to make it single channel (or a 1-column array).
        y = origData(:, 1) + origData(:, 2); %sum(y, 2) also accomplishes this
        peakAmp = max(abs(y)); 
        y = y/peakAmp;
        %  check the L/R channels for orig. peak Amplitudes
        peakL = max(abs(origData(:, 1)));
        peakR = max(abs(origData(:, 2))); 
        maxPeak = max([peakL peakR]);
        %apply x's original peak amplitude to the normalized mono mixdown 
        yMono = y*maxPeak;
    else
        yMono = origData; 
    end

    % 3.3. Play the sound in Matlab
    sound(yMono, fs);

    % 3.4. Write the sound to a new file.
    monoAudio = strcat(fileName, 'mono', '.wav');
    audiowrite(monoAudio,yMono,fs);

    % variable of all plots
    f1 = figure;
    
    % 3.5. Plot  sound waveform as a func of the sample number
    subplot(211);
    plot(yMono);
    title('Original Audio');
    xlabel('Number of Audio Samples');
    ylabel('Amplitude');    
    disp(fs);

    % 3.6 Downsample to 16kHz
    resampledAudio = resample(yMono, 16000, fs);
    % resampledAudio is an array of downsampled data
    
    % 3.5 (redo) after resampling
%     subplot(312);
%     plot(resampledAudio);
%     title('Resampled Audio');
%     xlabel('Number of Audio Samples');
%     ylabel('Amplitude');
    
    % 3.7 (struggle)
    % subplot: row col plot #
    subplot(212)
    dt = 1/fsNew; % pass in new FS
    F = 1000;
    T = 1/F;
    duration = m/fs; % legnth of original audio (sec)
    t = (0:dt:duration); 
    sig = cos(2*pi*F*t);
    plot(t, sig);
    xlim([0 2*T]); % display 2 cycles, or 2 periods
    title('2 cycles of Waveform');
    xlabel('Time (seconds)');
    ylabel('Amplitude')
    grid; grid minor;
    
    % save figure as .png + .fig
    saveas(f1,strcat(fileName, '.png'));
    savefig(strcat(fileName, '.fig'));
end


