% CI Project Phase 2

% Set the variable inputAudioName as the name of the audio file
inputAudioName = 'SPK_adult_male2';
samplingRate = 16000;

% "Main"
resampledAudio = task3(inputAudioName, samplingRate);
filters = makeBandPassBank;
lpfilter = makeLPF;
task4(resampledAudio, filters, samplingRate, lpfilter);

function filters = makeBandPassBank
    % For now, constant sampling rate and order
    Fs = 16000;
    N = 6; 
    
    % Construct FDESIGN objects
    h(1) = fdesign.bandpass('N,F3dB1,F3dB2', N, 125, 250, Fs);
    h(2) = fdesign.bandpass('N,F3dB1,F3dB2', N, 250, 500, Fs);
    h(3) = fdesign.bandpass('N,F3dB1,F3dB2', N, 500, 1000, Fs);
    h(4) = fdesign.bandpass('N,F3dB1,F3dB2', N, 1000, 2000, Fs);
    h(5) = fdesign.bandpass('N,Fp1,Fp2,Ap', 10, 2000, 4000, 1, Fs);
    h(6) = fdesign.bandpass('N,Fp1,Fp2,Ap', 10, 4000, 8000, 1, Fs);
    
    % Generate Filters - call BUTTER method
    for i=1:4
        filters(i) = design(h(i), 'butter');
    end
    
    filters(5) = design(h(5), 'cheby1');
    filters(6) = design(h(6), 'cheby1');

end

function lpfilter = makeLPF
    Fs = 16000;
    Nb   = 8;    % Numerator Order
    Na   = 8;    % Denominator Order
    F3dB = 400;  % 3-dB Frequency

    h_lpf = fdesign.lowpass('Nb,Na,F3dB', Nb, Na, F3dB, Fs);
    lpfilter = design(h_lpf, 'butter');
end

function resampledAudio = task3(fileName, fsNew)
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
    
    % 3.5. Plot  sound waveform as a func of the sample number
    %subplot(211);
    %plot(yMono);
%     title(fileName, 'Interpreter', 'none');
%     xlabel('Number of Audio Samples');
%     ylabel('Amplitude');

    % 3.6 Downsample to 16kHz
    resampledAudio = resample(yMono, fsNew, fs);
    % resampledAudio is an array of downsampled data
end

function task4(resampledAudio, filters, samplingRate, lpfilter)
    % variable of all plots
    f1 = figure;
    ch1 = filter(filters(1), resampledAudio);
    subplot(6,1,1)
    plotFFT(ch1, samplingRate);
    ch2 = filter(filters(2), resampledAudio);
    subplot(6,1,2)
    plotFFT(ch2, samplingRate);
    ch3 = filter(filters(3), resampledAudio);
    subplot(6,1,3)
    plotFFT(ch3, samplingRate);
    ch4 = filter(filters(4), resampledAudio);
    subplot(6,1,4)
    plotFFT(ch4, samplingRate);
    ch5 = filter(filters(5), resampledAudio);
    subplot(6,1,5)
    plotFFT(ch5, samplingRate);
    ch6 = filter(filters(6), resampledAudio);
    subplot(6,1,6)
    plotFFT(ch6, samplingRate);
    savefig(strcat('FFT of 6 Channels','.fig'));
   
    f2 = figure;
    abs_ch1 = abs(ch1);
    subplot(6,1,1)
    plot(abs_ch1);
    abs_ch2 = abs(ch2);
    subplot(6,1,2)
    plot(abs_ch1);
    abs_ch3 = abs(ch3);
    subplot(6,1,3)
    plot(abs_ch3);
    abs_ch4 = abs(ch4);
    subplot(6,1,4)
    plot(abs_ch4);
    abs_ch5 = abs(ch5);
    subplot(6,1,5)
    plot(abs_ch5);
    abs_ch6 = abs(ch6);
    subplot(6,1,6)
    plot(abs_ch6);
    savefig(strcat('Rectified Signals of 6 Channels','.fig'));
    
%     f3 = figure;
%     subplot(6,1,1)
%     plot(ch1);
%     subplot(6,1,2)
%     plot(ch2);
%     subplot(6,1,3)
%     plot(ch3);
%     subplot(6,1,4)
%     plot(ch4);
%     subplot(6,1,5)
%     plot(ch5);
%     subplot(6,1,6)
%     plot(ch6);
%     savefig(strcat('Signals wrt Fs of 6 Channels','.fig'));
%     
%     f4 = figure;
%     subplot(2,1,1)
%     plotFFT(ch1, samplingRate);
%     subplot(2,1,2)
%     plotFFT(ch6, samplingRate);
%     savefig(strcat('FFT Highest and Lowest Channels', '.fig'));
    
    
    % Low Pass Filter
    f5 = figure;
    lpf_ch1 = filter(lpfilter,abs_ch1);
    subplot(6,1,1)
    plot(lpf_ch1);
    
    lpf_ch2 = filter(lpfilter,abs_ch2);
    subplot(6,1,2)
    plot(lpf_ch2);
    
    lpf_ch3 = filter(lpfilter,abs_ch3);
    subplot(6,1,3)
    plot(lpf_ch3);
    
    lpf_ch4 = filter(lpfilter,abs_ch4);
    subplot(6,1,4)
    plot(lpf_ch4);
    
    lpf_ch5 = filter(lpfilter,abs_ch5);
    subplot(6,1,5)
    plot(lpf_ch5);
    
    lpf_ch6 = filter(lpfilter,abs_ch6);
    subplot(6,1,6)
    plot(lpf_ch6);
    
    % Envelopes of Lowest and Highest
    f6 = figure;
    lpf_ch1 = filter(lpfilter,abs_ch1);
    subplot(2,1,1)
    plot(ch1);
    hold on;
    plot(lpf_ch1);
    lpf_ch6 = filter(lpfilter,abs_ch6);
    subplot (2,1,2)
    plot (ch6);
    hold on;
    plot(lpf_ch6);
    savefig(strcat('Envelopes of Highest and Lowest Channels', '.fig'));
   
%     for i = 1:6
%         channelSignals(i) = filter(filters(i), resampledAudio);
%         subplot(6,1,i);
%         plotFFT(channelSignals(i), samplingRate);
%     end
end


function plotFFT(tdSignal,samplingRate)
    % PLOT FREQUENCY REPRESENTATION 
    % Get signal duration 
    [samplesRead, ~] = size(tdSignal);
    disp(samplesRead)
    L = samplesRead;
    
    % FFT two sided spectrum
    Y = fft(tdSignal);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = samplingRate*(0:(L/2))/L;
    plot(f,P1) 
    
    title('FFT')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')

end
