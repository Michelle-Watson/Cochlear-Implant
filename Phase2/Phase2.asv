% CI Project Phase 2
% BME 252 - Linear Systems and Signals, Spring 2020
% Hanaan Deen, Michelle Watson, Kayley Ting

% Name of Audio File
inputAudioName = 'SPK_adult_male2_cut';
samplingRate = 16000;

% Main 
% Input audio file, convert to mono, downsample to 16kHz
resampledAudio = task3(inputAudioName, samplingRate);
% Create bank of bandpass filters
filters = makeBandPassBank;
% Create lowpass filter for envelope extraction
lpfilter = makeLPF;
% Split audio file into channels, perform envelope extraction
task4(resampledAudio, filters, samplingRate, lpfilter);

% Task 3: (Phase I)
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
%    title(fileName, 'Interpreter', 'none');
%    xlabel('Number of Audio Samples');
%    ylabel('Amplitude');

    % 3.6 Downsample to 16kHz
    resampledAudio = resample(yMono, fsNew, fs);
    %plot(resampledAudio);
    % resampledAudio is an array of downsampled data
end

% Task 4: (Phase II) Bank of Bandpass Filters
function filters = makeBandPassBank
    % Sampling Rate of Signal
    Fs = 16000;
    % Filter Orders, N1(Ch1-4), N2(Ch5-6)
    N1 = 6; 
    N2 = 10;
    
    % Construct FDESIGN objects to store filter parameters
    % Lower bounds of each channel are 1 octave apart (frequency x 2) 
    h(1) = fdesign.bandpass('N,F3dB1,F3dB2', N1, 125, 250, Fs);
    h(2) = fdesign.bandpass('N,F3dB1,F3dB2', N1, 250, 500, Fs);
    h(3) = fdesign.bandpass('N,F3dB1,F3dB2', N1, 500, 1000, Fs);
    h(4) = fdesign.bandpass('N,F3dB1,F3dB2', N1, 1000, 2000, Fs);
    h(5) = fdesign.bandpass('N,Fp1,Fp2,Ap', N2, 2000, 4000, 1, Fs);
    h(6) = fdesign.bandpass('N,Fp1,Fp2,Ap', N2, 4000, 8000, 1, Fs);
    
    % Ch1-4 Generate Butterworth Filters
    for i=1:4
        filters(i) = design(h(i), 'butter');
    end
    
    % Ch5-6 Generate Chebyshev Filters
    filters(5) = design(h(5), 'cheby1');
    filters(6) = design(h(6), 'cheby1');
end

% Tasks 5-9: (Phase II) Filter into N=6 channels, Envelope Extraction
function task4(resampledAudio, filters, samplingRate, lpfilter)

    % Task 5: Apply filters to split signal into channels
    % Generate FFT plots to verify filters performed as expected
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
    suptitle('FFT of 6 Channels')
    savefig(strcat('FFT of 6 Channels','.fig'));
  
    % Plot outputs of bandpass bank 
    f2 = figure;
    subplot(6,1,1)
    plot(ch1);
    subplot(6,1,2)
    plot(ch2);
    subplot(6,1,3)
    plot(ch3);
    subplot(6,1,4)
    plot(ch4);
    subplot(6,1,5)
    plot(ch5);
    subplot(6,1,6)
    plot(ch6);
    suptitle('Outputs of Bandpass Bank')
    
    % Task 6: Plot output signals of Ch1 and Ch6 (lowest and highest freq)
    f3 = figure;
    subplot (2,1,1)
    plot(ch1);
    title('Lowest Frequency Channel 125-250Hz');
    xlabel('Number of Audio Samples');
    ylabel('Amplitude');
    subplot (2,1,2)
    plot(ch2);
    title('Highest Frequency Channel: 4000-8000Hz');
    xlabel('Number of Audio Samples');
    ylabel('Amplitude');
    suptitle('Task6 - Highest and Lowest Channels')
    savefig(strcat('Task6 - Highest and Lowest Channels', '.fig'));
    
    % Task 7: Envelope Extraction Step 1: Rectify channels
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
    suptitle('Task8 - Envelope Extraction Step 1 (Rectified Signals)')
    savefig(strcat('Rectified Signals of 6 Channels','.fig'));
    
    % Compare frequency concentrations in different channels
    figure;
    disp('sum')
    title('Amplitude Sums of Each Channel');
    xlabel('Channel #');
    ylabel('Amplitude Sum');
    %x = [1, 2, 3, 4, 5, 6];
    labels = {'ch1'; 'ch2'; 'ch3'; 'ch4';'ch5';'ch6'};
    y = [sum(abs_ch1,1), sum(abs_ch2,1), sum(abs_ch3,1)
            sum(abs_ch4,1), sum(abs_ch5,1), sum(abs_ch6,1)];
    bar(y);
    disp(y(1));
    disp(y(2));
    disp(y(3));
    disp(y(4));
    disp(y(5));
    disp(y(6));
    
    % Task 8: Envelope Extraction Step 2: Use LPF with 400Hz cutoff
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
    suptitle('Task8 - Envelope Extraction Step 2 (Low-Pass Filtered)')
    
    % Task 9: Plot envelopes of Lowest and Highest
    f6 = figure;
    lpf_ch1 = filter(lpfilter,abs_ch1);
    subplot(2,1,1)
    %plot(ch1);     % overlaying with original signal
    %hold on;
    plot(lpf_ch1);
    title('Lowest Frequency Channel: 125-250Hz');
    xlabel('Number of Audio Samples');
    ylabel('Amplitude');
    lpf_ch6 = filter(lpfilter,abs_ch6);
    subplot (2,1,2)
    %plot (ch6);
    %hold on;
    plot(lpf_ch6);
    title('Highest Frequency Channel: 4000-8000Hz');
    xlabel('Number of Audio Samples');
    ylabel('Amplitude');
    suptitle('Task9 - Envelopes of Highest and Lowest Channels')
    savefig(strcat('Task9 - Envelopes of Highest and Lowest Channels', '.fig')); 
end

% Helper Function: LPF for Envelope Extraction
function lpfilter = makeLPF
    Fs = 16000;
    Nb   = 8;    % Numerator Order
    Na   = 8;    % Denominator Order
    F3dB = 400;  % 3-dB Frequency

    h_lpf = fdesign.lowpass('Nb,Na,F3dB', Nb, Na, F3dB, Fs);
    lpfilter = design(h_lpf, 'butter');
end

% Helper Function: Plot frequency distribution
function plotFFT(tdSignal,samplingRate)
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
