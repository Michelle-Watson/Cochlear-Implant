% CI Project Phase 3
% BME 252 - Linear Systems and Signals, Spring 2020
% Hanaan Deen, Michelle Watson, Kayley Ting

% Name of Audio File
inputAudioName = '1.5b_justwhatIneeded-sarcasm';
samplingRate = 16000;

% Create cell arrays, accessible from different contexts
channels = {1,2,3,4,5,6};
rect_channels = {1,2,3,4,5,6};
env_channels = {1,2,3,4,5,6};

test = logspace(2,4,6);
disp(test);

% Main 
% (Phase I) Input audio file, convert to mono, downsample to 16kHz
resampledAudio = preprocess(inputAudioName, samplingRate);

% ITERATION: *** CHANGE PARAMETERS FOR TESTING HERE ***
% (Phase III) Iteration: set desired parameters (eg. octave vs. linear subband)
% Note: For version 1, set "OCTAVE" and 400
subBandType = "OCTAVE";
lpf_freq =400;
outputFileName = strcat('A_', inputAudioName,'.wav');

% (Phase II) Create bank of bandpass filters
[lowerfreqs, upperfreqs] = setCutoffs(subBandType, 6);
disp("CHECK THE FREQUENCIES");
disp(lowerfreqs);
disp(upperfreqs);

filters = makeBandPassBank(lowerfreqs, upperfreqs);
% Create lowpass filter for envelope extraction
lpfilter = makeLPF(lpf_freq);
% Split audio file into channels, perform envelope extraction
[channels, rect_channels, env_channels] = makeChannels(resampledAudio, filters, samplingRate, lpfilter);

% (Phase III) Synthesize audio output
cosines = genCosines(env_channels,samplingRate);
sound_output = generateSound(env_channels, cosines, outputFileName);
plot(sound_output);
title('Synthesized Output');

checkFreqDomain(channels, samplingRate, lowerfreqs, upperfreqs);

% Task 3: (Phase I)
function resampledAudio = preprocess(fileName, fsNew)
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
    % monoAudio = strcat(fileName, '_mono', '.wav');
    % audiowrite(monoAudio,yMono,fs);
    
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
function filters = makeBandPassBank(lowerfreqs, upperfreqs)
    % Sampling Rate of Signal
    Fs = 16000;
    % Filter Orders, N1(Ch1-4), N2(Ch5-6)
    N1 = 6; 
    N2 = 10;
   
    % Construct FDESIGN objects to store filter parameters. 
    % Uses vectors containing the desired frequency bounds
    % Note Ch5 and Ch6 parameters are different because they are chebyshev filters
    for i=1:4
        h(i) = fdesign.bandpass('N,F3dB1,F3dB2', N1, lowerfreqs(i), upperfreqs(i), Fs);
    end
    h(5) = fdesign.bandpass('N,Fp1,Fp2,Ap', N2, lowerfreqs(5), upperfreqs(5), 1, Fs);
    h(6) = fdesign.bandpass('N,Fp1,Fp2,Ap', N2, lowerfreqs(6), upperfreqs(6), 1, Fs);
    
    % Ch1-4 Generate Butterworth Filters
    for i=1:4
        filters(i) = design(h(i), 'butter');
    end
    
    % Ch5-6 Generate Chebyshev Filters
    filters(5) = design(h(5), 'cheby1');
    filters(6) = design(h(6), 'cheby1');
end

% Tasks 5-9: (Phase II) Filter into N=6 channels, Envelope Extraction
function [channels, rect_channels, env_channels] = makeChannels(resampledAudio, filters, samplingRate, lpfilter)
    for i=1:6
       % Task 5: Apply filters to split signal into channels
       channels{1,i} = filter(filters(i), resampledAudio);
       % Task 7: Envelope Extraction Step 1: Rectify channels
       rect_channels{1,i} = abs(channels{1,i});
       % Task 8: Envelope Extraction Step 2: Use LPF with 400Hz cutoff
       env_channels{1,i} = filter(lpfilter,rect_channels{1,i});
    end
    
    % Plot channel envelopes (Time Domain) 
    figure;
    for i=1:6
        subplot(6,1,i);
        plot(env_channels{1,i});
        title(['Envelope of ch' num2str(i)]);
    end      
end

% Tasks 10-12: (Phase III) Amplitude Modulation, Sum to create output sound
function sound_output = generateSound(env_channels, cosines, outputFileName)
    modAmp_channels = {1,2,3,4,5,6};
    
    % Amplitude Modulation
    for i=1:6
        % Pad the envelope with a 0, not sure why length diff
        chEnv = env_channels{1,i};
        chEnv(end+1)=0;
        cosine = transpose(cosines{1,i});
         % AmplitudeModulate
        modAmp_channels{1,i} = chEnv.*cosine;
    end

    % Plot Amplitude Modulation to Check
%     for i=1:6
%         figure;
%         subplot(2,1,1);
%         plot(env_channels{1,i});
%         title('Envelope');
% 
%         subplot(2,1,2);
%         plot(modAmp_channels{1,i});
%         title('Amplitude Modulated Signal');
%     end
    
    sound_output = modAmp_channels{1,i};
    for i=2:6
       sound_output = sound_output + modAmp_channels{1,i};
    end
    
    sound_output(end)=[];
    figure;
    plot(sound_output);
    title('Sound Output');
    sound(sound_output,16000);
    outputFileName = strcat(outputFileName, '.wav');
    audiowrite(outputFileName,sound_output,16000);
end

% Helper Function: Set cutoff bounds for filter bands
function [lowerfreqs, upperfreqs] = setCutoffs(type, numChannels)
    % Default, first lower bound is 100 + last upper bound is 8000
    lowerfreqs(1) = 100;
    upperfreqs(numChannels) = 8000;
    
    % Option 1: Octave Spacing
    if isequal(type, 'OCTAVE')
       lowerfreqs(1) = 100;  upperfreqs(1) = 250;
       lowerfreqs(2) = 250;  upperfreqs(2) = 500;
       lowerfreqs(3) = 500;  upperfreqs(3) = 1000;
       lowerfreqs(4) = 1000; upperfreqs(4) = 2000;
       lowerfreqs(5) = 2000; upperfreqs(5) = 4000;
       lowerfreqs(6) = 4000; upperfreqs(6) = 7999;
    end 
    
    % Option 2: Even linear spacing
    if isequal(type, 'LINEAR')
       lowerfreqs(1) = 100;  upperfreqs(1) = 1500;
       lowerfreqs(2) = 1500;  upperfreqs(2) = 2900;
       lowerfreqs(3) = 2900;  upperfreqs(3) = 4300;
       lowerfreqs(4) = 4300; upperfreqs(4) = 5700;
       lowerfreqs(5) = 5700; upperfreqs(5) = 7100;
       lowerfreqs(6) = 7100; upperfreqs(6) = 7999;
%         spacing = (8000-100)/numChannels;
%         for i=2:6
%             lowerfreqs(i) = lowerfreqs(i-1) + spacing;
%         end 
%         for i=1:6
%             upperfreqs(i) = 100 + (i*spacing);
%         end
    end
    
    % Option 3: Logarithmic Spacing of Corners
    if isequal(type, 'LOGCENTRE') % 100, 251, 631, 1585, 3981, *8000
       lowerfreqs(1) = 50;  upperfreqs(1) = 200; 
       lowerfreqs(2) = 200;  upperfreqs(2) = 315;
       lowerfreqs(3) = 315;  upperfreqs(3) = 1260;
       lowerfreqs(4) = 1260; upperfreqs(4) = 2000;
       lowerfreqs(5) = 2000; upperfreqs(5) = 7924;
       lowerfreqs(6) = 7924; upperfreqs(6) = 7999;
    end  
    
    % Option 4: Logarithmic Spacing of Corners
    if isequal(type, 'LOGARITHMIC')
       lowerfreqs(1) = 100;  upperfreqs(1) = 200;
       lowerfreqs(2) = 200;  upperfreqs(2) = 450;
       lowerfreqs(3) = 450;  upperfreqs(3) = 1000;
       lowerfreqs(4) = 1000; upperfreqs(4) = 2150;
       lowerfreqs(5) = 2150; upperfreqs(5) = 4650;
       lowerfreqs(6) = 4650; upperfreqs(6) = 7999;
    end  
    
    % Option 5: Overlapping bands
    if isequal(type, 'OVERLAP')
       lowerfreqs(1) = 100;  upperfreqs(1) = 400;
       lowerfreqs(2) = 350;  upperfreqs(2) = 600;
       lowerfreqs(3) = 500;  upperfreqs(3) = 1000;
       lowerfreqs(4) = 250; upperfreqs(4) = 6500;
       lowerfreqs(5) = 4000; upperfreqs(5) = 6000;
       lowerfreqs(6) = 4000; upperfreqs(6) = 7999;
    end  
end

% Helper Function: Make Cosines
function cosines = genCosines(env_channels, samplingRate)
    % Cosine Parameters
    dt = 1/samplingRate;
    [samplesRead, ~] = size(env_channels{1,1});
    duration = samplesRead/samplingRate; % length of resampled audio (sec) = samplesRead/samplingRate
    t = (0:dt:duration); 
    
    % Centre Frequencies
    fc(1) = sqrt(125*250);
    fc(2) = sqrt(250*500);
    fc(3) = sqrt(500*1000);
    fc(4) = sqrt(1000*2000);
    fc(5) = sqrt(2000*4000);
    fc(6) = sqrt(4000*8000);

    %Plot Cosines to Check
    figure;
    for i=1:6
       subplot(6,1,i);
       cosines{1,i} = cos(2*pi*fc(i)*t);
       plot(t, cosines{1,i})
       T = 1/200; % arbritrary period to display for
       %T = 1/fc(i);     % length of one period
       xlim([0 2*T]); % play 2 cycles, or 2 periods
       title(['2 cycles of cos Waveform' num2str(i)]);
    end
end

% Helper Function: LPF for Envelope Extraction
function lpfilter = makeLPF(lpf_freq)
    Fs = 16000;
    Nb   = 8;    % Numerator Order
    Na   = 8;    % Denominator Order
    F3dB = lpf_freq;  % 3-dB Frequency

    h_lpf = fdesign.lowpass('Nb,Na,F3dB', Nb, Na, F3dB, Fs);
    lpfilter = design(h_lpf, 'butter');
end

% Helper Function: Plot frequency distribution
function plotFFT(tdSignal,samplingRate)
    % Get signal duration 
    [samplesRead, ~] = size(tdSignal);
    %disp(samplesRead)
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

function checkFreqDomain(channels, samplingRate, lowerfreqs, upperfreqs)
    figure;
    for i=1:6
        subplot(6,1,i);
        plotFFT(channels{1,i}, samplingRate);
        plotName = strcat('FFT', 'Ch', num2str(i), ' (', num2str(lowerfreqs(i)), '-', num2str(upperfreqs(i)), ')');
        title(plotName);
    end
end

