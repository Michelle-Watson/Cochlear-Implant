function Hd = IIR-butterworth
%IIR-BUTTERWORTH Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.8 and Signal Processing Toolbox 8.4.
% Generated on: 01-Aug-2020 11:11:23

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N   = 10;    % Order
Fc1 = 4000;  % First Cutoff Frequency
Fc2 = 8000;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

% [EOF]
