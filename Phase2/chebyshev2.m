function Hd = chebyshev2
%CHEBYSHEV2 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.8 and Signal Processing Toolbox 8.4.
% Generated on: 13-Jul-2020 21:34:18

% Chebyshev Type I Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N      = 6;     % Order
Fpass1 = 4000;  % First Passband Frequency
Fpass2 = 8000;  % Second Passband Frequency
Apass  = 1;     % Passband Ripple (dB)

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandpass('N,Fp1,Fp2,Ap', N, Fpass1, Fpass2, Apass, Fs);
Hd = design(h, 'cheby1');

% [EOF]