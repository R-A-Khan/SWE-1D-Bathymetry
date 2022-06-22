function [smoothed]=grad_smooth_H2(signal,filt)
% Usage:
% smoothed = grad_smooth(signal,filt)
%
% Input:
% signal = signal to be smoothed
% fmin   = set to (1+L^2k^2) all Fourier modes greater than filt*max frequency
% (0<filt<1)
% Output:
% smoothed = smoothed sobolev gradient version of signal
 
N = numel(signal);
 
% Take fft of signal
fJ=fft(signal);
 
fs = 0:N/2; % Frequencies
 
fmin = round(filt*max(fs)); % Minimum filtering frequency (corresponding to 1/L in the notes)
fJ_H = fJ;
for k = 1:N/2+1
    fJ_H(k) = fJ(k)/(1+fs(k)^4/fmin^4);
end
 
% Use symmetry property of fft of real signals to deal with remaining Fourier modes 
% Frequencies associated with the elements 1:N/2+1 of the fft vector are 0:N/2
for k=N/2+2:N
    fJ_H(k)=conj(fJ_H(N-k+2));
end
 
% Inverse Fourier transform
smoothed=ifft(fJ_H);