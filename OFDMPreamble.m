function [sPreamble_t, lPreamble_t] = OFDMPreamble(overRate)

%% OFDM Parameter
fsMHz = 20;
Ts = 1/fsMHz;
nFFTSize = 64;
toneSize = 52;
% usedSubIndex = [2:27 39:64];    % central zero is DC index 1
usedSubIndex = [7:32 34:59];    % central zero is DC index 33
T = Ts * nFFTSize;

%% Oversample Param
Tcp = 1/4 * nFFTSize * overRate; 
fsMHz = 20;
fosMHz = overRate * fsMHz;

t = (0 : 1/fosMHz : T - 1/fosMHz).';       % OFDM symbol duration
FFTindex = -nFFTSize/2 : nFFTSize/2 - 1;        % subcarriers' index
% FFTindex = 0 : nFFTSize - 1;

%% Short Preamble
shortTrainingSymbol = sqrt(13/6)* ...
            [0 0  1+1j 0 0 0 -1-1j 0 0 0  1+1j 0 0 0 -1-1j 0 0 0 -1-1j 0 0 0  1+1j 0 0 0 ...    % subcarriers -32 : -1
             0 0 0 -1-1j 0 0 0 -1-1j 0 0 0  1+1j 0 0 0  1+1j 0 0 0  1+1j 0 0 0  1+1j 0 0].';    % subcarriers 0 : 31
sPreamble_f = zeros(nFFTSize, 1);
sPreamble_f(usedSubIndex) = shortTrainingSymbol;
sPreamble_t = 1/sqrt(nFFTSize) * exp(1j*2*pi*t*FFTindex/T) * sPreamble_f;

% a = fft(sPreamble_t);
% for k = 1: length(a)
%     if abs(a(k)) <0.01
%         a(k) = 0;
%     end
% end
% figure; stem(angle(a));

sPreamble_t = [sPreamble_t; sPreamble_t; sPreamble_t(1:end/2)];

%% Long Preamble
longTrainingSymbol = [1  1 -1 -1  1  1 -1  1 -1  1  1  1  1  1  1 -1 -1  1  1 -1  1 -1  1  1  1  1 ...      % subcarriers -32 : -1
                      1 -1 -1  1  1 -1  1 -1  1 -1 -1 -1 -1 -1  1  1 -1 -1  1 -1  1 -1  1  1  1  1].';      % subcarriers 0 : 31
lPreamble_f = zeros(nFFTSize, 1);
lPreamble_f(usedSubIndex) = longTrainingSymbol;
lPreamble_t = 1/sqrt(nFFTSize) * exp(1j*2*pi*t*FFTindex/T) * lPreamble_f;
lPreamble_t = [lPreamble_t(end- Tcp+ 1 : end); lPreamble_t(end- Tcp+ 1 : end); lPreamble_t; lPreamble_t];