function payloadTX_t = OFDMPayload(ModDataTX, NumSymbol, NumFrame, overRate)

%% OFDM Parameter
fsMHz = 20;
Ts = 1/fsMHz;
nFFTSize = 64;
cpSize = nFFTSize / 4;
toneSize = 52;
usedSubIndex = [7:32 34:59];    % central zero is DC index 33
T = Ts * nFFTSize;              % OFDM symbol legnth

%% oversampling Para
fosMHz = overRate * fsMHz; 
Tcp = cpSize * overRate;

t = (0 : 1/fosMHz : T-1/fosMHz).';  % OFDM symbol duration
FFTindex = -nFFTSize/2 : nFFTSize/2-1;     % subcarriers' index

%% payloadTX
ModDataTX = ifftshift(reshape(ModDataTX, toneSize, NumSymbol * NumFrame), 1);

payloadTX_f = zeros(nFFTSize, NumSymbol * NumFrame);
payloadTX_f(usedSubIndex, :) = reshape(ModDataTX, toneSize, NumSymbol * NumFrame);

payloadTX_t = 1/sqrt(nFFTSize) * exp(1j*2*pi*t*FFTindex/T) * payloadTX_f;
payloadTX_t = [payloadTX_t(end-Tcp+1 : end, :); payloadTX_t];
payloadTX_t = reshape(payloadTX_t, size(payloadTX_t, 1) * NumFrame * NumSymbol ,1);

