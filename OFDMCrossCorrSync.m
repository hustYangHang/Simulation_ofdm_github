function [SyncResult, lPreambleRX, SyncPeakIndex] = OFDMCrossCorrSync(frameRX_t, overRate)

nFFT = 64;
lPreambleLen = 128; 
cpSize = 1 / 4 * nFFT; 

[~,  lPreambleTX] = OFDMPreamble(overRate); 
lPreambleTX = lPreambleTX(2*cpSize*overRate+1: end);
LTF_Len = 128; 

SyncSize = LTF_Len * overRate; 
SyncResult = zeros(size(frameRX_t));
if overRate > 8
    output = 'Cross Correlation error, overRate is less than 8. '
    SyncPeakIndex = 0; 
    return;
end
% figure; plot(abs(frameRX_t));
for index = SyncSize : size(frameRX_t, 1)
    SyncResult(index, :) = lPreambleTX' * frameRX_t(index - SyncSize + 1: index, :);
end
SyncResult = SyncResult ./ (lPreambleTX' * lPreambleTX);    % Normalized to one

[~, SyncPeakIndex] = max(abs(SyncResult)); 
SyncPeakIndex = floor(SyncPeakIndex/overRate) * overRate; 
lPreambleRX = frameRX_t(SyncPeakIndex - lPreambleLen * overRate + 1: SyncPeakIndex, :); 

%{
figure; plot(abs(frameRX_t));
figure;  plot(abs(SyncResult));

%}