function ChannelState = OFDMChannelEstimation(lPreambleRX_t, overRate)

nFFT = 64;
toneSize = 52;
usedSubIndex = [2:27 39:64];

[~,  lPreambleTX_t] = OFDMPreamble(1); 
lPreambleTX_t = lPreambleTX_t(2 * 1/4* nFFT +1: end);



%% the first long preamble
lPreambleTX_t_1 = lPreambleTX_t(1: 1/2* end);
lPreambleRX_t_1 = lPreambleRX_t(1: 1/2* end);
lPreambleRX_t_1 = reshape(lPreambleRX_t_1, overRate, []).';

lPreambleTX_f_1 = fft(lPreambleTX_t_1, nFFT);
rep_lPreambleTX_f_1 = repmat(lPreambleTX_f_1, 1, overRate);
lPreambleRX_f_1 = fft(lPreambleRX_t_1, nFFT);

ChannelState_1 = zeros(nFFT, overRate);
ChannelState_1(usedSubIndex, :) = lPreambleRX_f_1(usedSubIndex, :) ./ rep_lPreambleTX_f_1(usedSubIndex, :);


%% the second long preamble
lPreambleTX_t_2 = lPreambleTX_t(1/2* end +1: end);
lPreambleRX_t_2 = lPreambleRX_t(1/2* end +1: end);
lPreambleRX_t_2 = reshape(lPreambleRX_t_2, overRate, []).';

lPreambleTX_f_2 = fft(lPreambleTX_t_2, nFFT);
rep_lPreambleTX_f_2 = repmat(lPreambleTX_f_2, 1, overRate);
lPreambleRX_f_2 = fft(lPreambleRX_t_2, nFFT);

ChannelState_2 = zeros(nFFT, overRate);
ChannelState_2(usedSubIndex, :) = lPreambleRX_f_2(usedSubIndex, :) ./ rep_lPreambleTX_f_2(usedSubIndex, :);

%%
ChannelState = 1/2 * (ChannelState_1 + ChannelState_2);
% ChannelState = ChannelState_1;

