function BER_1 = cal_through_bits(snr_data,csi_fre_data)
%% pack a ofdm symbol and implement modulation and demodulation 
% clear all;
% clc;

%% OFDM Parameter
fsMHz = 20;
Ts = 1/fsMHz;
nFFTSize = 64;
cpSize = nFFTSize/4;
toneSize = 52;
usedSubIndex = [2:27 39:64];    % central zero is DC index 1
T = Ts * nFFTSize;              % The length of OFDM symbol

%% oversampling Para
overRate = 1;    % oversample rate 

%% OFDM Preamble
[sPreambleTX_t, lPreambleTX_t] = OFDMPreamble(overRate); 

%% OFDM Payload
NumSymbol = 30;%6的倍数
NumFrame = 1;
NumData = toneSize * NumSymbol * NumFrame;
mod_order = 8;%1/bpsk-1/2 2/qpsk-1/2 3/16qam-1/2 4/64qam-2/3 5/bpsk-3/4 6/qpsk-3/4 7/16qam-3/4 8/64qam-3/4 
BER_1 = ones(8,1);
for method = 1:mod_order
% method = 4;
switch method
    case 1%BPSK 1/2
%         rawDataTX = randi(2, NumData, 1) - 1;
        rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 1);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);%扰码
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
%         rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(52*ceil(length(rawDataBinTX_conv)/52)-length(rawDataBinTX_conv),1)];
      
        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv,1);%交织
        interleaving_bit_length = length(rawDataBinTX_inter1);
%         
        rawDataBinTX_inter = [rawDataBinTX_inter1;zeros(52*(ceil(length(rawDataBinTX_inter1)/52))-length(rawDataBinTX_inter1),1)];
        
        ModDataTX = 1/sqrt(2)*step(comm.BPSKModulator, rawDataBinTX_inter);
    case 2%QPSK 1/2
        rawDataTX = randi(4, NumData, 1) - 1;
%         rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 2);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
%         rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(2*52*ceil(length(rawDataBinTX_conv)/(2*52))-length(rawDataBinTX_conv),1)];

        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv,1);%交织
        interleaving_bit_length = length(rawDataBinTX_inter1);
        
        rawDataBinTX_inter = [rawDataBinTX_inter1;zeros(2*52*(ceil(length(rawDataBinTX_inter1)/(2*52)))-length(rawDataBinTX_inter1),1)];        

        ModDataTX =1/sqrt(4)*step(comm.QPSKModulator('BitInput', true), rawDataBinTX_inter);
    case 3%16-QAM 1/2
%         rawDataTX = randi(16, NumData, 1) - 1;
        rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 4);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
%         rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(4*52*ceil(length(rawDataBinTX_conv)/(4*52))-length(rawDataBinTX_conv),1)];
 
        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv,1);%交织

        interleaving_bit_length = length(rawDataBinTX_inter1);
        
        rawDataBinTX_inter = [rawDataBinTX_inter1;rawDataBinTX_inter1;zeros(4*52*(2*ceil(length(rawDataBinTX_inter1)/(4*52)))-2*length(rawDataBinTX_inter1),1)];
        
        ModDataTX = 1/sqrt(10)*step(comm.RectangularQAMModulator('BitInput', true), rawDataBinTX_inter);
    case 4%64-QAM 2/3
%         rawDataTX = randi(64, NumData, 1) - 1;
        rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 6);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
        rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(6*52*ceil(length(rawDataBinTX_conv)/(6*52))-length(rawDataBinTX_conv),1)];

        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv1,1);%交织
        interleaving_bit_length = length(rawDataBinTX_inter1);    
        rawDataBinTX_inter = [rawDataBinTX_inter1;rawDataBinTX_inter1;zeros(6*52*(ceil(2*length(rawDataBinTX_inter1)/(6*52)))-2*length(rawDataBinTX_inter1),1)];

        
        ModDataTX = 1/sqrt(43)*step(comm.RectangularQAMModulator('ModulationOrder', 64, 'BitInput', true) , rawDataBinTX_inter);
    case 5%BPSK 3/4
%         rawDataTX = randi(2, NumData, 1) - 1;
        rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 1);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);%扰码
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
%         rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(52*ceil(length(rawDataBinTX_conv)/52)-length(rawDataBinTX_conv),1)];

        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv,1);%交织
        interleaving_bit_length = length(rawDataBinTX_inter1);      
        rawDataBinTX_inter = [rawDataBinTX_inter1;rawDataBinTX_inter1;zeros(52*(ceil(2*length(rawDataBinTX_inter1)/52))-2*length(rawDataBinTX_inter1),1)];
        
        ModDataTX = 1/sqrt(2)*step(comm.BPSKModulator, rawDataBinTX_inter);
    case 6%QPSK 3/4
%         rawDataTX = randi(4, NumData, 1) - 1;
        rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 2);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
%         rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(2*52*ceil(length(rawDataBinTX_conv)/(2*52))-length(rawDataBinTX_conv),1)];
  
        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv,1);%交织
        interleaving_bit_length = length(rawDataBinTX_inter1);   
        rawDataBinTX_inter = [rawDataBinTX_inter1;rawDataBinTX_inter1;zeros(2*52*(ceil(2*length(rawDataBinTX_inter1)/(2*52)))-2*length(rawDataBinTX_inter1),1)];        
        
        ModDataTX =1/sqrt(4)*step(comm.QPSKModulator('BitInput', true), rawDataBinTX_inter);
    case 7%16-QAM 3/4
%         rawDataTX = randi(16, NumData, 1) - 1;
        rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 4);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
%         rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(4*52*ceil(length(rawDataBinTX_conv)/(4*52))-length(rawDataBinTX_conv),1)];

        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv,1);%交织
        interleaving_bit_length = length(rawDataBinTX_inter1);     
        rawDataBinTX_inter = [rawDataBinTX_inter1;rawDataBinTX_inter1;zeros(4*52*(ceil(2*length(rawDataBinTX_inter1)/(4*52)))-2*length(rawDataBinTX_inter1),1)];
        
        ModDataTX = 1/sqrt(10)*step(comm.RectangularQAMModulator('BitInput', true), rawDataBinTX_inter);
    case 8%64-QAM 3/4
%         rawDataTX = randi(64, NumData, 1) - 1;
        rawDataTX = zeros(NumData,1);
        rawDataBinTX = Dec2BinVector(rawDataTX, 6);
        rawDataBinTX_scr = scramble(rawDataBinTX,1);
        
        rawDataBinTX_conv = tx_conv_encoder(method,rawDataBinTX_scr);%卷积
        conv_bit_length = length(rawDataBinTX_conv);
%         rawDataBinTX_conv1 = [rawDataBinTX_conv;zeros(6*52*ceil(length(rawDataBinTX_conv)/(6*52))-length(rawDataBinTX_conv),1)];

        rawDataBinTX_inter1 = interleaving(rawDataBinTX_conv,1);%交织
        interleaving_bit_length = length(rawDataBinTX_inter1);    
        rawDataBinTX_inter = [rawDataBinTX_inter1;rawDataBinTX_inter1;zeros(6*52*(ceil(2*length(rawDataBinTX_inter1)/(6*52)))-2*length(rawDataBinTX_inter1),1)];

        
        ModDataTX = 1/sqrt(43)*step(comm.RectangularQAMModulator('ModulationOrder', 64, 'BitInput', true) , rawDataBinTX_inter);

end
switch method
    case 1%BPSK 1/2
        NumSymbol1 = length(rawDataBinTX_inter)/52;
    case 2%QPSK 1/2
        NumSymbol1 = length(rawDataBinTX_inter)/(2*52);
    case 3%16-QAM 1/2
        NumSymbol1 = length(rawDataBinTX_inter)/(4*52);
    case 4%64-QAM 2/3
        NumSymbol1 = length(rawDataBinTX_inter)/(6*52);
    case 5%BPSK 3/4
        NumSymbol1 = length(rawDataBinTX_inter)/52;
    case 6%QPSK 3/4
        NumSymbol1 = length(rawDataBinTX_inter)/(2*52);
    case 7%16-QAM 3/4
        NumSymbol1 = length(rawDataBinTX_inter)/(4*52);
    case 8%64-QAM 3/4
        NumSymbol1 = length(rawDataBinTX_inter)/(6*52); 
end

% NumSymbol1 = NumSymbol;
payloadTX_t = OFDMPayload(ModDataTX, NumSymbol1, NumFrame, overRate);
% figure(method);
% plot(abs(reshape(payloadTX_t,numel(payloadTX_t),1)));
%% frame
% frameTX_t = zeros(32560,1,'gpuArray');
frameTX_t = [repmat(sPreambleTX_t, 1, NumFrame); ...	NumFrame frames
    repmat(lPreambleTX_t, 1, NumFrame); ...
    reshape(payloadTX_t, [], NumFrame)];

frameLen = length(frameTX_t);

%% add frequency offset
% FreqOffset = 0.0001;
% 
% frameTX_t = reshape(frameTX_t, [], 1);
% for k = 1: length(frameTX_t)
%     frameTX_t(k) = frameTX_t(k) * exp(1j * (k-1) * FreqOffset);
% end
% frameTX_t = reshape(frameTX_t, [], NumFrame);


%% Add channel and noise
csi_fre_data_choose(:,1) = csi_fre_data;
for i = 5:length(frameTX_t)/80-1
    frameTX_t(17+80*(i-1):80*i,1) = fft(frameTX_t(17+80*(i-1):80*i,1));
    frameTX_t(17+80*(i-1):80*i,1) = frameTX_t(17+80*(i-1):80*i,1).*csi_fre_data_choose(:,1);
    frameTX_t(17+80*(i-1):80*i,1) = ifft(frameTX_t(17+80*(i-1):80*i,1));
end
for i = 1:2
    frameTX_t(160+32+1+64*(i-1):160+32+64*i,1) = fft(frameTX_t(160+32+1+64*(i-1):160+32+64*i,1));
    frameTX_t(160+32+1+64*(i-1):160+32+64*i,1) = frameTX_t(160+32+1+64*(i-1):160+32+64*i,1).*csi_fre_data_choose(:,1);
    frameTX_t(160+32+1+64*(i-1):160+32+64*i,1) = ifft(frameTX_t(160+32+1+64*(i-1):160+32+64*i,1));
end

% switch method
%     case 1
%         frameRX_t = awgn(frameTX_t,esnr(1,1),'measured');   
%     case 2
%         frameRX_t = awgn(frameTX_t,esnr(2,1),'measured');
%     case 3
%         frameRX_t = awgn(frameTX_t,esnr(3,1),'measured');
%     case 4
%         frameRX_t = awgn(frameTX_t,esnr(4,1),'measured');
%     case 5
%         frameRX_t = awgn(frameTX_t,esnr(1,1),'measured');
%     case 6
%         frameRX_t = awgn(frameTX_t,esnr(2,1),'measured');
%     case 7
%         frameRX_t = awgn(frameTX_t,esnr(3,1),'measured');
%     case 8
%         frameRX_t = awgn(frameTX_t,esnr(4,1),'measured');        
% end
frameRX_t = awgn(frameTX_t,snr_data,'measured');
% frameRX_t = frameTX_t;
%% %%
payloadRX_f = zeros(nFFTSize, NumSymbol * NumFrame, overRate);

for FrameIndex = 1: NumFrame

%% time sync
[result, lPreambleRX_t_x, SyncPeakIndex_x] = OFDMCrossCorrSync(frameRX_t(:, FrameIndex), overRate);    % time sync
% SyncPeakIndex_x
% figure; plot(abs(result));
% SyncPeakIndex_x = 320 * overRate + offset * overRate;
SyncPeakIndex_x = 320;
lPreambleRX_t_x = frameRX_t(SyncPeakIndex_x - 160 * overRate + 32 * overRate + 1: SyncPeakIndex_x, FrameIndex);

%% Channel Estimation
ChannelState = OFDMChannelEstimation(lPreambleRX_t_x, overRate);
% ChannelState = csi_fre_data;

%% get payload
if SyncPeakIndex_x > 320 * overRate
    payloadRX_t = [frameRX_t(SyncPeakIndex_x + 1 : end, FrameIndex); zeros(max(0, SyncPeakIndex_x - 320 * overRate), 1)];
else
    payloadRX_t = frameRX_t(SyncPeakIndex_x + 1 : end - (320 * overRate - SyncPeakIndex_x), FrameIndex);
end
payloadRX_t = reshape(payloadRX_t, (nFFTSize + cpSize) * overRate, NumSymbol1);
payloadRX_t = payloadRX_t(cpSize * overRate + 1: end, :);  % Remove cp

FFTindex = (-nFFTSize/2 : nFFTSize/2 - 1).';
for OSIndex = 1 : overRate
    payloadRX_f(:, (FrameIndex - 1)* NumSymbol1 +1: FrameIndex * NumSymbol1, OSIndex) = ... 
        1/sqrt(nFFTSize) * fft(payloadRX_t(OSIndex: overRate: end, :));    
    payloadRX_f(:, (FrameIndex - 1)* NumSymbol1 +1: FrameIndex * NumSymbol1, OSIndex) = ...
        payloadRX_f(:, (FrameIndex - 1)* NumSymbol1 +1: FrameIndex * NumSymbol1, OSIndex) ...
        ./ repmat(ChannelState(:, OSIndex), 1, NumSymbol1);
end
end
% 
% hold on
% plot(abs(reshape(payloadRX_f,numel(payloadRX_f),1)));

%% 1x; BPSK, QPSK, QAM decode
payloadRX_f_1 = payloadRX_f(:, :, 1);
ModDataRX_1 = payloadRX_f_1(usedSubIndex, 1:NumSymbol1);
% NumSymbol = NumSymbol-1;
ModDataRX_1 = reshape(ModDataRX_1, toneSize * NumSymbol1 * NumFrame, 1);
switch method
    case 1%BPSK 1/2
        rawDataRX_1 = step(comm.BPSKDemodulator, sqrt(2)*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 1);
        
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织
       
        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积
        
        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);
        BER_1(1,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*1,1) - rawDataBinTX(1:52*(NumSymbol-1)*1,1))) / (NumData-52);

%         [1;length(rawDataBinRX_1(1:52*(NumSymbol-1)*1,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*1,1) - rawDataBinTX(1:52*(NumSymbol-1)*1,1))==1)]
%% 计算正确解包概率    
%         correct_pkt = 0;
%         for i = 1:NumSymbol-1
%             if biterr(rawDataBinRX_1(52*(i-1)+1:52*i,1),rawDataBinTX(52*(i-1)+1:52*i,1)) == 0
%                 correct_pkt = correct_pkt + 1;
%             end
%         end
%         correct_poss(1,1) = correct_pkt/i;
%         
%         through_bits(1,1) = (NumData-52) *(1-BER_1);
    case 2%QPSK 1/2
        rawDataRX_1 = step(comm.QPSKDemodulator, 2*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 2);
        
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织     
        
        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积      

        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);

        BER_1(3,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*2,1) - rawDataBinTX(1:52*(NumSymbol-1)*2,1))) / ((NumData-52) * 2);
        
%         [3;length(rawDataBinRX_1(1:52*(NumSymbol-1)*2,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*2,1) - rawDataBinTX(1:52*(NumSymbol-1)*2,1))==1)]
    case 3%16-QAM 1/2
        rawDataRX_1 = step(comm.RectangularQAMDemodulator, sqrt(10)*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 4);      
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织   
        
        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积  

        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);

        BER_1(5,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*4,1) - rawDataBinTX(1:52*(NumSymbol-1)*4,1))) / ((NumData-52) * 4);
        
%         [5;length(rawDataBinRX_1(1:52*(NumSymbol-1)*4,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*4,1) - rawDataBinTX(1:52*(NumSymbol-1)*4,1))==1)]
    case 4%64-QAM 2/3
        rawDataRX_1 = step(comm.RectangularQAMDemodulator(64), sqrt(43)*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 6);
         
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织   
        
        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积       

        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);

        BER_1(7,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*6,1) - rawDataBinTX(1:52*(NumSymbol-1)*6,1))) / ((NumData-52) * 6);
%         [7;length(rawDataBinRX_1(1:52*(NumSymbol-1)*6,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*6,1) - rawDataBinTX(1:52*(NumSymbol-1)*6,1))==1)]
     case 5%BPSK 3/4
        rawDataRX_1 = step(comm.BPSKDemodulator, sqrt(2)*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 1);       
      
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织    

        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积

        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);

        BER_1(2,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*1,1) - rawDataBinTX(1:52*(NumSymbol-1)*1,1))) / (NumData-52);
%         [2;length(rawDataBinRX_1(1:52*(NumSymbol-1)*1,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*1,1) - rawDataBinTX(1:52*(NumSymbol-1)*1,1))==1)]
    case 6%QPSK 3/4
        rawDataRX_1 = step(comm.QPSKDemodulator, 2*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 2);
        
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织    

        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积   

        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);

        BER_1(4,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*2,1) - rawDataBinTX(1:52*(NumSymbol-1)*2,1))) / ((NumData-52) * 2);
%         [4;length(rawDataBinRX_1(1:52*(NumSymbol-1)*2,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*2,1) - rawDataBinTX(1:52*(NumSymbol-1)*2,1))==1)]
    case 7%16-QAM 3/4
        rawDataRX_1 = step(comm.RectangularQAMDemodulator, sqrt(10)*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 4);      
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织   

        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积

        
        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);

        BER_1(6,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*4,1) - rawDataBinTX(1:52*(NumSymbol-1)*4,1))) / ((NumData-52) * 4);
%         [6;length(rawDataBinRX_1(1:52*(NumSymbol-1)*2,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*2,1) - rawDataBinTX(1:52*(NumSymbol-1)*2,1))==1)]
    case 8%64-QAM 3/4
        rawDataRX_1 = step(comm.RectangularQAMDemodulator(64), sqrt(43)*ModDataRX_1);
        rawDataBinRX_1 = Dec2BinVector(rawDataRX_1, 6);       
        rawDataBinRX_1_deinter = rawDataBinRX_1(1:interleaving_bit_length,1);
        rawDataBinRX_2 = deinterleaving(rawDataBinRX_1_deinter,1);%解交织      
        
        rawDataBinRX_1_deconv1 = rawDataBinRX_2(1:conv_bit_length,1);
        rawDataBinRX_1_deconv = tx_conv_decoder(method,rawDataBinRX_1_deconv1);%解卷积           

        rawDataBinRX_1 = scramble(rawDataBinRX_1_deconv,2);

        BER_1(8,1) = sum(abs(rawDataBinRX_1(1:52*(NumSymbol-1)*6,1) - rawDataBinTX(1:52*(NumSymbol-1)*6,1))) / ((NumData-52) * 6);
%         [8;length(rawDataBinRX_1(1:52*(NumSymbol-1)*6,1));find((rawDataBinRX_1(1:52*(NumSymbol-1)*6,1) - rawDataBinTX(1:52*(NumSymbol-1)*6,1))==1)]
end

end

