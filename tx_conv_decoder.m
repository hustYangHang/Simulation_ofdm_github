function decode_bit = tx_conv_decoder(MOD_ORDER,codeword)
%% ½â¾í»ý±àÂë
switch MOD_ORDER
    case 1%BPSK 1/2
        traceBack = 1;
        tPoly = poly2trellis(7,[133,171]);
        decDelay = traceBack;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
    case 2%QPSK 1/2
        traceBack = 1;
        tPoly = poly2trellis(7,[133,171]);
        decDelay = traceBack;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
    case 3%16-QAM 1/2
        traceBack = 1;
        tPoly = poly2trellis(7,[133,171]);
        decDelay = traceBack;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
    case 4%64-QAM 2/3
        traceBack = 2;
        tPoly = poly2trellis([5 4],[23 35 0; 0 5 13]);
        decDelay = traceBack*2;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
    case 5%BPSK 3/4
        traceBack = 3;
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        decDelay = traceBack*3;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
    case 6%QPSK 3/4
        traceBack = 3;
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        decDelay = traceBack*3;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
    case 7%16-QAM 3/4
        traceBack = 3;
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        decDelay = traceBack*3;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
    case 8%64-QAM 3/4
        traceBack = 3;
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        decDelay = traceBack*3;
        dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
        decode_bit = dataOut(decDelay+1:end);
end





% traceBack = 1;
% tPoly = poly2trellis(7,[133,171]);
% dataOut = vitdec(codeword,tPoly,traceBack,'cont','hard');
% decDelay = traceBack;
% decode_bit = dataOut(decDelay+1:end);

% [numErrors,ber] = biterr(dataIn(1:end-decDelay),dataOut(decDelay+1:end));
% coded_bits = dataOut;