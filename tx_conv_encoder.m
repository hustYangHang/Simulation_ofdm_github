function coded_bits = tx_conv_encoder(MOD_ORDER,in_bits)
% % ¾í»ý±àÂë
dataIn = in_bits;
switch MOD_ORDER
    case 1%BPSK 1/2
        traceBack = 1;
        dataIn = [dataIn;zeros(1,traceBack)'];
        tPoly = poly2trellis(7,[133,171]);
        codeword = convenc(dataIn,tPoly);        
    case 2%QPSK 1/2
        traceBack = 1;
        dataIn = [dataIn;zeros(1,traceBack)'];
        tPoly = poly2trellis(7,[133,171]);
        codeword = convenc(dataIn,tPoly);      
    case 3%16-QAM 1/2
        traceBack = 1;
        dataIn = [dataIn;zeros(1,traceBack)'];
        tPoly = poly2trellis(7,[133,171]);
        codeword = convenc(dataIn,tPoly);       
    case 4%64-QAM 2/3
        traceBack = 2;
        dataIn = [dataIn;zeros(1,2*traceBack)'];
        tPoly = poly2trellis([5 4],[23 35 0; 0 5 13]);
        codeword = convenc(dataIn,tPoly);
    case 5%BPSK 3/4
        traceBack = 3;
        dataIn = [dataIn;zeros(1,2*traceBack)'];
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        codeword = convenc(dataIn,tPoly);
    case 6%QPSK 3/4
        traceBack = 3;
        dataIn = [dataIn;zeros(1,2*traceBack)'];
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        codeword = convenc(dataIn,tPoly);
    case 7%16-QAM 3/4
        traceBack = 3;
        dataIn = [dataIn;zeros(1,2*traceBack)'];
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        codeword = convenc(dataIn,tPoly);
    case 8%64-QAM 3/4
        traceBack = 3;
        dataIn = [dataIn;zeros(1,2*traceBack)'];
        tPoly = poly2trellis([5 5 5],[23 35 11 0; 0 5 13 21;27 31 0 15]);
        codeword = convenc(dataIn,tPoly);        
end

coded_bits = codeword;


