function BinData = Dec2BinVector(DecData, BinLen)

BinData = zeros(BinLen, size(DecData, 1));
for k = 1: BinLen
    BinData(BinLen + 1 - k, :) = mod(DecData, 2);
    DecData = floor(DecData / 2);
end
BinData = reshape(BinData, [], 1);
