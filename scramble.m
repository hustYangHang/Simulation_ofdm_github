function scrData = scramble(input_data,choice)
%% 18.1.24

scrambler = comm.Scrambler(2,[1 0 1 1 0 1 0 1],[0 1 0 0 1 0 1]); 
descrambler = comm.Descrambler(2,[1 0 1 1 0 1 0 1],[0 1 0 0 1 0 1]);
switch choice
    case 1
        scrData = scrambler(input_data);%scramble the data 
    case 2
        scrData = descrambler(input_data);%descramble the data
end
