function choose_seq = Choose_Seq(File_name)
correct_poss = importdata(File_name);
[m,n] = size(correct_poss);
choose_seq = zeros(m,1);
for i = 1:m
    for j = 1:n
        if correct_poss(i,j) <= 0.01
            choose_seq(i,1) = j-1;
        end
    end
end