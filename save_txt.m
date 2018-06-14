function save_txt(Array,file_path)
[m,n]=size(Array);
fid = fopen(file_path,'a');
for i=1:n
    for j=1:m
        fprintf(fid,'%f',Array(j,i));
        fprintf(fid,'\t');
    end
    fprintf(fid,'\r\n');
end
fclose(fid);