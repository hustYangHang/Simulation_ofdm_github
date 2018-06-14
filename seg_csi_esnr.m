clear all;
clc;
%% 将csi分段
nphdf5_file_name = {
                    'D:\uav_data\hdf5data_180502\np_rx_ofdm_entries_01.hdf5';...
            }; 
csi_fre_data_all = Pro_nphdf5(nphdf5_file_name{1,1},2);
seg_num = 6;
csi_fre_data_all_new = csi_fre_data_all(:,1:seg_num*floor(length(csi_fre_data_all)/seg_num));
for i = 1:seg_num
    save_txt(csi_fre_data_all_new(:,floor(length(csi_fre_data_all)/seg_num)*(i-1)+1:...
                floor(length(csi_fre_data_all)/seg_num)*i)',...
                ['D:\simulation_ofdm\csi_seg\hs_csi_re' num2str(i) '.txt']);
    save_txt((csi_fre_data_all_new(:,floor(length(csi_fre_data_all)/seg_num)*(i-1)+1:...
                floor(length(csi_fre_data_all)/seg_num)*i)*(-1j))',...
                ['D:\simulation_ofdm\csi_seg\hs_csi_im' num2str(i) '.txt']);
end


%% 将esnr分段
% esnr_file_name = {
%                   'D:\simulation_ofdm\esnr\esnr_ls.txt';...
%                   'D:\simulation_ofdm\esnr\esnr_ms.txt';...
%                   'D:\simulation_ofdm\esnr\esnr_hs.txt';
%                   };
% esnr_data_all = importdata(esnr_file_name{3,1});
% seg_num = 3;
% esnr_data_all_new = esnr_data_all(1:seg_num*floor(length(esnr_data_all)/seg_num),:);
% for i = 1:seg_num
%     save_txt(esnr_data_all_new(floor(length(esnr_data_all)/seg_num)*(i-1)+1:...
%                 floor(length(esnr_data_all)/seg_num)*i,:)',...
%                 ['D:\simulation_ofdm\esnr\hs_esnr' num2str(i) '.txt']);
% end