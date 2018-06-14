clear all;
clc;
tic
%% 18.1.28
% the simulation function
%% following file path should be the corresponding path in your PC
% file_name = '...\np_rx_ofdm_entries3.2.hdf5';
csi_file_name = {
            'D:\simulation_ofdm\csi_seg\ls_csi_re6.txt';...
            'D:\simulation_ofdm\csi_seg\ls_csi_im6.txt';...
            }; 
esnr_file_name = {
                    'D:\simulation_ofdm\esnr\ls_esnr1.txt';...
                    'D:\simulation_ofdm\esnr\ls_esnr2.txt';...
                    'D:\simulation_ofdm\esnr\ls_esnr3.txt';...
                    'D:\simulation_ofdm\esnr\ls_esnr4.txt';...
                    'D:\simulation_ofdm\esnr\ls_esnr5.txt';...
                    'D:\simulation_ofdm\esnr\ls_esnr6.txt';
                  };
save_file = {
             'D:\simulation_ofdm\out_put\correct_poss_ls6.txt';...
             'D:\simulation_ofdm\out_put\choose_seq_ls6.txt';
             };
csire = importdata(csi_file_name{1,1});
csiim = importdata(csi_file_name{2,1});
csi_data = csire + 1j.*csiim;
esnr_data_all = importdata(esnr_file_name{6,1});
for i = 1:length(csi_data)
% for i = 1:2
    channel_choice = i;
    esnr_data = esnr_data_all(channel_choice,:);
    a = cal_through_bits(esnr_data,csi_data(:,channel_choice));
    save_txt(a,save_file{1,1});
end
t = toc
% save_txt(Choose_Seq(save_file{1,1})',save_file{2,1});






