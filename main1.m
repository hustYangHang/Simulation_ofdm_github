clear all;
clc;
%% 18.1.28
% the simulation function
%% following file path should be the corresponding path in your PC
% file_name = '...\np_rx_ofdm_entries3.2.hdf5';
csi_file_name = {
            'D:\Simple_Use\csi_real_data_002.txt';...
            'D:\Simple_Use\csi_imag_data_002.txt';...
            }; 
esnr_file_name = {
                    'D:\Simple_Use\snr_002.txt';...
                  };
save_file = {
             'D:\Simple_Use\output\correct_pass_test002.txt';...
%              'D:\simulation_ofdm\out_put\choose_seq_001.txt';
             };
csire = importdata(csi_file_name{1,1})';
csiim = importdata(csi_file_name{2,1})';
csi_data = csire + 1j.*csiim;
esnr_data_all = importdata(esnr_file_name{1,1});
data_num = length(esnr_data_all);

ber = zeros(100,8);
h = waitbar(0,'0','Name','Run schedule...');
for i = floor(data_num/5):floor(data_num/5) + 99
    channel_choice = i;
    esnr_data = esnr_data_all(channel_choice,1);
%     esnr_data = 50;
    ber(i+1-floor(data_num/5),:) = cal_through_bits(esnr_data,csi_data(:,channel_choice))';
    waitbar((i+1-floor(data_num/5))/(100),h,sprintf('%12.4f',(i+1-floor(data_num/5))/(100)));
%     save_txt(ber,save_file{1,1});
end
save_txt(ber','D:\simulation_ofdm\ber_csi.txt');
% save_txt(Choose_Seq(save_file{1,1})',save_file{2,1});
%% determine esnr model figure data
% h = waitbar(0,'0','Name','Run schedule...');
% temp = 0;
% snr_range = 50;
% repeat_num = 100;
% for i = 1:snr_range
%     channel_choice = floor(data_num/5);
%     temp = (i-1)*repeat_num;
%     for k = 1:repeat_num
%         esnr_data = i;
%         ber = cal_through_bits(esnr_data,csi_data(:,channel_choice));
%         save_txt(ber,save_file{1,1});
%         waitbar((temp+k)/(snr_range*repeat_num),h,sprintf('%12.4f',(temp+k)/(snr_range*repeat_num)));
%     end
% end
% close(h);




