function esnr = get_esnr(file_name,channel_choice) 
%% calculate the esnr data corresponding the channel
ofdm_data = h5read(file_name,'/RX_OFDM');
ofdm_data.power = ofdm_data.power-min(ofdm_data.power)+10;
double(ofdm_data.power);
% figure(1);
% plot(ofdm_data.power);
esnr = ofdm_data.power(channel_choice,1);


