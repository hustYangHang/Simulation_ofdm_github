function csi_ifft = Pro_nphdf5(file_name,choice)
%% 取出对应信道的csi数据
ofdm_data = h5read(file_name,'/RX_OFDM');
ofdm_csi = ofdm_data.chan_est;

%%  plot_csi
ofdm_csi  = double((ofdm_csi));
[m,n,k] = size(ofdm_csi);
csi_data = zeros(n,k);

for i = 1:k
    csi_data(:,i) = (ofdm_csi(1,:,i)/2^15) + 1j*(ofdm_csi(2,:,i)/2^15);
end

% cis_ifft = zeros(n,k);
switch choice
    case 1
        csi_ifft = ifft(csi_data,64,1);%time region
    case 2
        csi_ifft = csi_data;%frequence region
end




