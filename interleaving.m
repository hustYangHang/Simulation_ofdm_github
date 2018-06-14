function y=interleaving(x,Nbpsc)
%% ½»Ö¯

x = [x;zeros(16*ceil(length(x)/16)-length(x),1)];
y=x; % To make y have the same size as x
Ncbps = length(x);
k=0:Ncbps-1;
i = Ncbps/16*mod(k,16) + floor(k/16); % 1st permutation
z=x(i+1);
i=k;
s= max(Nbpsc/2,1); 
j=s*floor(i/s) + mod(i+Ncbps-floor(16*i/Ncbps),s); % 2nd permutation
y=z(j+1);