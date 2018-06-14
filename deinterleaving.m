function x=deinterleaving(y,Nbpsc)
%% �⽻֯

Ncbps = length(y);
x=y; % To make x have the same size as y
i=0:Ncbps-1; s= max(Nbpsc/2,1); 
j=s*floor(i/s) + mod(i+Ncbps-floor(16*i/Ncbps),s); % 1st permutation
z(j+1)=y;
k=i;
i = Ncbps/16*mod(k,16) + floor(k/16); % 2nd permutation
x(i+1)=z;
