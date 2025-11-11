function Error = ComputeStreeError_A(ck,E11,S11,E22,S22)


lambda11=sqrt(1+2*E11);
lambda22=sqrt(1+2*E22);

I4=lambda11.^2;
I4p=lambda22.^2; % perpendicular

E_par=(1-ck(4))*(I4-1); % parallele
E_per=ck(4)*(I4p-1); % perpendicular

W1=ck(1)/2;
W4=ck(2)*(1-ck(4))*(E_par).*exp(ck(3)*(E_par.^2+E_per.^2));
W4p=ck(2)*(ck(4))*(E_per).*exp(ck(3)*(E_par.^2+E_per.^2));

S11_th=ck(1)*(1-1./(lambda11.^4.*lambda22.^2))+2*W4;
S22_th=ck(1)*(1-1./(lambda11.^2.*lambda22.^4))+2*W4p;

ER11=abs(S11_th-S11);
ER22=abs(S22_th-S22);

Error=[ER11; ER22];  




end