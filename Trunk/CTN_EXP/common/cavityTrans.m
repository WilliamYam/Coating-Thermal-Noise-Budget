function trans = cavityTrans(R1,R2,Loss,phi)

%cavityTrans:Cavity Transmission coefficient
%
%trans = cavityTrans(R1,R2,Loss,phi)

r1 = sqrt(R1);
r2 = sqrt(R2);
t1 = sqrt(1-R1);
t2 = sqrt(1-R2);

%Cavity transmission response
trans = (-t1*t2*(1-Loss)^(1/4)*exp(-1i*phi))./(1-r1*r2*sqrt(1-Loss)*exp(-1i*phi));




