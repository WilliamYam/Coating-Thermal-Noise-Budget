function refl = cavityRefl(R1,R2,Loss,phi)

%cavityTrans:Cavity Reflection coefficient
%
%trans = cavityRefl(R1,R2,Loss,phi)

r1 = sqrt(R1);
r2 = sqrt(R2);

%Cavity reflectivity response
refl = (r1-r2*sqrt(1-Loss)*exp(-1i*phi))./(1-r1*r2*sqrt(1-Loss)*exp(-1i*phi));




