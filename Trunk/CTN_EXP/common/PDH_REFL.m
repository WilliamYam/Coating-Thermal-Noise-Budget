function out = PDH_REFL(R1,R2,Loss,Pinc,phiMod,deltaMod,demodPhase,phi)

%PDH_TRANS: calculates error signal in transmission and transmitted power
%
%[ErrorSig, Ptrans] = PDH_TRANS(R1,R2,Loss,Pinc,phiMod,deltaMod,demodPhase,phi)

%Vectorise the carrier and sideband components
phi_v = [phi; phi+phiMod; phi-phiMod];
% 
% 
% %Cavity reflectivity response
% Frefl = (r1-r2*sqrt(1-Loss)*exp(-1i*phi_v))./(1-r1*r2*sqrt(1-Loss)*exp(-1i*phi_v));

%Transmission coefficient
Frefl = cavityRefl(R1,R2,Loss,phi_v);

%Power in the Carrier and the Sidebands
Pc = Pinc*(besselj(0,deltaMod))^2;
Ps = Pinc*(besselj(1,deltaMod))^2;

Err = 2*sqrt(Pc*Ps)*(Frefl(1,:).*conj(Frefl(2,:))-conj(Frefl(1,:)).*Frefl(3,:));

out = real(Err)*cos(demodPhase)+imag(Err)*sin(demodPhase);
