function [ErrorSig, Ptrans] = PDH_TRANS(R1,R2,Loss,Pinc,phiMod,deltaMod,demodPhase,phi)

%PDH_TRANS: calculates error signal in transmission and transmitted power
%
%[ErrorSig, Ptrans] = PDH_TRANS(R1,R2,Loss,Pinc,phiMod,deltaMod,demodPhase,phi)

%Vectorise the carrier and sideband components
phi_v = [phi; phi+phiMod; phi-phiMod];

%Transmission coefficient
fTrans = cavityTrans(R1,R2,Loss,phi_v);

%Power in the Carrier and the Sidebands
Pc = Pinc*(besselj(0,deltaMod))^2;
Ps = Pinc*(besselj(1,deltaMod))^2;

%Error signal
Err = 2*sqrt(Pc*Ps)*(fTrans(1,:).*conj(fTrans(2,:))-conj(fTrans(1,:)).*fTrans(3,:));
ErrorSig = real(Err)*cos(demodPhase)+imag(Err)*sin(demodPhase);

%Transmitted power
Ptrans = Pc*(abs(cavityTrans(R1,R2,Loss,0))).^2 + Ps*(abs(cavityTrans(R1,R2,Loss,phiMod))).^2 +...
    Ps*(abs(cavityTrans(R1,R2,Loss,-phiMod))).^2;



