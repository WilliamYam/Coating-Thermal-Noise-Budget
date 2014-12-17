
function ctnParams

%ctnParms  Defines parameters for the CTN Simulink model

% It outputs to a variable called "fcsParams" in the matlab workspace.


%% Common Parameters
fcsParams.paramsFileName = mfilename('fullpath');
%fcsParams.name = 'ctnFCS';

fcsParams.freq = logspace(log10(2), log10(1e5), 1000);   

%%  Data loading
load('ctnData.mat') 

%% Common Parameters
% 
fcsParams.common.c = 299792458;                                             % [m / s] speed of light  
fcsParams.common.wavelength = 1064e-9;                                      % [m] laser wavelength 
% fcsParams.common.e_ch          = 1.602e-19;     % charge per electron [C/electron]
fcsParams.common.h = 6.626e-34;                                             % [J / s] Planck's constant

%% PD and Mixer parameters

%%% Incident light PDA100A detector (incident onto cavity)
PDtranimpInc = 1135;                                                        % [V / W]

%%% 1811
PDtranimp00 = 0.8 * 4e4;                                                    % [A / W * Vrf / A]   PDgain for 00

% %%% LSC RF PD 
% PDtranimp02 = 4.0e4;                 % [V/W] transimpedance at 17.23MHz 
% PDtranimp20 = 2.0e4;                 % [V/W] transimpedance at 27.61MHz 

%%% QPD Trans 
PDtranimp02 = 6.882e3;                                                      % [V / W] (measured) 
PDtranimp20 = 6.882e3;                                                      % [V / W] (measured)

%%%Mixer gain
fcsParams.PD.gMixer00 = 1;
fcsParams.PD.gMixer02 = 0.5 * 100;                                          % [Vdc / Vref] mixer gain + SR560 gain
fcsParams.PD.gMixer20 = 0.5 * 100;                                          % [Vdc / Vref] mixer gain + SR560 gain

%transimpedance after mixer
fcsParams.PD.PDrespTot00 = PDtranimp00 * fcsParams.PD.gMixer00;             % [V / W]
fcsParams.PD.PDrespTot02 = PDtranimp02 * fcsParams.PD.gMixer02;             % [V / W]
fcsParams.PD.PDrespTot20 = PDtranimp20 * fcsParams.PD.gMixer20;             % [V / W]

%% Cavity Parameters
% 
fcsParams.cavity.Length = 0.095;                                            % [m] cavity length 
fcsParams.cavity.RoC1 = 0.050;                                              % [m] input coupler radius of curvature 
fcsParams.cavity.RoC2 = 0.050;                                              % [m] output coupler radius of curvature
fcsParams.cavity.RoC102 = (1 + 1.7e-4) * fcsParams.cavity.RoC1;             % [m] input coupler radius of curvature for 02
fcsParams.cavity.RoC202 = (1 + 1.7e-4) * fcsParams.cavity.RoC2;             % [m] output coupler radius of curvature for 02
fcsParams.cavity.RoC120 = (1 - 1.7e-4) * fcsParams.cavity.RoC1;             % [m] input coupler radius of curvature for 20
fcsParams.cavity.RoC220 = (1 - 1.7e-4) * fcsParams.cavity.RoC2;             % [m] output coupler radius of curvature for 20
fcsParams.cavity.R1 = 1 - 200e-6;                                           % Reflectivity of input coupler
fcsParams.cavity.R2 = 1 - 200e-6;                                           % Reflectivity of end coupler
fcsParams.cavity.FSR = fcsParams.common.c /...
                                    (2 * fcsParams.cavity.Length);          % FSR [Hz]
% fcsParams.cavity.LightTransitTime = fcsParams.cavity.Length /...
%                                     fcsParams.common.c ;                       % transit time per bounce [s]
%%%% Finesse, FWHM, Pole
            %-------- Finesse function def.-------------------------
                Fin = @(R1, R2)((pi * (R1 * R2)^(1 / 4)) / (1 - sqrt(R1 * R2)));
            %-------------------------------------------------------
            
fcsParams.cavity.Finesse = Fin(fcsParams.cavity.R1, fcsParams.cavity.R2);   % expected cavity finesse
%   
fcsParams.cavity.FWHM02 = (225.355e6 - 225.256e6) * 4;                      % [Hz] measured value of FWHM
fcsParams.cavity.FWHM20 = (226.455e6 - 226.385e6) * 4;                      % [Hz] measured value of FWHM
fcsParams.cavity.Finesse02 = fcsParams.cavity.FSR / fcsParams.cavity.FWHM02;   
fcsParams.cavity.Finesse20 = fcsParams.cavity.FSR / fcsParams.cavity.FWHM20;   
%%%%%%%
fcsParams.cavity.Loss02    = 2 * pi / fcsParams.cavity.Finesse02 +...       % cavity loss per round trip         
                         log(fcsParams.cavity.R1 * fcsParams.cavity.R2);
%%%%%%%
fcsParams.cavity.Loss20    = 2 * pi / fcsParams.cavity.Finesse20 +...       % cavity loss per round trip         
                         log(fcsParams.cavity.R1 * fcsParams.cavity.R2);                       
% fcsParams.cavity.Pole    = fcsParams.cavity.FWHM02/2;                             % [Hz] cavity pole


%%%% Reflection and Transmission coefficients
fcsParams.cavity.fRefl02 = cavityRefl(fcsParams.cavity.R1, fcsParams.cavity.R2,...
    fcsParams.cavity.Loss02, 0);                                            % Amplitude reflection coefficient at resonance
fcsParams.cavity.fRefl20 = cavityRefl(fcsParams.cavity.R1, fcsParams.cavity.R2,...
    fcsParams.cavity.Loss20, 0); 
fcsParams.cavity.Refl02 = abs(fcsParams.cavity.fRefl02)^2;
fcsParams.cavity.Refl20 = abs(fcsParams.cavity.fRefl20)^2;

fcsParams.cavity.fTrans02 = cavityTrans(fcsParams.cavity.R1, fcsParams.cavity.R2,...
    fcsParams.cavity.Loss02, 0);                                            % Amplitude transmission coefficient at resonance
fcsParams.cavity.fTrans20 = cavityTrans(fcsParams.cavity.R1, fcsParams.cavity.R2,...
    fcsParams.cavity.Loss20, 0); 
fcsParams.cavity.Trans02 = abs(fcsParams.cavity.fTrans02)^2;
fcsParams.cavity.Trans20 = abs(fcsParams.cavity.fTrans20)^2;

%%%% Gouy Phase

  %-------- Oneway Gouy phase def.----------------------------------------
   Gouy = @(n, m, L, RoC1, RoC2)((n + m + 1) * acos(-sqrt((1 - L / RoC1) * (1 - L / RoC1))));
  %-----------------------------------------------------------------------

fcsParams.cavity.Gouy00 = 2 * Gouy(0, 0, fcsParams.cavity.Length,...
                                     fcsParams.cavity.RoC1,...
                                     fcsParams.cavity.RoC2);                % [rad]  Round trip Gouy phase for 00
fcsParams.cavity.Gouy02 = 2 * Gouy(0, 2, fcsParams.cavity.Length,...
                                     fcsParams.cavity.RoC102,...
                                     fcsParams.cavity.RoC202);                % [rad]  Round trip Gouy phase for 02                                 
fcsParams.cavity.Gouy20 = 2 * Gouy(2, 0, fcsParams.cavity.Length,...
                                     fcsParams.cavity.RoC120,...
                                     fcsParams.cavity.RoC220);                % [rad]  Round trip Gouy phase for 20


%% Cavity beam 

fcsParams.beam.laserFreq00 = fcsParams.common.c /...
                                 fcsParams.common.wavelength;               % laser light freq.
fcsParams.beam.laserFreq02 = fcsParams.beam.laserFreq00 + 226e6 * 2;
fcsParams.beam.laserFreq20 = fcsParams.beam.laserFreq00 + 225e6 * 2;

powerFactor = 1;                                                            % Ideally scales budget with current HOM power
fcsParams.beam.PwrInc00 = 0.355e-3;                                         % [W] incident Gaussian on input coupler         
fcsParams.beam.PwrInc02 = (0.910 - 0.0242) / PDtranimpInc * powerFactor;    % [W] incident HOM02 on input coupler
fcsParams.beam.PwrInc20 = (1.19 - 0.0242) / PDtranimpInc * powerFactor;     % [W] incident HOM20 on input coupler
fcsParams.beam.Coupling00 = 0.20;                                           % cavity coupling of 00 
fcsParams.beam.Coupling02 = 0.182 / 4.12;                                   % cavity coupling of 02
fcsParams.beam.Coupling20 = 0.28 / 4.72;                                    % cavity coupling of 20


PwrCpl00 = fcsParams.beam.PwrInc00 * fcsParams.beam.Coupling00;             % [W]  power coupled to cavity 00
PwrCpl02 = fcsParams.beam.PwrInc02 * fcsParams.beam.Coupling02 / ...
    (1 - fcsParams.cavity.Refl02);                                          % [W]  power coupled to cavity 02
PwrCpl20 = fcsParams.beam.PwrInc20 * fcsParams.beam.Coupling20 / ...
    (1 - fcsParams.cavity.Refl20);                                          % [W]  power coupled to cavity 20

% fcsParams.beam.PwrCir00 = fcsParams.beam.PwrInc00 * ...
%                           fcsParams.beam.Coupling00 *...
%                           fcsParams.cavity.Finesse/pi;                  % [W] circulating power of 00 in the cavity
%                       
% fcsParams.beam.PwrCir02 = fcsParams.beam.PwrInc02 * ...
%                           fcsParams.beam.Coupling02 *...
%                           fcsParams.cavity.Finesse/pi;                  % [W] circulating power of 02 in the cavity 
%                       
% fcsParams.beam.PwrCir20 = fcsParams.beam.PwrInc20 * ...
%                           fcsParams.beam.Coupling20 *...
%                           fcsParams.cavity.Finesse/pi;                  % [W] circulating power of 20 in the cavity

%% Opt Gain
                     
%%%laser controller gain (fast channel)                      
fcsParams.OptGain.laserFast = 5e6;                                          % [Hz / V]

%%%Marconi gain
fcsParams.OptGain.PhaseDetectorGain02 = 4e-7;                               % [V / Hz] (for the Marconi noise measurement)
fcsParams.OptGain.PhaseDetectorGain20 = 4e-7;  
fcsParams.OptGain.MarconiGain02 = 10e3 / sqrt(2);                           % [Hz / V] note: Marconi external modulation per Vrms
fcsParams.OptGain.MarconiGain20 = 10e3 / sqrt(2);                           % [Hz / V] note: Marconi external modulation per Vrms
fcsParams.OptGain.PhaseDetectorGainBN = (0.233 + 0.256) / (2e5);            % [V / Hz]
fcsParams.OptGain.PhaseDetectorGainBN_refl = 9.4e-7;                        % [V / Hz] From 11/6/2014

%%%Frequency 2 roundtripPhase
fcsParams.OptGain.Hz2rad  = 4 * pi * fcsParams.cavity.Length /...
                                   fcsParams.common.c;                      % [rad / Hz] round trip phase (longitudonal wave)                                   fcsParams.common.c;  
%%%Length 2 roundtripPhase
fcsParams.OptGain.m2rad    =  4 * pi / fcsParams.common.wavelength;         % [rad / m]  Length2phase conversion

%%%Frequency 2 Length
fcsParams.OptGain.Hz2m00 = fcsParams.cavity.Length / fcsParams.beam.laserFreq00;    % [m / Hz]
fcsParams.OptGain.Hz2m02 = fcsParams.cavity.Length / fcsParams.beam.laserFreq02;    % [m / Hz]
fcsParams.OptGain.Hz2m20 = fcsParams.cavity.Length / fcsParams.beam.laserFreq20;    % [m / Hz]


%%%% GouyPhase contribution (to be completed) 
fcsParams.OptGain.GouyContr02 = fcsParams.cavity.Gouy02 /...                     
                                   fcsParams.OptGain.Hz2rad;                % [Hz] Gouy phase contribution to 02   
fcsParams.OptGain.GouyContr20 = fcsParams.cavity.Gouy20 /...
    fcsParams.OptGain.Hz2rad;                                               % [Hz] Gouy phase contribution to 20
fcsParams.OptGain.beatFreq = abs(fcsParams.OptGain.GouyContr02 -...
    fcsParams.OptGain.GouyContr20);

%% PDH on REFLECTION

%%% Phase Modulation Properties
fcsParams.errSig.beta00 = 0.46;                                             % Modulation depth for 00
% fcsParams.errSig.beta02 = 0.35;                  % Modulation depth for 02
% fcsParams.errSig.beta20 = 0.35;                  % Modulation depth for 20
%%% Modulation Frequency
fcsParams.errSig.modFreq00 = 29.00e6;                                       % Modulation frequency for 00
% fcsParams.errSig.modFreq02 = 27.61e6;            % Modulation frequency for 02
% fcsParams.errSig.modFreq20 = 17.23e6;            % Modulation frequency for 20
%%% Round trip phase for sidebands
phi_mod00 = 2 * pi * 2 * fcsParams.cavity.Length * ...
    fcsParams.errSig.modFreq00 / fcsParams.common.c;   %------------------check this!!!!
% phi_mod02 = 2*pi*2*fcsParams.cavity.Length*fcsParams.errSig.modFreq02/fcsParams.common.c;   % made as roundtrip phase
% phi_mod20 = 2*pi*2*fcsParams.cavity.Length*fcsParams.errSig.modFreq20/fcsParams.common.c;


% %%% Cavity Roundtrip Phase [rad] (vector for plot)
phi_cav = linspace(-pi / 1000, pi / 1000, 10001);

%%% The Error Signal [W] (demod phase = 90deg, sidebands outside linewidth, want I)
Err00 = PDH_REFL(fcsParams.cavity.R1, fcsParams.cavity.R2,...
                 fcsParams.cavity.Loss02, PwrCpl00, phi_mod00,...
                 fcsParams.errSig.beta00, pi / 2, phi_cav);                 % PDH error signal for 00 used for plot
% Err02 = PDH_Esig(fcsParams.cavity.R1,fcsParams.cavity.R2,...
%                  fcsParams.cavity.Loss,PwrCpl02,phi_mod02,...
%                  fcsParams.errSig.beta02,pi/2,phi_cav);                      %PDH error signal for 02 used for plot
% Err20 = PDH_Esig(fcsParams.cavity.R1,fcsParams.cavity.R2,...
%                  fcsParams.cavity.Loss,PwrCpl20,phi_mod20,...
%                  fcsParams.errSig.beta20,pi/2,phi_cav);                      %PDH error signal for 20 used for plot
% 
% errSig00_pk2pk = 2*max(Err00);        % [W] ppk-pk error signal for 00
% errSig02_pk2pk = 2*max(Err02);        % [W] ppk-pk error signal for 02
% errSig20_pk2pk = 2*max(Err20);        % [W] ppk-pk error signal for 20
% 
%  % [W / m] gradient of the Error Signal 
grad_Err00 = fcsParams.OptGain.m2rad * abs(gradient(Err00) ./ gradient(phi_cav));                        
% grad_Err02 = fcsParams.OptGain.m2rad*abs(gradient(Err02)./gradient(phi_cav));
% grad_Err20 = fcsParams.OptGain.m2rad*abs(gradient(Err20)./gradient(phi_cav));

%%% Cavity Response (Error Signal slope) 
fcsParams.errSig.m2W00e = max(grad_Err00);                                  % [W / m] (esimated)
% fcsParams.errSig.m2W02e = max(grad_Err02) * fcsParams.beam.Coupling00;    % [W/m]
% fcsParams.errSig.m2W20e = max(grad_Err20)*fcsParams.beam.Coupling00;      % [W/m]

fcsParams.errSig.m2W00m = 6.1913e+06 * fcsParams.beam.Coupling00;           % [W / m] (measured)
% fcsParams.errSig.m2W02m = 3.3599e+05 * powerFactor * fcsParams.beam.Coupling02;                   % [W/m]
% fcsParams.errSig.m2W20m = 2.4796e+05 * powerFactor * fcsParams.beam.Coupling20;                   % [W/m]

 

%% PDH on TRANSMISSION [W] (demod phase = 90deg, sidebands outside linewidth, want I)

%%% Phase Modulation Properties
fcsParams.errSig.beta02 = (78 / 250) * pi;                                  % Modulation depth for 02
fcsParams.errSig.beta20 = (83 / 250) * pi;                                  % Modulation depth for 20
%%% Modulation Frequency
fcsParams.errSig.modFreq02 = 3.2e4;                                         % Modulation frequency for 02
fcsParams.errSig.modFreq20 = 5.0e4;                                         % Modulation frequency for 20
%%% Round trip phase for sidebands
phi_mod02 = 2 * pi * 2 * fcsParams.cavity.Length *...
    fcsParams.errSig.modFreq02 / fcsParams.common.c;                        % made as roundtrip phase
phi_mod20 = 2 * pi * 2 * fcsParams.cavity.Length *...
    fcsParams.errSig.modFreq20 / fcsParams.common.c;

[Err02, Ptrans02] = PDH_TRANS(fcsParams.cavity.R1, fcsParams.cavity.R2,...
                    fcsParams.cavity.Loss02, PwrCpl02, phi_mod02,...
                    fcsParams.errSig.beta02, pi / 2, phi_cav);              % PDH error signal for 02 used for plot
[Err20, Ptrans20] = PDH_TRANS(fcsParams.cavity.R1, fcsParams.cavity.R2,...
                    fcsParams.cavity.Loss20, PwrCpl20, phi_mod20,...
                    fcsParams.errSig.beta20, pi / 2, phi_cav);              % PDH error signal for 20 used for plot


% errSig02_pk2pk = 2*max(Err02(:,1));                                         % [W] ppk-pk error signal for 02
% errSig20_pk2pk = 2*max(Err20(:,1));                                         % [W] ppk-pk error signal for 20


% [W / m] gradient of the Error Signal                    
grad_Err02 = fcsParams.OptGain.m2rad * abs(gradient(Err02) ./ gradient(phi_cav));   
grad_Err20 = fcsParams.OptGain.m2rad * abs(gradient(Err20) ./ gradient(phi_cav));   

%%% Cavity Response (Error Signal slope) 
fcsParams.errSig.m2W02e = max(grad_Err02);                                  % [W / m]  (estimated)
fcsParams.errSig.m2W20e = max(grad_Err20);                                  % [W / m]  (estimated)
fcsParams.errSig.m2W02m = (0.110 / 4e4) / fcsParams.OptGain.Hz2m02 / ...
    fcsParams.PD.PDrespTot02 * powerFactor;                                 % [W / m]  (measured)
fcsParams.errSig.m2W20m = (0.220 / 4e4) / fcsParams.OptGain.Hz2m20 / ...
    fcsParams.PD.PDrespTot20 * powerFactor;                                 % [W / m]  (measured)

fcsParams.errSig.Ptrans02e = Ptrans02 * powerFactor;                        % [W] Power on transmission (esimated)
fcsParams.errSig.Ptrans20e = Ptrans20 * powerFactor;                        % [W] Power on transmission (esimated)
fcsParams.errSig.Ptrans02m = 27.0e-3 / PDtranimp02 * powerFactor;           % [W] Power on transmission (measured)
fcsParams.errSig.Ptrans20m = 65.5e-3 / PDtranimp20 * powerFactor;           % [W] Power on transmission (measured)



%% Servo Parameters
   
%%%00
%when integrator and boost off- used for test
gain00 = -37.7;                                                             % [dB] measured gain at 
freq00 = 8e4;                                                               % [Hz] this trequency
fcsParams.servo00off.zeros = [1e6];                                         % [Hz]
fcsParams.servo00off.poles = [30];                                          % [Hz]
fcsParams.servo00off.gain = find_K(fcsParams.servo00off.zeros,...
                                       fcsParams.servo00off.poles,...
                                       gain00, freq00);         
gain00 = -37.7;                                                             % [dB] measured gain at 
freq00 = 8e4;                                                               % [Hz] this frequency
fcsParams.servo00.zeros = [10e3];                                           % [Hz]
fcsParams.servo00.poles = [0, 0];                                           % [Hz]
fcsParams.servo00.gain = find_K(fcsParams.servo00.zeros,...
                                       fcsParams.servo00.poles,...
                                       gain00, freq00);         
%%%HOM02
gain02 = 71 - mag2db(powerFactor);                                        % [dB] measured gain at 
freq02 = 30;                                                                % [Hz] this frequency
fcsParams.servo02.zeros = [];                                               % [Hz] Boost On
fcsParams.servo02.poles = [30, 11e3];                                       % [Hz]
fcsParams.servo02.gain = find_K(fcsParams.servo02.zeros,...
                                       fcsParams.servo02.poles,...
                                       gain02, freq02);         
%%%HOM20
gain20 = 66 - mag2db(powerFactor);                                          % [dB] measured gain at 
freq20 = 30;                                                                % [Hz] this frequency
fcsParams.servo20.zeros = [];                                               % [Hz] Boost On
fcsParams.servo20.poles = [30, 11e3];                                       % [Hz]
fcsParams.servo20.gain = find_K(fcsParams.servo20.zeros,...
                                       fcsParams.servo20.poles,...
                                       gain20, freq20);         

%% Transfer Functions (measured)

% servo TFs

%%%00 servo
fcsParams.TF.TFmagS00 = interp1(h00ServoMag(:, 1), h00ServoMag(:, 2), fcsParams.freq);
fcsParams.TF.TFphS00 = interp1(h00ServoPh(:,1), h00ServoPh(:, 2), fcsParams.freq);
%%%02 servo
fcsParams.TF.TFmagS02 = interp1(h02ServoMag(:, 1), h02ServoMag(:, 2), fcsParams.freq);
fcsParams.TF.TFphS02 = interp1(h02ServoPh(:,1), h02ServoPh(:, 2), fcsParams.freq);
%%%20 servo
fcsParams.TF.TFmagS20 = interp1(h20ServoMag(:, 1), h20ServoMag(:, 2), fcsParams.freq);
fcsParams.TF.TFphS20 = interp1(h20ServoPh(:, 1), h20ServoPh(:, 2), fcsParams.freq);

% marconi in to mixer out TFs

%%%02
fcsParams.TF.TFmagM2M02 = interp1(h02M2MMag(:, 1), h02M2MMag(:, 2), fcsParams.freq);
fcsParams.TF.TFphM2M02 = interp1(h02M2MPh(:,1), h02M2MPh(:, 2), fcsParams.freq);
%%%20
fcsParams.TF.TFmagM2M20 = interp1(h20M2MMag(:, 1), h20M2MMag(:, 2), fcsParams.freq);
fcsParams.TF.TFphM2M20 = interp1(h20M2MPh(:, 1), h20M2MPh(:, 2), fcsParams.freq);


% open loop TFs

%%%Gaussian
fcsParams.TF.TFmag00 = interp1(h00OLMag(:, 1), h00OLMag(:, 2), fcsParams.freq);
fcsParams.TF.TFph00 = interp1(h00OLPh(:, 1), h00OLPh(:, 2), fcsParams.freq);
%%%HOM02
fcsParams.TF.TFmag02 = interp1(h02OLMag(:, 1), h02OLMag(:, 2), fcsParams.freq);
fcsParams.TF.TFph02 = interp1(h02OLPh(:, 1), h02OLPh(:, 2), fcsParams.freq);
%%%HOM20
fcsParams.TF.TFmag20 = interp1(h20OLMag(:, 1), h20OLMag(:, 2), fcsParams.freq);
fcsParams.TF.TFph20 = interp1(h20OLPh(:, 1), h20OLPh(:, 2), fcsParams.freq);


assignin('base', 'fcsParams', fcsParams);

