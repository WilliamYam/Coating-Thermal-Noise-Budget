
function ctnParams

%ctnParms  Defines parameters for the HOM20.mdl Simulink model
% mode = 'data' or 'model'
% 'data'refers to measured quantity whereas 'model' to expected value of
% the same quantity.
% It outputs to a variable called "fclParams" in the matlab workspace.
% 
% based on the data by William, September 2014


%% Common Parameters
fcsParams.paramsFileName = mfilename('fullpath');
%fcsParams.name = 'ctnFCS';

fcsParams.freq = logspace(log10(2),log10(1e5),1000);   %%%---------------


%%  Data loading
load('ctnData.mat') 

%% Common Parameters
% 
fcsParams.common.c             = 299792458;     % speed of light  [m/s]
fcsParams.common.wavelength    = 1064e-9;        % laser wavelength [m]
% fcsParams.common.e_ch          = 1.602e-19;     % charge per electron [C/electron]
fcsParams.common.h             = 6.626e-34;     % Planck's constant [J/s]
fcsParams.common.laserFreq00   = fcsParams.common.c /...
                                 fcsParams.common.wavelength;   %laser light freq.
fcsParams.common.laserFreq02   = fcsParams.common.laserFreq00 + 226e6*2;
fcsParams.common.laserFreq20   = fcsParams.common.laserFreq00 + 225e6*2;


%% Cavity Parameters
% 
fcsParams.cavity.Length           = 0.095;                                    % [m] cavity length 
% fcsParams.cavity.RoC1             = 0.050;                                    % [m] input coupler radius of curvature 
% fcsParams.cavity.RoC2             = 0.050;                                    % [m] output coupler radius of curvature 
 fcsParams.cavity.R1               = 1-200e-6;                                 % Reflectivity of input coupler
 fcsParams.cavity.R2               = 1-200e-6;                                 % Reflectivity of end coupler
fcsParams.cavity.FSR              = fcsParams.common.c /...
                                    (2 * fcsParams.cavity.Length);              %FSR [Hz]
% fcsParams.cavity.LightTransitTime = fcsParams.cavity.Length /...
%                                     fcsParams.common.c ;                       % transit time per bounce [s]
% %%%% Finesse, FWHM, Pole
%             %-------- Finesse function def.-------------------------
%                 Fin = @(R1,R2)((pi*(R1*R2)^(1/4))/(1-sqrt(R1*R2)));
%             %-------------------------------------------------------
%             
% fcsParams.cavity.Finesse = Fin(fcsParams.cavity.R1,fcsParams.cavity.R2);       % expected cavity finesse
% 
%   
fcsParams.cavity.FWHM02    =  (225.389e6-225.284e6)*4;                                              % [Hz] measured value of FWHM
fcsParams.cavity.FWHM20    =  (226.496e6-226.424e6)*4;                                              % [Hz] measured value of FWHM
fcsParams.cavity.Finesse02 =  fcsParams.cavity.FSR/fcsParams.cavity.FWHM02;   
fcsParams.cavity.Finesse20 =  fcsParams.cavity.FSR/fcsParams.cavity.FWHM20;   
%%%%%%%
fcsParams.cavity.Loss02    = 2 * pi/ fcsParams.cavity.Finesse02 +...                % cavity loss per round trip         
                         log(fcsParams.cavity.R1*fcsParams.cavity.R2);
%%%%%%%
fcsParams.cavity.Loss20    = 2 * pi/ fcsParams.cavity.Finesse20 +...                % cavity loss per round trip         
                         log(fcsParams.cavity.R1*fcsParams.cavity.R2);                       
fcsParams.cavity.Pole    = fcsParams.cavity.FWHM02/2;                             % [Hz] cavity pole


% %%%% Gouy Phase
% 
%   %-------- Oneway Gouy phase def.----------------------------------------
%    Gouy = @(n,m,L,RoC1, RoC2)((n+m+1)*acos(-sqrt((1-L/RoC1)*(1-L/RoC1))));
%   %-----------------------------------------------------------------------
% 
%   fcsParams.cavity.Gouy00 = 2 * Gouy(0, 0, fcsParams.cavity.Length,...
%                                      fcsParams.cavity.RoC1,...
%                                      fcsParams.cavity.RoC1);                        % [rad]  Rounttrip Gout phase for 00
%   fcsParams.cavity.Gouy02 = 2 * Gouy(0, 2, fcsParams.cavity.Length,...
%                                      fcsParams.cavity.RoC1,...
%                                      fcsParams.cavity.RoC1);                        % [rad]  Rounttrip Gout phase for 02
%                                  
%    fcsParams.cavity.Gouy20 = 2 * Gouy(2, 0, fcsParams.cavity.Length,...
%                                      fcsParams.cavity.RoC1,...
%                                      fcsParams.cavity.RoC1);                        % [rad]  Rounttrip Gout phase for 20

%% Cavity beam 
%
powerFactor = 1;                                                % Ideally scales budget with current HOM power
fcsParams.beam.PwrInc00   =    0.355e-3;                        % [W] incident Gaussian on input coupler
%fcsParams.beam.PwrInc02   =    .508/1135 * powerFactor;          % [W] incident HOM02 on input coupler
%fcsParams.beam.PwrInc20   =    .824/1135 * powerFactor;          % [W] incident HOM20 on input coupler
fcsParams.beam.PwrInc02   =    .7e-3 * powerFactor;
fcsParams.beam.PwrInc20   =    .7e-3 * powerFactor;
fcsParams.beam.Coupling00   =   0.20;                           % cavity couppling of 00 
fcsParams.beam.Coupling02   =   1*0.182/4.12;                         % cavity couppling of 02
fcsParams.beam.Coupling20   =   1*0.28/4.72;                         % cavity couppling of 20


PwrCpl00 = fcsParams.beam.PwrInc00 * fcsParams.beam.Coupling00;    % [W]  power coupled to cavity 00
PwrCpl02 = fcsParams.beam.PwrInc02 * fcsParams.beam.Coupling02;    % [W]  power coupled to cavity 02
PwrCpl20 = fcsParams.beam.PwrInc20 * fcsParams.beam.Coupling20;    % [W]  power coupled to cavity 20

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
fcsParams.OptGain.laserFast      = 5e6;                                            % [Hz/V]

%%%Marconi gain
fcsParams.OptGain.PhaseDetectorGain02  =  4e-7;                                    % [V/Hz] ???????? (for the Marconi noise measurement)
fcsParams.OptGain.PhaseDetectorGain20  =  4e-7;  
fcsParams.OptGain.MarconiGain        = 10e3 / sqrt(2);                             % [Hz/V], note: Marconi external modulation per Vrms
fcsParams.OptGain.PhaseDetectorGainBN = 0.412/(2e5);                                    % [V/Hz]
fcsParams.OptGain.PhaseDetectorGainBN_refl = 9.4e-7;                                % [V/Hz] From 11/6/2014

% %%%Frequency 2 roundtripPhase
% fcsParams.OptGain.Hz2rad  = 4 * pi * fcsParams.cavity.Length /...
%                                    fcsParams.common.c;                             % [rad/Hz] roundtrip phase (longitudonal wave)                                   fcsParams.common.c;  
%%%Length 2 roundtripPhase
fcsParams.OptGain.m2rad    =  4 * pi / fcsParams.common.wavelength;                % [rad/m]  Length2phase conversion

%%%Frequency 2 Length
fcsParams.OptGain.Hz2m00 = fcsParams.cavity.Length / fcsParams.common.laserFreq00;  % [m/Hz]
fcsParams.OptGain.Hz2m02 = fcsParams.cavity.Length / fcsParams.common.laserFreq02;  % [m/Hz]
fcsParams.OptGain.Hz2m20 = fcsParams.cavity.Length / fcsParams.common.laserFreq20;  % [m/Hz]


%%%% GouyPhase contribution (to be completed) 
% fcsParams.OptGain.GouyContr02 = fcsParams.cavity.Gouy02 /...                     % [rad/Hz] Gouyphase contribution to 02
%                                    fcsParams.common.c;   
% fcsParams.OptGain.GouyContr20 = fcsParams.cavity.Gouy20 /...                     % [rad/Hz] Gouyphase contribution to 02
% 

%% PDH on REFLECTION

%%% Phase Modulation Properties
fcsParams.errSig.beta00 = 0.46;                  % Modulation depth for 00
% fcsParams.errSig.beta02 = 0.35;                  % Modulation depth for 02
% fcsParams.errSig.beta20 = 0.35;                  % Modulation depth for 20
%%% Modulation Frequency
fcsParams.errSig.modFreq00 = 29.00e6;            % Modulation frequency for 00
% fcsParams.errSig.modFreq02 = 27.61e6;            % Modulation frequency for 02
% fcsParams.errSig.modFreq20 = 17.23e6;            % Modulation frequency for 20
%%% Round trip phase for sidebands
phi_mod00 = 2*pi*2*fcsParams.cavity.Length*fcsParams.errSig.modFreq00/fcsParams.common.c;   %------------------check this!!!!
% phi_mod02 = 2*pi*2*fcsParams.cavity.Length*fcsParams.errSig.modFreq02/fcsParams.common.c;   % made as roundtrip phase
% phi_mod20 = 2*pi*2*fcsParams.cavity.Length*fcsParams.errSig.modFreq20/fcsParams.common.c;


% %%% Cavity Roundtrip Phase [rad] (vector for plot)
phi_cav = linspace(-pi/1000,pi/1000,10001);

%%% The Error Signal [W] (demod phase = 90deg, sidebands outside linewidth, want I)
Err00 = PDH_REFL(fcsParams.cavity.R1,fcsParams.cavity.R2,...
                 fcsParams.cavity.Loss02,PwrCpl00,phi_mod00,...
                 fcsParams.errSig.beta00,pi/2,phi_cav);                      %PDH error signal for 00 used for plot
% Err02 = PDH_Esig(fcsParams.cavity.R1,fcsParams.cavity.R2,...
%                  fcsParams.cavity.Loss,PwrCpl02,phi_mod02,...
%                  fcsParams.errSig.beta02,pi/2,phi_cav);                      %PDH error signal for 02 used for plot
% Err20 = PDH_Esig(fcsParams.cavity.R1,fcsParams.cavity.R2,...
%                  fcsParams.cavity.Loss,PwrCpl20,phi_mod20,...
%                  fcsParams.errSig.beta20,pi/2,phi_cav);                      %PDH error signal for 20 used for plot
% 
errSig00_pk2pk = 2*max(Err00);        % [W] ppk-pk error signal for 00
% errSig02_pk2pk = 2*max(Err02);        % [W] ppk-pk error signal for 02
% errSig20_pk2pk = 2*max(Err20);        % [W] ppk-pk error signal for 20
% 
%  % [W/m] gradient of the Error Signal 
grad_Err00 = fcsParams.OptGain.m2rad*abs(gradient(Err00)./gradient(phi_cav));                        
% grad_Err02 = fcsParams.OptGain.m2rad*abs(gradient(Err02)./gradient(phi_cav));
% grad_Err20 = fcsParams.OptGain.m2rad*abs(gradient(Err20)./gradient(phi_cav));

%%% Cavity Response (Error Signal slope) 
fcsParams.errSig.m2W00e = max(grad_Err00);                                % [W/m] (esimated)
% fcsParams.errSig.m2W02e = max(grad_Err02) * fcsParams.beam.Coupling00;    % [W/m]
% fcsParams.errSig.m2W20e = max(grad_Err20)*fcsParams.beam.Coupling00;      % [W/m]

fcsParams.errSig.m2W00m = 6.1913e+06 * fcsParams.beam.Coupling00;   % [W/m] (measured)
% fcsParams.errSig.m2W02m = 3.3599e+05 * powerFactor * fcsParams.beam.Coupling02;                   % [W/m]
% fcsParams.errSig.m2W20m = 2.4796e+05 * powerFactor * fcsParams.beam.Coupling20;                   % [W/m]

 

%% PDH on TRANSMISSION [W] (demod phase = 90deg, sidebands outside linewidth, want I)

%%% Phase Modulation Properties
fcsParams.errSig.beta02 = (75/250)*pi;                  % Modulation depth for 02
fcsParams.errSig.beta20 = (80/250)*pi;                  % Modulation depth for 20
%%% Modulation Frequency
fcsParams.errSig.modFreq02 = 5.7e4;            % Modulation frequency for 02
fcsParams.errSig.modFreq20 = 5.7e4;            % Modulation frequency for 20
%%% Round trip phase for sidebands
phi_mod02 = 2*pi*2*fcsParams.cavity.Length*fcsParams.errSig.modFreq02/fcsParams.common.c;   % made as roundtrip phase
phi_mod20 = 2*pi*2*fcsParams.cavity.Length*fcsParams.errSig.modFreq20/fcsParams.common.c;

[Err02, Ptrans02] = PDH_TRANS(fcsParams.cavity.R1,fcsParams.cavity.R2,...
                    fcsParams.cavity.Loss02,PwrCpl02,phi_mod02,...
                    fcsParams.errSig.beta02,pi/2,phi_cav);                  %PDH error signal for 02 used for plot
[Err20, Ptrans20] = PDH_TRANS(fcsParams.cavity.R1,fcsParams.cavity.R2,...
                    fcsParams.cavity.Loss20,PwrCpl20,phi_mod20,...
                    fcsParams.errSig.beta20,pi/2,phi_cav);                  %PDH error signal for 20 used for plot


% errSig02_pk2pk = 2*max(Err02(:,1));                                         % [W] ppk-pk error signal for 02
% errSig20_pk2pk = 2*max(Err20(:,1));                                         % [W] ppk-pk error signal for 20


 % [W/m] gradient of the Error Signal                    
grad_Err02 = fcsParams.OptGain.m2rad*abs(gradient(Err02)./gradient(phi_cav));   
grad_Err20 = fcsParams.OptGain.m2rad*abs(gradient(Err20)./gradient(phi_cav));   

%%% Cavity Response (Error Signal slope) 
fcsParams.errSig.m2W02e = max(grad_Err02);                                %[W/m]  (estimated)
fcsParams.errSig.m2W20e = max(grad_Err20);                                %[W/m]  (estimated)

fcsParams.errSig.m2W02m = (0.02/4e4)/fcsParams.OptGain.Hz2m02/(0.8 * 10.0e3)/50 * powerFactor;                   % [W/m]
fcsParams.errSig.m2W20m = (0.0472/4e4)/fcsParams.OptGain.Hz2m20/(0.8 * 10.0e3)/50 * powerFactor;                   % [W/m]

fcsParams.errSig.Ptrans02e = Ptrans02 * powerFactor;                                      %[W] Power on transmission (esimated)
fcsParams.errSig.Ptrans20e = Ptrans20 * powerFactor;                                      %[W] Power on transmission (esimated)

%fcsParams.errSig.Ptrans02m = 13.7e-3/6882 * powerFactor;                                        %[W] Power on transmission (measured)

%%%%%%%%%%%
fcsParams.errSig.Ptrans02m = 25.2e-3/6882 * powerFactor; 
fcsParams.errSig.Ptrans20m = 45.5e-3/6882 * powerFactor;                                       %[W] Power on transmission (measured)


%% PD and Mixer parameters

%%%Mixer gain
fcsParams.PD.G_Mixer00 = 1;
fcsParams.PD.G_Mixer02 = 50;            % mixer gain Vdc/Vref + SR560 gain of 50
fcsParams.PD.G_Mixer20 = 50;            % mixer gain Vdc/Vref + SR560 gain of 50
%%% 1811
PDtranimp00 = 0.8 * 4e4;             % [A / W * Vrf / A]   PDgain for 00
% %%% LSC RF PD 
% PDtranimp02 = 4.0e4;                 % [V/W] transimpedance at 17.23MHz 
% PDtranimp20 = 2.0e4;                 % [V/W] transimpedance at 27.61MHz 
%%% QPD Trans 
PDtranimp02 = 0.8 * 10.0e3;                 % [V/W]  
PDtranimp20 = 0.8 * 10.0e3;                 % [V/W] 

%transimpedance after mixer
fcsParams.PD.PDrespTot00 = PDtranimp00 * fcsParams.PD.G_Mixer00;    % [V/W]
fcsParams.PD.PDrespTot02 = PDtranimp02 * fcsParams.PD.G_Mixer02;    % [V/W]
fcsParams.PD.PDrespTot20 = PDtranimp20 * fcsParams.PD.G_Mixer20;    % [V/W]



%% Servo Parameters
   
%%%00
%when integrator and boost off- used for test
gain00        = -37.7;                       % dB measured gain at 
freq00        = 8e4;                          % Hz this trequency
fcsParams.servo00off.zeros       = [1e6];       % Hz
fcsParams.servo00off.poles       = [30];         % Hz
fcsParams.servo00off.gain        = find_K(fcsParams.servo00off.zeros,...
                                       fcsParams.servo00off.poles,...
                                       gain00, freq00);         % dB
gain00        = -37.7;                          % dB measured gain at 
freq00        = 8e4;                            % Hz this frequency
fcsParams.servo00.zeros       = [10e3];         % Hz
fcsParams.servo00.poles       = [0,0];          % Hz
fcsParams.servo00.gain        = find_K(fcsParams.servo00.zeros,...
                                       fcsParams.servo00.poles,...
                                       gain00, freq00);         % dB
%%%HOM02
gain02        = 93 - mag2db(powerFactor);  % dB measured gain at 
freq02        = 30;                              % Hz this frequency
fcsParams.servo02.zeros       = [];             % Hz
fcsParams.servo02.poles       = [8e5,8e5,1e6,0];    % Hz Boost On
fcsParams.servo02.gain        = find_K(fcsParams.servo02.zeros,...
                                       fcsParams.servo02.poles,...
                                       gain02, freq02);         % dB
%%%HOM20
gain20        = 94 - mag2db(powerFactor);    % dB measured gain at 
freq20        = 30;                              % Hz this frequency
fcsParams.servo20.zeros       = [];             % Hz Boost On
fcsParams.servo20.poles       = [8e5,8e5,1e6,0];    % Hz
fcsParams.servo20.gain        = find_K(fcsParams.servo20.zeros,...
                                       fcsParams.servo20.poles,...
                                       gain20, freq20);         % dB

%% Transfer Function (measured)

% servo TFs

%%%00 servo
fcsParams.TF.TFmagS00 = interp1(h00ServoMag(:,1),h00ServoMag(:,2),fcsParams.freq);
fcsParams.TF.TFphS00 = interp1(h00ServoPh(:,1),h00ServoPh(:,2),fcsParams.freq);
%%%02 servo
fcsParams.TF.TFmagS02 = interp1(h02ServoMag(:,1),h02ServoMag(:,2),fcsParams.freq);
fcsParams.TF.TFphS02 = interp1(h02ServoPh(:,1),h02ServoPh(:,2),fcsParams.freq);
%%%00 servo
fcsParams.TF.TFmagS20 = interp1(h20ServoMag(:,1),h20ServoMag(:,2),fcsParams.freq);
fcsParams.TF.TFphS20 = interp1(h20ServoPh(:,1),h20ServoPh(:,2),fcsParams.freq);


% measured at the closed loop

%%%Gaussian
fcsParams.TF.TFmag00 = interp1(h00OLMag(:,1),h00OLMag(:,2),fcsParams.freq);
fcsParams.TF.TFph00 = interp1(h00OLPh(:,1),h00OLPh(:,2),fcsParams.freq);
%%%HOM02
fcsParams.TF.TFmag02 = interp1(h02OLMag(:,1),h02OLMag(:,2),fcsParams.freq);
fcsParams.TF.TFph02 = interp1(h02OLPh(:,1),h02OLPh(:,2),fcsParams.freq);
%%%HOM20
fcsParams.TF.TFmag20 = interp1(h20OLMag(:,1),h20OLMag(:,2),fcsParams.freq);
fcsParams.TF.TFph20 = interp1(h20OLPh(:,1),h20OLPh(:,2),fcsParams.freq);
assignin('base', 'fcsParams', fcsParams);


fcsParams.errSig.Ptrans02m/fcsParams.errSig.Ptrans20m
fcsParams.errSig.Ptrans02e/fcsParams.errSig.Ptrans20e
