function ctnNbParams
%ctnNbParams:  Defines noises for the ctnFCS Simulink NoiseBudget
%
% It outputs to a variable called "fcsParams" in the matlab workspace.
%% Setup and meta-parameters

fcsParams = evalin('base', 'fcsParams');

% fcsParams.ctnNb.meta.paramsFileName = mfilename('fullpath');
% fcsParams.ctnNb.meta.name = 'CTN';


%%  Data loading
load('ctnData.mat') 

%% Laser Noise

% units: Hz/rtHz
% typical value for NPRO 1/f noise

fcsParams.ctnNb.LaserFreqNoise = 100;                                       % [Hz/rtHz] @100Hz - laser frequency noise 

%% Laser Intensity Noise

% units: 1/rtHz
% Relative intensity noise (for transmission, where NPRO is not shot noise
% limited

fcsParams.ctnNb.IntensityNoise = db2mag(-135);                        % [1/rtHz] relative intensity noise 

%% PD shot noise
% units: Hz/rtHz
% estimated from the incident power on PD when cavity locked

     %%%-------definition of PD shot noise-------------------------------
     ShotNoise = @(pPD, hPlanck, fLight)(sqrt(2 * pPD * hPlanck * fLight));    % W/rtHz at PD
     %--------------------------------------------------------------------


%%% 1811 (Gaussian)
pPD00 = fcsParams.beam.PwrInc00 * (1 - fcsParams.beam.Coupling00);              % W
fcsParams.ctnNb.PDshotNoise00 = ShotNoise(pPD00, fcsParams.common.h,...
                                        fcsParams.beam.laserFreq00);
                                                                     
%%%LSCPD (HOM02)
% pPD02 = fcsParams.beam.PwrInc02 * (1 - fcsParams.beam.Coupling02);              % W
pPD02 = fcsParams.errSig.Ptrans02m;                                             %Measured Transmitted Power [W]
fcsParams.ctnNb.PDshotNoise02 = ShotNoise(pPD02, fcsParams.common.h,...
                                        fcsParams.beam.laserFreq02);
%%%LSCPD (HOM20)
% pPD20 = fcsParams.beam.PwrInc20 * (1 - fcsParams.beam.Coupling20);              % W
pPD20 = fcsParams.errSig.Ptrans02m;                                             %Measured Transmitted Power [W]
fcsParams.ctnNb.PDshotNoise20 = ShotNoise(pPD20, fcsParams.common.h,...
                                        fcsParams.beam.laserFreq20);
%%%Total HOM
% fcsParams.ctnNb.PDshotNoiseHOM = sqrt(PDshotNoise02.^2 + PDshotNoise20.^2);     % W
                                    
%% PD dark noise
% deatils in ctnData.m
fcsParams.ctnNb.n00PDOut = interp1(n00PDServoIn(:,1), n00PDServoIn(:,2), fcsParams.freq);
fcsParams.ctnNb.n02PDOut = interp1(n02PDServoIn(:,1), n02PDServoIn(:,2), fcsParams.freq);
fcsParams.ctnNb.n20PDOut = interp1(n20PDServoIn(:,1), n20PDServoIn(:,2), fcsParams.freq);

%% Servo
% deatils in ctnData.m
fcsParams.ctnNb.n00sServoOut = interp1(n00sServoOut(:,1), n00sServoOut(:,2), fcsParams.freq); 
fcsParams.ctnNb.n02sServoOut = interp1(n02sServoOut(:,1), n02sServoOut(:,2), fcsParams.freq); 
fcsParams.ctnNb.n20sServoOut = interp1(n20sServoOut(:,1), n20sServoOut(:,2), fcsParams.freq); 

%% VCO (Marconi)
% deatils in ctnData.m
%fcsParams.ctnNb.n02Marconi = n02Marconi * fcsParams.OptGain.PhaseDetectorGain;           % Hz/rtHz
fcsParams.ctnNb.n02Marconi = interp1(n02Marconi(:,1), n02Marconi(:,2), fcsParams.freq) /...
                             fcsParams.OptGain.PhaseDetectorGain02;                       % Hz/rtHz
fcsParams.ctnNb.n20Marconi = interp1(n02Marconi(:,1), n02Marconi(:,2), fcsParams.freq) /...
                             fcsParams.OptGain.PhaseDetectorGain20;                       % Hz/rtHz
                         
%% Length Noise (model)

fcsParams.ctnNb.noiseLength = 1e-14 ./ (1 + (fcsParams.freq / 100).^2);     % m / rtHz cavity length noise

%%CTN
fcsParams.ctnNb.CTN = 7.2e-18 .* sqrt(100 ./ fcsParams.freq) * 100 /...
    (fcsParams.beam.waist * 1e6);                                              % [m/rtHz] predicted coating thermal noise

%% Gas Phase Noise

%phase noise due to forward scattering under the assumption that the mean
%free path of the gas is much smaller than the beam waist. Based on a
%document written by Weiss, "Considerations in operating an
%interferometric..."

fcsParams.ctnNb.noiseGas = 8 * pi * fcsParams.common.alpha / fcsParams.common.wavelength *...
    sqrt(2 * fcsParams.cavity.Length * fcsParams.common.nDensity / fcsParams.common.diffusion) ./ ...
    sqrt(1 + fcsParams.beam.waist^4 * fcsParams.freq.^2 / fcsParams.common.diffusion);  % [rad / rtHz]
                         
%% Total Noise
% deatils in ctnData.m

%%%Gaussian
fcsParams.ctnNb.nTotServoIn00 = interp1(n00totServoIn(:,1), n00totServoIn(:,2), fcsParams.freq);
fcsParams.ctnNb.nTotServoOut00 = interp1(n00totServoOut(:,1), n00totServoOut(:,2), fcsParams.freq);
%%%HOM02
fcsParams.ctnNb.nTotServoIn02 = interp1(n02totServoIn(:,1), n02totServoIn(:,2), fcsParams.freq);
fcsParams.ctnNb.nTotServoOut02 = interp1(n02totServoOut(:,1), n02totServoOut(:,2), fcsParams.freq);
%%%HOM20
fcsParams.ctnNb.nTotServoIn20 = interp1(n20totServoIn(:,1), n20totServoIn(:,2), fcsParams.freq);
fcsParams.ctnNb.nTotServoOut20 = interp1(n20totServoOut(:,1), n20totServoOut(:,2), fcsParams.freq);
%%%BeatNote
% fcsParams.ctnNb.beatNote = interp1(nHOMBeatNoteC(:,1), nHOMBeatNoteC(:,2), fcsParams.freq) *...
%     2 * fcsParams.OptGain.MarconiGain02;
fcsParams.ctnNb.beatNote = interp1(nHOMBeatNote(:,1), nHOMBeatNote(:,2), fcsParams.freq) /...
    fcsParams.OptGain.PhaseDetectorGainBN;
% fcsParams.ctnNb.beatNote_refl = interp1(nHOMBeatNote_refl(:,1), nHOMBeatNote_refl(:,2), fcsParams.freq) /...
%     fcsParams.OptGain.PhaseDetectorGainBN_refl;
fcsParams.ctnNb.beatNote_air = interp1(nHOMBeatNote_air(:,1), nHOMBeatNote_air(:,2), fcsParams.freq) /...
    fcsParams.OptGain.PhaseDetectorGainBN_air;

%% Output

assignin('base', 'fcsParams', fcsParams);