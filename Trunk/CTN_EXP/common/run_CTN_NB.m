
%RUN_CTN_NB  Script to run the CTN NoiseBudget demo

%% Path setup

parent = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath([parent filesep 'Common' filesep 'Utils']));
addpath(genpath([parent filesep 'CTN_EXP']));

%% Load parameters, linearize the model, and extract noise terms

disp('Loading parameters for the CTN Simulink model')
ctnData;
ctnParams;
ctnNbParams;
% [noises, sys] = nbFromSimulink('ctnFCS', fcsParams.freq, 'dof', 'ServoOut_00');
% [noises, sys] = nbFromSimulink('ctnFCS', fcsParams.freq, 'dof', 'ServoOut_02');
%[noises, sys] = nbFromSimulink('ctnFCS', fcsParams.freq, 'dof', 'ServoOut_20');
% [noises, sys] = nbFromSimulink('ctnFCS', fcsParams.freq, 'dof', 'beat02');
%[noises, sys] = nbFromSimulink('ctnFCS', fcsParams.freq, 'dof', 'TF2');
% [noises, sys] = nbFromSimulink('ctnFCS', fcsParams.freq, 'dof', '00');
[noises, sys] = nbFromSimulink('ctnFCS', fcsParams.freq, 'dof', 'BeatNote');
%% Make a quick NB plot

close all

disp('Plotting noises')
nb = nbGroupNoises('ctnFCS', noises, sys);
nb.sortModel();

%%%appending measured noises
noiseCal = bode(sys(1),2*pi*fcsParams.freq);    %freq vector in rad/s
noiseCal = squeeze(noiseCal)';

totalNoiseBN.asd = fcsParams.ctnNb.beatNote ./ noiseCal;
totalNoiseBN.f = fcsParams.freq;
totalNoiseBN.name = 'Measured Total';

%%% From 11/6/2014, 0.105e-3 W power on each HOM in reflection
totalNoiseBN_refl.asd = fcsParams.ctnNb.beatNote_refl ./ noiseCal;
totalNoiseBN_refl.f = fcsParams.freq;
totalNoiseBN_refl.name = 'Measured Total Refl (11/6/2014)';

%%%coating thermal noise
coatingNoise.asd = fcsParams.ctnNb.CTN * sqrt(2);
coatingNoise.f = fcsParams.freq;
coatingNoise.name = 'Coating';

measuredNoiseBN = Noise(totalNoiseBN);
measuredNoiseBN_refl = Noise(totalNoiseBN_refl);
coatingNoise = Noise(coatingNoise);
nb.referenceNoises = {measuredNoiseBN, measuredNoiseBN_refl, coatingNoise};
matlabNoisePlot(nb);



%% bode plots
% bops = bodeoptions;
% bops.FreqUnits = 'Hz';
% bops.Xlim = [1e3, 100e3];
% 
% % TF for calibration
% figure(101)
% bodeplot(sys(1),bops)
% 
% % TFs for each noise source to noise sink
% for n = 2:length(sys.d)
%     figure(100+n)
%     bodeplot(sys(n),bops)
% end

