%Saves data for CTN noise curves into the file ctnData

clear

root = fileparts(fileparts(mfilename('fullpath')));

rowDataDir = [root filesep 'DATA' filesep 'RecordedData' filesep];          %dir of recorded data (asd & tf)
saveDataDir = [root filesep 'DATA' filesep];                                %dir where ctnData.mat is saved
%% Transfer Functions

% units: [dB]
% open loop
%%% mode00 (William servo)
h00OLMag = importdata([rowDataDir, '00tolm.asc'], ' ', 17);                 % [dB] Mag   
h00OLPh  = importdata([rowDataDir, '00told.asc'], ' ', 17);                 % [deg] Phase
%%% HOM02 (Hang servo)
h02OLMag = importdata([rowDataDir, '02tfm.asc'], ' ', 17);
h02OLPh  = importdata([rowDataDir, '02tfd.asc'], ' ', 17);
%%% HOM20 (Hang servo)
h20OLMag = importdata([rowDataDir, '20tfm.asc'], ' ', 17);
h20OLPh  = importdata([rowDataDir, '20tfd.asc'], ' ', 17);

% units: [dB]
% servo
%%% 00
h00ServoMag = importdata([rowDataDir, '00stm.asc'], ' ', 17);
h00ServoPh  = importdata([rowDataDir, '00std.asc'], ' ', 17);
%%% 02
h02ServoMag = importdata([rowDataDir, '02tfsm.asc'], ' ', 17);
h02ServoPh  = importdata([rowDataDir, '02tfsd.asc'], ' ', 17);
%%% 20
h20ServoMag = importdata([rowDataDir, '20tfsm.asc'], ' ', 17);
h20ServoPh  = importdata([rowDataDir, '20tfsd.asc'], ' ', 17);

% units: [dB]
% marconi in to mixer out
%%% 02
h02M2MMag = importdata([rowDataDir, '02m2mm.asc'], ' ', 17);
h02M2MPh  = importdata([rowDataDir, '02m2md.asc'], ' ', 17);
%%% 20
h20M2MMag = importdata([rowDataDir, '20m2mm.asc'], ' ', 17);
h20M2MPh  = importdata([rowDataDir, '20m2md.asc'], ' ', 17);


h00OLMag = h00OLMag.data;
h00OLPh = h00OLPh.data;
h02OLMag = h02OLMag.data;
h02OLPh = h02OLPh.data;
h20OLMag = h20OLMag.data;
h20OLPh = h20OLPh.data;
h00ServoMag = h00ServoMag.data;
h00ServoPh = h00ServoPh.data;
h02ServoMag = h02ServoMag.data;
h02ServoPh = h02ServoPh.data;
h20ServoMag = h20ServoMag.data;
h20ServoPh = h20ServoPh.data;
h02M2MMag = h02M2MMag.data;
h02M2MPh = h02M2MPh.data;
h20M2MMag = h20M2MMag.data;
h20M2MPh = h20M2MPh.data;


%% Servo Noise

% units: V/rtHz
% measure at the output of the fast channel
n00sServoOut = importdata([rowDataDir, '00sn.asc'], ' ', 14);               % 00
n02sServoOut = importdata([rowDataDir, '02sn.asc'], ' ', 14);               % 02
n20sServoOut = importdata([rowDataDir, '20sn.asc'], ' ', 14);               % 20

n00sServoOut = n00sServoOut.data;
n02sServoOut = n02sServoOut.data;
n20sServoOut = n20sServoOut.data;

%% PD dark noise
% units: V/rtHz
% measured after Mixer
n00PDServoIn = importdata([rowDataDir, '00pd.asc'], ' ', 14);               % 00 
n02PDServoIn = importdata([rowDataDir, '02pd.asc'], ' ', 14);               % 02
n20PDServoIn = importdata([rowDataDir, '20pd.asc'], ' ', 14);               % 20

n00PDServoIn = n00PDServoIn.data;
n02PDServoIn = n02PDServoIn.data;
n20PDServoIn = n20PDServoIn.data;

%% Marconi Noise
% units: V/rtHz
% measured at resonant freqs of HOM02/20 (phase detector ...???)
% AOM driver noise not included
n02Marconi = importdata([rowDataDir, '02mn.asc'], ' ', 14);                 % 02
n20Marconi = importdata([rowDataDir, '20mn.asc'], ' ', 14);                 % 20

n02Marconi = n02Marconi.data;
n20Marconi = n20Marconi.data;

%% Total Loop Noise 
% units: V/rtHz
% Measured at closed loop at the input and output of the servo, respectively
%%% 00
n00totServoIn = importdata([rowDataDir, '00snin.asc'], ' ', 14);    
n00totServoOut = importdata([rowDataDir, '00snout.asc'], ' ', 14);
%%% 02
n02totServoIn = importdata([rowDataDir, '02nsin.asc'], ' ', 14);
n02totServoOut = importdata([rowDataDir, '02nsout.asc'], ' ', 14);
%%% 20
n20totServoIn = importdata([rowDataDir, '20nsin.asc'], ' ', 14);
n20totServoOut = importdata([rowDataDir, '20nsout.asc'], ' ', 14);
%%% Beat Note
nHOMBeatNote = importdata([rowDataDir, 'bn.asc'], ' ', 14);
nHOMBeatNoteI = importdata([rowDataDir, 'bni.asc'], ' ', 14);
nHOMBeatNote_refl = importdata([rowDataDir, 'bn_refl.asc'], ' ', 14);
nHOMBeatNoteC = importdata([rowDataDir, 'bnc.asc'], ' ', 14);
nHOMBeatNote_air = importdata([rowDataDir, 'bn_air.asc'], ' ', 14);


n00totServoIn  = n00totServoIn.data;
n00totServoOut = n00totServoOut.data;
n02totServoIn  = n02totServoIn.data;
n02totServoOut = n02totServoOut.data;
n20totServoIn  = n20totServoIn.data;
n20totServoOut = n20totServoOut.data;
nHOMBeatNote  = nHOMBeatNote.data;
nHOMBeatNoteI  = nHOMBeatNoteI.data;
nHOMBeatNote_refl  = nHOMBeatNote_refl.data;
nHOMBeatNoteC  = nHOMBeatNoteC.data;
nHOMBeatNote_air  = nHOMBeatNote_air.data;


f = n00PDServoIn(:, 1);

save([saveDataDir, 'ctnData']);