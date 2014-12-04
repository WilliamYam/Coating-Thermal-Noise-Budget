%% open loop TFs

ctnParams;
ctnNbParams;
load_system('ctnFCS')

omega = 2 * pi * fcsParams.freq;

%%% 00 servo

%linearize model
io(1)=linio('ctnFCS/00 FCL/Servo/Sum', 1, 'openinput');
io(2)=linio('ctnFCS/00 FCL/Servo/Servo00', 1, 'openoutput');
setlinio('ctnFCS',io);
sysmodel = linearize('ctnFCS', io);
[magS00, phS00]=bode(sysmodel, omega);%bode freq vector in radians/sec
dbS00 = mag2db(squeeze(magS00));
phS00 = squeeze(phS00);

figure(200)

subplot(2, 1, 1)

semilogx(fcsParams.freq, [dbS00, fcsParams.TF.TFmagS00']);

grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([1e3, 100e3])
legend('Model', 'Data')
title('00 Servo')

subplot(2, 1, 2)

semilogx(fcsParams.freq, [phS00, fcsParams.TF.TFphS00'])

grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
xlim([1e3, 100e3])

%%% 02 servo

%linearize model
% io(1)=linio('ctnFCS/02 FCL/Servo/In1', 1, 'openinput');
% io(2)=linio('ctnFCS/02 FCL/Servo/Servo', 1, 'openoutput');
io(1)=linio('ctnFCS/02 FCL/IntensityNoise', 1, 'openinput');
io(2)=linio('ctnFCS/02 FCL/Sum1', 1, 'openoutput');
setlinio('ctnFCS',io);
sysmodel = linearize('ctnFCS', io);
[magS02, phS02]=bode(sysmodel, omega);%bode freq vector in radians/sec
dbS02 = mag2db(squeeze(magS02));
phS02 = squeeze(phS02);

figure(201)

subplot(2, 1, 1)

semilogx(fcsParams.freq, [dbS02, fcsParams.TF.TFmagS02']);

grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([1e3, 100e3])
legend('Model', 'Data')
title('02 Servo')

subplot(2, 1, 2)

semilogx(fcsParams.freq, [phS02, fcsParams.TF.TFphS02'])

grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
xlim([1e3, 100e3])

%%% 20 servo

%linearize model
% io(1)=linio('ctnFCS/20 FCL/Servo/In1', 1, 'openinput');
% io(2)=linio('ctnFCS/20 FCL/Servo/Servo', 1, 'openoutput');
io(1)=linio('ctnFCS/20 FCL/IntensityNoise', 1, 'openinput');
io(2)=linio('ctnFCS/20 FCL/Sum1', 1, 'openoutput');
setlinio('ctnFCS',io);
sysmodel = linearize('ctnFCS', io);
[magS20, phS20]=bode(sysmodel, omega);%bode freq vector in radians/sec
dbS20 = mag2db(squeeze(magS20));
phS20 = squeeze(phS20);

figure(202)

subplot(2, 1, 1)

semilogx(fcsParams.freq, [dbS20, fcsParams.TF.TFmagS20']);

grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([1e3, 100e3])
legend('Model', 'Data')
title('20 Servo')

subplot(2, 1, 2)

semilogx(fcsParams.freq, [phS20, fcsParams.TF.TFphS20'])

grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
xlim([1e3, 100e3])


%%% 00 loop

%linearize model
io=linio('ctnFCS/00 FCL/Sum1', 1, 'looptransfer');
setlinio('ctnFCS',io);
sysmodel = linearize('ctnFCS', io);
[magOL00, phOL00]=bode(sysmodel, omega);
dbOL00 = mag2db(squeeze(magOL00));
phOL00 = squeeze(phOL00);


figure(203)

subplot(2, 1, 1)

semilogx(fcsParams.freq, [dbOL00, fcsParams.TF.TFmag00' + 6]); %6dB needed for this OL measurement

grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([1e3, 100e3])

legend('Model', 'Data')
title('00 OL')

subplot(2, 1, 2)

semilogx(fcsParams.freq, [phOL00, fcsParams.TF.TFph00'])

grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
xlim([1e3, 100e3])

%%% 02 loop

%linearize model
io=linio('ctnFCS/02 FCL/Servo/Sum', 1, 'looptransfer');
setlinio('ctnFCS', io);
sysmodel = linearize('ctnFCS', io);
[magOL02, phOL02]=bode(sysmodel, omega);
dbOL02 = mag2db(squeeze(magOL02));
phOL02 = squeeze(phOL02);


figure(204)

subplot(2, 1, 1)

semilogx(fcsParams.freq, [dbOL02, fcsParams.TF.TFmag02']);

grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([1e3, 100e3])

legend('Model', 'Data')
title('02 OL')

subplot(2, 1, 2)

semilogx(fcsParams.freq, [phOL02, fcsParams.TF.TFph02'])

grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
xlim([1e3, 100e3])

%%% 20 loop

%linearize model
io=linio('ctnFCS/20 FCL/Sum', 1, 'looptransfer');
setlinio('ctnFCS',io);
sysmodel = linearize('ctnFCS', io);
[magOL20, phOL20]=bode(sysmodel, omega);
dbOL20 = mag2db(squeeze(magOL20));
phOL20 = squeeze(phOL20);


figure(205)

subplot(2, 1, 1)

semilogx(fcsParams.freq, [dbOL20, fcsParams.TF.TFmag20']);

grid on
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([1e3, 100e3])

legend('Model', 'Data')
title('20 OL')

subplot(2, 1, 2)

semilogx(fcsParams.freq, [phOL20, fcsParams.TF.TFph20'])

grid on
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
xlim([1e3, 100e3])



% 
% clear
% open_system('ctnFCS')

% % io(1)=linio('ctnFCS/00 FCL/Sum1',1,'openinput');
% % io(2)=linio('ctnFCS/00 FCL/Servo',1,'openoutput');
% io(1)=linio('ctnFCS/00 FCL/Sum1',1,'looptransfer');
% %io(1)=linio('ctnFCS/00 FCL/Sum1',1,'compsensitivity');
% 
% setlinio('ctnFCS',io);
% %setlinio('ctnFCS/00 FCL',io)
% systest = linearize('ctnFCS',io);
% figure(200)
% bode(systest)
% 
% 
% figure(100)
% % semilogx(fcsParams.freq,fcsParams.TF.TFmag00)
% grid on
% % hold on
% 
% % gain00        = -37.7;                         % dB measured gain at 
% % freq00        = 8e4;                            % Hz this trequency
% % fcsParams.servo00.zeros       = [10e3];         % Hz
% % fcsParams.servo00.poles       = [0,0];          % Hz
% % fcsParams.servo00.gain        = find_K(fcsParams.servo00.zeros,...
% %                                        fcsParams.servo00.poles,...
% %                                        gain00, freq00);  
% 
% k = 10^(fcsParams.servo00.gain/20);   % unitless
% G=freqresp(zpk(-2 * pi * fcsParams.servo00.zeros, -2 * pi * fcsParams.servo00.poles,k),fcsParams.freq,'Hz');
% 
% for ii= 1:size(fcsParams.freq,2)    
%     g(ii)= 20*log10(abs(G(1,1,ii))); 
% end
% subplot(2,1,1)
% semilogx(fcsParams.freq,g,'r')
% xlim([1e3,100e3])
% subplot(2,1,2)
% semilogx(fcsParams.freq,angle(squeeze(G(1,1,:)))*180/pi,'r')
% xlim([1e3,100e3])


% 
% 
% 
% 
% 
% 
% clear
% open_system('ctnFCS')
% ctnParams;
% ctnNbParams;
% %  io(1)=linio('ctnFCS/SystemNoise/Sum1',1,'openinput');
% %  io(2)=linio('ctnFCS/02 FCL/Cavity02/Sum1',1,'openoutput');
% %io(2)=linio('ctnFCS/02 FCL/Cavity02/Sum2',1,'openoutput');
% io(1)=linio('ctnFCS/02 FCL/Sum1',1,'looptransfer');
% 
% %io(1)=linio('ctnFCS/02 FCL/Cavity02/Sum1',1,'compsensitivity');   % 02 
% 
% setlinio('ctnFCS',io);
% %setlinio('ctnFCS/00 FCL',io)
% systest = linearize('ctnFCS',io);
% figure(200)
% bode(systest)
% 
% 
% 
% 
% 
% clear
% open_system('ctnFCS')
% ctnParams;
% ctnNbParams;
% %  io(1)=linio('ctnFCS/SystemNoise/Sum1',1,'openinput');
% %  io(2)=linio('ctnFCS/20 FCL/Cavity20/Sum1',1,'openoutput');
% %io(2)=linio('ctnFCS/20 FCL/Cavity20/Sum2',1,'openoutput');
% io(1)=linio('ctnFCS/20 FCL/Sum1',1,'looptransfer');
% 
% %io(1)=linio('ctnFCS/02 FCL/Cavity02/Sum1',1,'compsensitivity');   % 02 
% 
% setlinio('ctnFCS',io);
% systest = linearize('ctnFCS',io);
% figure(200)
% bode(systest)
% 
%  
% 
% 
% 
% clear
% open_system('test')
% 
%   %io(1)=linio('test/Gain1',1,'openinput');
%   %io(2)=linio('test/Gain2',1,'openoutput');
% %io(1)=linio('test/Gain',1,'sensitivity');
% io(1)=linio('test/Gain',1,'looptransfer');
% %io(1)=linio('test/Gain',1,'compsensitivity');
% %io(1)=linio('ctnFCS/02 FCL/Cavity02/Sum1',1,'compsensitivity');   % 02 
% 
% setlinio('test',io);
% %setlinio('ctnFCS/00 FCL',io)
% systest = linearize('test',io);
% figure(103)
% bode(systest)
% 
% 


 
