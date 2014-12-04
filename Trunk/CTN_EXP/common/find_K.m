
%
% k = findGain(z, p, g, f)
% returns correct k-param in for known gain g at the frequency f----
% 

function k = find_K(z, p, g, f)
    
gain0 = 20*log10(abs(evalfr(zpk(-2 * pi * z, -2 * pi * p, 1),2*pi*1i*f)));    %[dB]

kDB = g-gain0;              % [dB]

k = db2mag(kDB);

% %abs(evalfr(zpk(2 * pi * z, 2 * pi * p, k),2*pi*1i*100))
% 
% %%% plot 
% freq = 0:2:1600;
% k = 10^(k/20);   % unitless
% G=freqresp(zpk(-2 * pi * z, -2 * pi * p,k),freq,'Hz');
% 
% for ii= 1:size(freq,2)    
%     g(ii)= 20*log10(abs(G(1,1,ii))); 
% end
% loglog(freq,g)
% %bode(zpk(-2 * pi * z, -2 * pi * p,k))

    
    

  
  
  
  
  
  