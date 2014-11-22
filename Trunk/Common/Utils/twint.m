function [twint_ss] = twint(f0,Q);

% twint  state space model of a twin-T notch filter with 
% positive feedback, described by:
%
%                   s^2 + w0^2
%	T(s)  =    ----------------------
%              s^2 + (w0/Q)s + w0^2
%
%  See page 6.37 of the electronic filter design book
%
%         sys = twint(f0,Q)   returns a state-space
%         model of the filter, with f0 the resonant frequency
%         in Hertz (=w0/2pi).
%

w0 = 2*pi*f0;
sys = tf([1 0 w0^2],[1 w0/Q w0^2]);
twint_ss = ss(sys);
