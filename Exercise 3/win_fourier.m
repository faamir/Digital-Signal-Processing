% Inputs: N -- filter length (odd integer)
%         f_type -- 1 for lowpass, 2 for highpass, 3 for bandpass,
%                   and 4 for bandstop.
%         f_para -- omega_c for types 1 and 2, [omega_c1 omega_c2] for
%                   types 3 and 4.
%         w_type -- 1 for von Hann window, 2 for Hamming window, 3 for
%                   Blackman window and 4 for Kaiser window.
%         b -- beta for Kaiser window (optional).
% Output: h -- filter's impulse response (coefficients).
% Written by W.-S. Lu, University of Victoria
% Last modified: Dec. 18, 2002.
% Examples: h = win_fourier(71,1,0.6,2); %lowpass
%           h = win_fourier(71,2,0.5,3);  %highpass
%           h = win_fourier(71,3,[0.35 0.65],4,5); %bandpass
%           h = win_fourier(71,4,[0.3 0.7],4,5);  %bandstop

%problem 4.4 : win_fourier(25,1,0.25,2) 
%12.5/(100/2)=0.25
function h = win_fourier(N,f_type,f_para,w_type,b)

M = (N-1)/2;
n = 0:M-1;

% Generate impulse response of the ideal filter

if f_type == 1, % lowpass
   omec = pi*f_para;
   hd = sin((n-M)*omec)./((n-M)*pi);
   hd = [hd omec/pi fliplr(hd)];
elseif f_type == 2, % highpass
   omec = pi*f_para;
   hd = -sin((n-M)*omec)./((n-M)*pi);
   hd = [hd 1-omec/pi fliplr(hd)];
elseif f_type == 3, % bandpass
   omec1 = pi*f_para(1);
   omec2 = pi*f_para(2);
   hd = (sin((n-M)*omec2) - sin((n-M)*omec1))./((n-M)*pi);
   hd = [hd (omec2-omec1)/pi fliplr(hd)];
elseif f_type == 4, % bandstop
   omec1 = pi*f_para(1);
   omec2 = pi*f_para(2);
   hd = (sin((n-M)*omec1) - sin((n-M)*omec2))./((n-M)*pi);
   hd = [hd 1-(omec2-omec1)/pi fliplr(hd)];
end

% Geneate window function

if w_type == 1, % von Hann
   w = 0.5*(1-cos(pi*n/M));
   w = [w 1 fliplr(w)];
elseif w_type == 2, % Hamming
   w = 0.54 - 0.46*cos(pi*n/M);
   w = [w 1 fliplr(w)];
elseif w_type == 3, % Blackman
   w = 0.42 - 0.5*cos(pi*n/M) + 0.08*cos(2*pi*n/M);
   w = [w 1 fliplr(w)];
elseif w_type == 4, % Kaiser
   w = besseli(0,b*sqrt(1-(n/M-1).^2))/besseli(0,b);
   w = [w 1 fliplr(w)];
end

% Obtain impulse response

h = w.*hd;
figure(1)
%subplot(121)
plot(h)
axis([1 N 1.25*min(h) 1.25*max(h)])
axis('square')
grid
title('Impulse response of the FIR filter')

% Evaluate frequency response of the filter and plot

[H,w] = freqz(h,1,1024);
figure(2)
%subplot(122)
%plot(w*50/pi,20*log10(abs(H)))
plot(w/pi,20*log10(abs(H)))
%axis([0 pi -100 10])
axis('square')
grid
title('Amplitude response of the FIR filter')
xlabel('Frequency (MHz)')

