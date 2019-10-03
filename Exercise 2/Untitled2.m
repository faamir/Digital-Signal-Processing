I=imread('1.png');
DD = I(200:200,:);
C=fft(DD);
n=256;
W=1:1:n;
L1=n;
Fs=2;
P1 = abs(C/L1); P1 = P1(1:L1/2+1); P1(2:end-1) = 2*P1(2:end-1); 
f = Fs*(0:(L1/2))/L1;
figure(7);
plot(f,P1);
title('Magnitude Spectrum a line Image');
whos I