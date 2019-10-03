% [PSNR, SNR, error_before,error_after]=LAB_3(381,0.171,0.19,4,5)
function [PSNR, SNR, error_before,error_after]=LAB_3(N,wc1,wc2,w,b)
I=load("homework_if.mat");
I=I.homework_if;
F=load("homework.mat");
F=F.homework;
[m,n] = size(I);
J=mat2gray(I);
K=mat2gray(F);
%I = uint8(255 * mat2gray(I));
m = squeeze(I); 
whos m;
XX = m(200:200,:);
C=fft(XX);
W=1:1:n;
L1=n;
Fs=2;
P1 = abs(C/L1); P1 = P1(1:L1/2+1); P1(2:end-1) = 2*P1(2:end-1); 
f = Fs*(0:(L1/2))/L1;
figure(4);
plot(f,P1);
title('Magnitude Spectrum a line of Noisy Image');

FFT = fft2(double(I));
FFT1=fftshift(FFT);
figure(5);
imshow(log10(FFT1),[]);title('Noisy Image 2D Fouried Spectrum');
set(gcf, 'Position',  [460, 160, 500, 500])
N=N;wc1=wc1;wc2=wc2;w=w;b=b;
BC = win_fourier(N,4,[wc1 wc2],w,b);
%BC= fir1(2000,[0.1789 0.1852],'stop');

IMG_Out = conv2(BC,I);
IMG_Out(IMG_Out < 0) = 0;
IMG_Out(IMG_Out > 255) = 255;
IMG_Out(:, 1:round(N/2)-1) = [];
IMG_Out(:, n+1:end) = [];
%IMG_Out=mat2gray(IMG_Out);
CC = IMG_Out(200:200,:);
X=(fft(CC));
P1 = abs(X/L1); P1 = P1(1:L1/2+1); P1(2:end-1) = 2*P1(2:end-1); 
f = Fs*(0:(L1/2))/L1;
figure(6);
plot(f,P1);
title('Magnitude Spectrum of a line of filtered image');
%Plot
figure(7);
subplot(2,1,1);imshow(J);title('Noisy Image');
subplot(2,1,2);imshow((IMG_Out),[0,255]);title('Restored Image with Bandstop Filter');
H=double(uint8(255 * mat2gray(IMG_Out)));
J=double(uint8(255 * mat2gray(J)));
HH=double(uint8(255 * mat2gray(F)));
H=IMG_Out;
HH=F;
before=J-HH;
after=H-HH;
%FRO_out=norm(IMG_Out,'fro');
%FRO_=norm(F,'fro');
FRO_out2=norm(H,'fro');
FRO_2=norm(HH,'fro');
zigma= mean(mean((double(H) - double(HH)).^2));
PSNR = 10*log10((255.^2)./zigma);
SNR = snr(double(HH),double(H));
error_before=norm(before,'fro')/norm(HH,'fro');
error_after=norm(after,'fro')/norm(HH,'fro');
