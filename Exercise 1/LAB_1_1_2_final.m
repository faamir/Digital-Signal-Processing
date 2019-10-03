%DFT Averaging

clc;
clear all;
close all;

%(a)noise corrupted signal
filename = 'xn.mat';
s = load(filename)
% Extract array from structure.
m = s.xn;
m=transpose(m);
N=1792; J=128;
Fs = 128; T = 1/Fs;
L1 = 128; L2=256; L3=512; L4=1024; L5=1792;

%(b)Construct several subsets by taking the first 128, 256, 512, 1024, and 1792 samples 
%from sequence {x[n]}, and denote them by s(1), s(2), s(3), s(4), and s(5), respectively.
m=m(1:1792);
s1=m(1:J);s2=m(1:2*J);s3=m(1:4*J);s4=m(1:8*J);s5=m(1:14*J);

%(c)Apply DFT to each subset of samples
S1=fft(s1);S2=fft(s2);S3=fft(s3);S4=fft(s4);S5=fft(s5);
f1=0:1:J/2-1;

%DFT plots
P1 = abs(S1/L1); P1 = P1(1:L1/2+1); P1(2:end-1) = 2*P1(2:end-1); 
f = Fs*(0:(L1/2))/L1;
figure();plot(f,P1) 
title('Magnitude Spectrum of x - first 128 samples'); xlabel('f (Hz)'); ylabel('|P1(f)|');

P1 = abs(S2/L2); P1 = P1(1:L2/2+1); P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L2/2))/L2;
figure(); plot(f,P1) 
title('Magnitude Spectrum of x - first 256 samples'); xlabel('f (Hz)'); ylabel('|P1(f)|');

P1 = abs(S3/L3); P1 = P1(1:L3/2+1); P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L3/2))/L3;
figure(); plot(f,P1) 
title('Magnitude Spectrum of x - first 512 samples'); xlabel('f (Hz)'); ylabel('|P1(f)|');

P1 = abs(S4/L4); P1 = P1(1:L4/2+1); P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L4/2))/L4;
figure(); plot(f,P1) 
title('Magnitude Spectrum of x - first 1024 samples'); xlabel('f (Hz)'); ylabel('|P1(f)|');

P1 = abs(S5/L5); P1 = P1(1:L5/2+1); P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L5/2))/L5;
figure(); plot(f,P1) 
title('Magnitude Spectrum of x - first 1792 samples'); xlabel('f (Hz)'); ylabel('|P1(f)|');


% %DFT Plots
% f2=0:1:2*J/2-1;f3=0:1:4*J/2-1;f4=0:1:8*J/2-1;f5=0:1:14*J/2-1;
% figure();plot(f1,abs(S1(1:J/2)));title('Magnitude of s1: first 128')
% figure();plot(f2,abs(S2(1:J/2*2)));title('Magnitude of s1: first 256')
% figure();plot(f3,abs(S3(1:J/2*4)));title('Magnitude of s1: first 512')
% figure();plot(f4,abs(S4(1:J/2*8)));title('Magnitude of s1: first 1024')
% figure();plot(f5,abs(S5(1:J/2*14)));title('Magnitude of s1: first 1792')


%(d.1-3)apply the DFT averaging method where the length of each subset is taken to be K = 128,
%and the number of subsets is taken to be L = 14
K=128;
L=14;
v=linspace(128,128,14);
Out = mat2cell(m, 1, [v]); 
%Out{1}
X=0;
 for k = 1:L;
     X =  X+fft(Out{k})/L;
 end

%(d.4) Plot the magnitude of the spectrum 
P1 = abs(X/L1); P1 = P1(1:L1/2+1); P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L1/2))/L1;
figure();plot(f,P1) 
caption = sprintf('Restored signal (Magnitude) with L= %d', L);
title(caption, 'FontSize', 12); xlabel('f (Hz)'); ylabel('|P1(f)|');

figure()
plot(f1,abs(X(1:J/2)));
grid
caption = sprintf('Restored signal (Magnitude) with L= %d', L);
title(caption, 'FontSize', 12);
xlabel('Frequency in Hz');ylabel('|X(f)|');
figure()
%Original signal
x = ifft(X);
plot(f1,x(1:J/2))
caption = sprintf('Restored signal with L= %d', L);
title(caption, 'FontSize', 12);
