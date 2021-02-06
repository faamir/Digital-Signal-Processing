%DFT Averaging
 
clc;
clear all;
close all;
 
filename = 'xn.mat';
s = load(filename)
% Extract array from structure.
m = s.xn;
m=transpose(m);
 
%Construct several subsets by taking the first 128, 256, 512, 1024, and 1792 samples 
%from sequence {x[n]}, and denote them by s(1), s(2), s(3), s(4), and s(5), respectively.
m=m(1:1792);
s1=m(1:128);s2=m(1:256);s3=m(1:512);s4=m(1:1024);s5=m(1:1792);
%Apply DFT to each subset of samples
Xs1=fft(s1);Xs2=fft(s2);Xs3=fft(s3);Xs4=fft(s4);Xs5=fft(s5);
f1=1:1:128;f2=1:1:256;f3=1:1:512;f4=1:1:1024;f5=1:1:1792;
%Plots
figure();plot(f1,abs(Xs1(1:128)));title('Magnitude of s1: first 128')
figure();plot(f2,abs(Xs2(1:128*2)));title('Magnitude of s1: first 256')
figure();plot(f3,abs(Xs3(1:128*4)));title('Magnitude of s1: first 512')
figure();plot(f4,abs(Xs4(1:128*8)));title('Magnitude of s1: first 1024')
figure();plot(f5,abs(Xs5(1:128*14)));title('Magnitude of s1: first 1792')
 
 
%apply the DFT averaging method where the length of each subset is taken to be K = 128,
%and the number of subsets is taken to be L = 14
K=128;
L=14;
v=linspace(1792/L,1792/L,L)
%v=linspace(1792/L,1792/L,L)
Out = mat2cell(m, 1, [v]); 
%Out{1}
 
X=0;
 for k = 1:L;
     X = X + fft(Out{k})/L;
 end
 
fx=1:1:1792/L
figure()
plot(fx,abs(X(1:1792/L)));
grid
caption = sprintf('Restored signal (Magnitude) with L= %d', L);
title(caption, 'FontSize', 12);
xlabel('Frequency in Hz')
figure()
%accurate number of original signal
x = ifft(X)
plot(fx,x)
caption = sprintf('Restored signal with L= %d', L);
title(caption, 'FontSize', 12);

