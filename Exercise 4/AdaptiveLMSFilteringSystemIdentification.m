clc;
clear all;
close all;
N = 1824;
N1=1800;
sig2 = 1e-6;
M1=21;
mu=0.01;
coeff=0.5;
NumberCases = 200;
InputLength = N + M1 -1;
w_avg = zeros(M1,1);
e_avg = zeros(N,1);
wn = zeros(M1,1);
u1 = zeros(M1,1);
y = ones(N,1);
e = ones(N,1);
%Constructing an FIR unknown system and obtaining unknown system response and desired response from it
window=hamming(M1);
h=fir1(M1-1,coeff,'low',window);
%Generation of input signal and noise signal
%a white noise sequence with zero mean and unity variance of length N can
for i = 1 : NumberCases
randn('state', i) % choose a state for subsequent construction
u = randn(N,1);
m = mean(u);
u = u - m; % to ensure a zero mean value
vr = var(u);
u = u/sqrt(vr); % to ensure an unity variance

%Generate a noise signal v(n) which has a zero mean
randn('state', i + NumberCases) % choose a state for subsequent construction
v = randn(N1,1);
m = mean(v);
v = v - m;
vr = var(v);
v = sqrt(sig2/vr)*v; % where sig2 denote the variance of v(n)
z = zeros(N-N1,1);
v = [v; z];
u = [0; u];
W1 = wn;

%adaptive system identification
for i=M1:1:N
    X=h*u(i-M1+1:1:i);
    d=X+v(i)';
    y=W1'*u(i-M1+1:1:i);
    E=d-y';
    W1=W1+mu*u(i-M1+1:1:i)*E;
    J(i)=E^2;
    W(:,i) = W1;
    
end

w_avg(:,i)=W1;
w_avg = sum(w_avg,2);
e_avg(:,i)=J;
e_avg = sum(e_avg,2);
end

e=e_avg/NumberCases;
wn=w_avg/NumberCases;
format long
h1=h
w=wn'
Norm=norm(h-w,'fro');
[H,w1] = freqz(h,1,1024);
[H1,w2] = freqz(w,1,1024);
%Plot the frequency response
figure(1)
plot(w1/pi,20*log10(abs(H)),'.')
grid
hold on
plot(w2/pi,20*log10(abs(H1)),'--')
title("Amplitude Responses of h (dotted line) and w (dashed line)")

%Plot the impulse responses
hfigure=figure(2);
t=1:1:M1;
pt1=plot(t,h);
hold on;
pt2=plot(t,w,'--');
legend(pt2, 'w');      
title("impulse response of h (solid line) and w (dashed line)")

%Plot the average squared error
figure(3)
t=1:1:N;
semilogy(t',e)
axis([0 1800 1e-7 1])
xlabel('epochs iterations')
ylabel('Ensemble Mean squared error (MSE)')
title("Error")
