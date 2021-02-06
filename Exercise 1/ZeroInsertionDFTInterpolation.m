%Zero Insertion - DFT-Based Interpolation
clc;
clear all;
close all;

%(a)Load the music signal handel
load handel;
%player = audioplayer(y, Fs);
%play(player);

%(b)Compute the DFTs of the above generated signals and store them as X2,X3,X4.
N = 20000;
x = y(1:N);
x2 = x(1:2:N); % produces a sequence with even length
x3 = x(1:3:N); % produces a sequence with odd length
x4 = x(1:4:N);

%(c)Apply DFT based method for interpolation (x2 with K = 1, x3 with K = 2,
%x4 with K = 3)-Perform zero insertion
n = input('Enter a 1, 2, 3 for x2, x3, x4: ');
switch n
    case 1
        %DFT
        X2 = fft(x2);
        N2=length(X2);
        %Even Zero insertion
        N2=N2/2;
        K=1;
        X_pad = [X2(1:N2); X2(N2+1)/2; zeros(K*length(x2)-1, 1); X2(N2+1)/2; X2((N2+2):length(x2))];

    case 2
        X3 = fft(x3);
        N3=length(X3);
        %Odd
        N1=(N3+1)/2;
        K=2;
        X_pad = [X3(1:N1);zeros(K*length(x3) ,1); X3((N1+1):length(x3))];
        
    case 3
        X4 = fft(x4);
        N4=length(X4);
        %Even
        N2=N4/2;
        K=3;
        X_pad = [X4(1:N2); X4(N2+1)/2; zeros(K*length(x4)-1, 1); X4(N2+1)/2; X4((N2+2):length(x4))];
end

%Plot x and X_pad zero inserted
fn2 = 0: 0.5/256: 0.5;
a2 = abs(X_pad(1:257));
subplot(2,1,1)
switch n
    case 1
        plot(x2, '-')
    case 2
        plot(x3, '-')
    case 3
        plot(x4, '-')
end

axis ([0 256 -1.2 1.2])
legend('Original Signal N=256')
xlabel('n')
ylabel('x[n]')

subplot(2,1,2)
plot(fn2, a2, '-')
grid
axis ([0 0.6 0 25])
caption = sprintf(' K= %d', K);
title(caption, 'FontSize', 12);
xlabel('Normalized frequency fft 256')
ylabel('|X(k)|')

%(c.2)Take inverse discrete fourier transform which will give the time domain representation of the interpolated signal.
%(c.3)Rescale the amplitude of the time domain signal by multiplying with K+1
x_i=ifft(X_pad)*(K+1);

%(c.4)Take its first (K+1)(N-1) + 1 samples
switch n
    case 1
        xo =x_i(1:(K+1)*((length(x2)-1)+1));
    case 2
        xo =x_i(1:(K+1)*((length(x3)-1)+1));
    case 3
        xo =x_i(1:(K+1)*((length(x4)-1)+1));
end

%(d)Calculate the difference between original signal and interpolated signal using 2-norm
x_z = [x; zeros(100,1)];
xo_z = [xo; zeros(100,1)];
Diffnorm=norm(x_z(1:20090)-xo_z(1:20090))
       
%(e)Take first 50 samples of the original and interpolated signal and plot them together on the same graph
t=1:1:50;
figure()
plot(t,(xo(1:50)), '--')
xlabel('n')
ylabel('x[n]')
hold on
plot(x(1:50), '-')
hold off
caption = sprintf('Comparision between interpolated signal x with K= %d and original signal ', K);
title(caption, 'FontSize', 12);
legend({'Interpolated Signal N=50','Original Signal N=50'},'Location','northeast')
xlabel('time in seconds')
ylabel('signal amplitude')

player = audioplayer(xo, Fs);
play(player);

%diff=xo(1:50)-x(1:50);
%sum(diff)
