%Image compression using DCT method 
%[peaksnr, pzeros, I_comp, Time, ORIG_filesize, COMPRESS_filesize] = dct_compression(1,10)
function [peaksnr, pzeros, I_comp, Time, ORIG_filesize, COMPRESS_filesize] = dct_compression(n,level) 
%n: 1.cameraman, 2.boat, 3.goldhill, 4.peppers and 5.lena
%level:  enter a quality level;

%(a) Download the following test images from the course website
%camera256.mat, boat512.mat, goldhill512.mat, and peppers512.mat
%Load images
switch n
    case 1
        xn=load('camera256.mat');
        I=xn.camera256;
    case 2
        xn=load('boat512.mat');
        I=xn.boat512;   
    case 3
        xn=load('goldhill512.mat');
        I=xn.goldhill512;
    case 4
        xn=load('peppers512.mat');
        I=xn.peppers512;
    case 5
        xn=load('lena256.mat');
        I=xn.lena256;
   
    
    
        
end

%Plot original image
J=mat2gray(I);
disp('Original Image entropy: ');disp(entropy(J));
figure(1)
%subplot(1,2,1);
imagesc(I)
title("Original Image");
colormap(gray);
%(b) Perform DCT compression as follows:
%(1) First generate an all zero matrix with size same as the test image. 
%This will store the final compressed image. Let this be I_comp.
%Get size of image
[M,N]=size(I);
m=M./8;
n=N./8;
I_comp=zeros(M,N);
%Generate a variable to store the quantization matrix Q50
%Use a variable eg. Level to store the desired quality level for the compressed image level=(10,40 or 50)
Q50 = [16, 11, 10, 16, 24,  40,  51,  61;
      12, 12, 14, 19, 26,  58,  60,  55;
      14, 13, 16, 24, 40,  57,  69,  56;
      14, 17, 22, 29, 51,  87,  80,  62;
      18, 22, 37, 56, 68,  109, 103, 77;
      24, 35, 55, 64, 81,  104, 113, 92;
      49, 64, 78, 87, 103, 121, 120, 101;
      72, 92, 95, 98, 112, 100, 103, 99];
level;
tic
%Split the test image I into 8 x 8 block matrices
for i=1:m
    for j=1:n
        tStart = tic;
        %extract 8*8 block
        B=I((i-1)*8+1:(i-1)*8+8,(j-1)*8+1:(j-1)*8+8); 
        %Step 1: Level-Off and 2-D DCT
        dB=double(B);
        %B is leveled off by subtracting 128 from each entry
        Bt=dB-128; 
        %Applying 2-D DCT
        C=dct2(Bt); 
        
        %Step 2: Quantization
        %T is a scaling factor determined by the quality level
        %Compute the scaling factor for quantization
        if level < 50
            T = floor(50 ./ level);
        else
            T = (100-level) ./ 50;
        end
        %use this scaling factor to get the quantization matrix for the desired quality level
        %The scaled quantization matrix is then rounded and clipped to integer values between 0 and 255.
        %Perform Quantization by pointwise dividing C by the quantization matrix Q obtained in step 3)
        %and rounding off the resulting matrix. Let this matrix be S.
        Q=double(T.*Q50);
        S=round(C./Q);
        %put every computed block in image matrix for compressed
        compressed((i-1)*8+1:(i-1)*8+8,(j-1)*8+1:(j-1)*8+8)=S;
              
        %Step 4: Decompression
        %Pointwise multiplication of matrix S with the quantization matrix Q to get an image block Ct = Q*S
        Ct=S.*Q; 
        %Apply 2-D inverse DCT to matrix Ct
        Bt=idct2(Ct); 
        %The effect of the “level-off” operation in Step 1 is taken into account by adding 128 to the entries of matrix Bt obtained above
        Bt=Bt+128; 
        %Place each processed block in its correct location
        I_comp((i-1)*8+1:(i-1)*8+8,(j-1)*8+1:(j-1)*8+8)=Bt;
        %get zeros of S
        Rzero((i-1)*8+1:(i-1)*8+8,(j-1)*8+1:(j-1)*8+8)=S;
    end
end
toc
Time = toc
%compute zeros
%idx=Rzero==0;
%pzeros=(sum(idx(:))/(M*N))*100
pzeros=nnz(~Rzero)/(M*N)*100;
%plot compressed image
figure(2)
%subplot(2,1,2);
imagesc(I_comp);
colormap(gray);
caption = sprintf('Compressed Image, quality level= %d, zeros=%.2f', level, pzeros);
title(caption);
disp('Changed entropy: ');disp(entropy(mat2gray(compressed)));
%Compute peak signal-to-noise ratio (PSNR)
zigma= mean(mean((double(I_comp) - double(I)).^2));
peaksnr = 10*log10((255.^2)./zigma)
%peaksnr=psnrs(I_comp, I)
%figure(2);subplot(121);imhist(J,255);title('histogram of original Image');
%figure(2);subplot(122);imhist(mat2gray(I_comp),255);title('histogram of Compressed Image');
%image_contrast = max(I(:)) - min(I(:))

m = squeeze(I); 
XX = m(128:128,:);
YY = m(129:129,:);
E1=sum(abs(XX).^2);
figure(6)
subplot(221)
plot(XX)
title('8-bit gray level of 128th row')
%Plot auto correlation
cor = xcorr(XX,YY, 'biased');
subplot(222)
plot(abs(cor))
%axis([256 500 0 10000])
title('Normalized auto correlation of signal ')


dct1 = dct(XX);
subplot(223)
E=sum(abs(dct1).^2);
plot(dct1)
title('1-D DCT of the 128 row of image')

% cor2 = xcorr(dct1,dct1,'unbiased');
% subplot(224)
% plot(cor2)
% axis([256 500 -3000 10000])
% title('Normalized auto correlation of DCT ')
%Plot Magnitude response
%whos decompressed
DD = I(128:128,:);
C=fft(DD);
W=1:1:n;
L1=n;
Fs=2;
P1 = abs(C/L1); P1 = P1(1:L1/2+1); P1(2:end-1) = 2*P1(2:end-1); 
f = Fs*(0:(L1/2))/L1;
figure(7);
plot(f,P1, "*");
title('Magnitude Spectrum 128th row of Image');
imwrite(uint8(I), 'ORIG.jpeg');
imwrite(uint8(compressed), 'COMPRESS.jpeg');
FileInfoO = dir([pwd,'\ORIG.jpeg']);
ORIG_filesize = FileInfoO.bytes;
FileInfoC = dir([pwd,'\COMPRESS.jpeg']);
COMPRESS_filesize = FileInfoC.bytes;
%5
%3
%flat freq+high time+so much correlation
