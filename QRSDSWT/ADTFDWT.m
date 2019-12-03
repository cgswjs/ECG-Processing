function [Y,SNR0ADTF,RMSE]=ADTFDWT(fs,x,N1,N2,w1,w2)
S=x+w1*N1+w2*N2;%synthesis noisy signal 115
xn=w1*N1+w2*N2;%noise signal

sig1 = highpass(S,0.5,fs);%noisy ecg without baseline
sig2 = highpass(x,0.5,fs);%clean ecg without baseline
LenS1 = length(sig1);

power1=abs(fft(x)).^2/LenS1;%power of clean signal
power2=abs(fft(xn)).^2/LenS1;%power of noise signal
SNR0=10*log10(power1/power2);%SNR of corrupted signal
%% Step 1 Eliminate D1 and D2 in DWT
ecgdata=sig1;
% divide into 9 levels
wavename='db6';
level=9;
[C,L]=wavedec(ecgdata,level,wavename);

%approx coefficients at level 9
A9=appcoef(C,L,wavename,9);
%detailed coefficients for each level
D1=detcoef(C,L,1);
D2=detcoef(C,L,2);
D3=detcoef(C,L,3);
D4=detcoef(C,L,4);
D5=detcoef(C,L,5);
D6=detcoef(C,L,6);
D7=detcoef(C,L,7);
D8=detcoef(C,L,8);
D9=detcoef(C,L,9);
D1=D1.*0;
D2=D2.*0;

ecgdata=waverec(C,L,wavename);
%% Step 2 Dual thresholding algorithm for ST wave
beta=0.75;%coefficient control factor
m=15;%window length due
n=length(ecgdata);
count=1;
for i=1:length(ecgdata)/m    
    S_block{i}(1:m)=ecgdata(count:count+m-1);%window
    g(i)=mean(S_block{i});%mean of window
    Mx(i)=max(S_block{i});%window max
    Mi(i)=min(S_block{i});%window min
    Ht(i)= g(i)+((Mx(i)-g(i))*beta);%window high th
    Lt(i)= g(i)-((g(i)-Mi(i))*beta);%window low th
    
    %check each sample in this window
    for j=count:count+m-1
        if ecgdata(j)>Lt(i) && ecgdata(j)<Ht(i)
            ecgrec(j)=ecgdata(j);%noisy free if lies in the range
        else
            ecgrec(j)=g(i);%noisy otherwise
        end
    end
    count=count+m;
end

%% Step 3 Peak correction
t3=0.1;
m3=t3*fs;%window length in step 3
phi=max(abs(ecgrec));%maxima of ecgrec absolute
alpha=phi*0.1;%alpha coefficient
n3=10;%QRS correct window
for i=1:length(ecgrec)
        if abs(ecgrec(i))>=alpha
            if i>n3 && i+n3<n
                Y(i-n3:i+n3)=sig1(i-n3:i+n3);%set peaks to origin
            elseif i>n3 && i+n3>n
                Y(i-n3:n)=sig1(i-n3:n);
            elseif i<n3
                Y(1:i+n3)=sig1(1:i+n3);
            end
        else
            Y(i)=ecgrec(i);%leave non-peak
        end
end

if length(Y)>length(sig2)
    SNR0ADTF=snr(sig2,Y(1:length(sig2))-sig2);
else
    SNR0ADTF=snr(sig2(1:length(Y)),Y-sig2(1:length(Y)));
end

%% MSE AND RMSE
if length(Y)>length(sig2)
    n=length(sig2);
    MSE=1/n*sum((Y(1:n)-sig2).^2);
    RMSE=sqrt(MSE);
else
    n=length(Y);
    MSE=1/n*sum((Y-sig2(1:n)).^2);
    RMSE=sqrt(MSE);
end


