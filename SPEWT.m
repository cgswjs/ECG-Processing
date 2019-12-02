function [Denozsig1,SNR0spewt,RMSE,SNR0]=SPEWT(fs,x,N1,N2,w1,w2,f1,f2)
%f1 and f2 are spectrum boundary conditions for EWT


S=x+w1*N1+w2*N2;%synthesis noisy signal 115
xn=w1*N1+w2*N2;%Noisy additives

sig1 = highpass(S,0.5,fs);%noisy ecg without baseline
sig2 = highpass(x,0.5,fs);%clean ecg without baseline
LenS1 = length(sig1);

t0 = 0:1/fs:(LenS1-1)/fs;
power1=abs(fft(x)).^2/LenS1;%power of clean signal
power2=abs(fft(xn)).^2/LenS1;%power of noise signal
SNR0=10*log10(power1/power2);%SNR of corrupted signal

% figure()
% plot(t0,sig1)
% hold on
% plot(t0,sig2,'r')
% hold off
% legend('Corrupted','Clean')
%Spectrum processing
%Corrupted
SIG1=fft(sig1);%spectrum
n1 = length(sig1);          % number of samples
fr1 = (0:n1-1)*(fs/n1);     % frequency range
power1 = abs(SIG1).^2/n1;    % power of the DFT

%Clean
SIG2=fft(sig2);%spectrum
n2 = length(sig2);          % number of samples
fr2 = (0:n2-1)*(fs/n2);     % frequency range
power2 = abs(SIG2).^2/n2;    % power of the DFT

%Noise
XN=fft(xn);%spectrum
nXN = length(xn);          % number of samples
fXN = (0:nXN-1)*(fs/nXN);     % frequency range
powerXN = abs(XN).^2/nXN;    % power of the DFT

%EWT
params.globtrend='none';
params.log=0;
params.preproc='none';
params.method='locmaxmin';
params.detect='locmaxmin';
params.reg='none';
% params.lengthFilter=2;
% params.sigmaFilter=5;
params.N=2;
params.completion=0;
params.preproc='none';
params.reg='none';
params.InitBounds=[f1 f2];
[ewt,mfb]=EWT1D(sig1',params);
ewt1=ewt{1};
ewt2=ewt{2};
% figure()
% subplot(2,1,1)
% plot(ewt1)
% legend('mode I')
% subplot(2,1,2)
% plot(ewt2)
% legend('Mode II')

%Wavelet Thresholding
thrHeu = thselect(ewt1,'heursure');
thrRig = thselect(ewt1,'rigrsure');
n_ewt1=length(ewt1);
alpha=(sum(ewt1.^2)-n_ewt1)/n_ewt1;
beta=sqrt(1/n_ewt1*(log(n_ewt1))^3/log(2));
if alpha<beta
    thr=thrRig;
else  
    thr=thrHeu;
end
for i=1:n_ewt1
    if ewt1(i)>-thr && ewt1(i)<thr
        ewt1_th(i)=0;
    else
        ewt1_th(i)=ewt1(i);
    end
end
%Inverse EWT
ewt_th=cell(2,1);
ewt_th{1}=ewt1_th';
ewt_th{2}=ewt2;
% ewt_th=ewt1_th+ewt2;

Denozsig1=iEWT1D(ewt_th,mfb);

% %% Improvement evaluation
if length(Denozsig1)>length(sig2)
    SNR0spewt=snr(sig2,Denozsig1(1:length(sig2))'-sig2);
else
    SNR0spewt=snr(sig2(1:length(Denozsig1)),Denozsig1'-sig2(1:length(Denozsig1)));
end

%% MSE AND RMSE
if length(Denozsig1)>length(sig2)
    n=length(sig2);
    MSE=1/n*sum((Denozsig1(1:n)'-sig2).^2);
    RMSE=sqrt(MSE);
else
    n=length(Denozsig1);
    MSE=1/n*sum((Denozsig1'-sig2(1:n)).^2);
    RMSE=sqrt(MSE);
end
