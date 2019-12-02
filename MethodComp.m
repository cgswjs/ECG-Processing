clear;clc
close all

%% Loading data
load 'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\Clean ECG and noise sources\MIT_115_clean.csv'
load 'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\Clean ECG and noise sources\MIT_116_clean.csv'
load 'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\Clean ECG and noise sources\MIT_117_clean.csv'
load 'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\Clean ECG and noise sources\MIT_234_clean.csv'
load 'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\Clean ECG and noise sources\ma.csv'
load 'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\Clean ECG and noise sources\em.csv'


fs=360;%sampling rate
A=MIT_234_clean;%clean noise
B=ma;%ma
C=em;%em

x=A(1:10*fs,1)';%processed clean data
n1=B(1:10*fs,1)';%processed noisy MA
n2=C(1:10*fs,1)';%EM

w1_1=1;%weight of ma for S1
w2_1=0.3;%weight of em for S1

w1_2=0.3;%weight of ma for S2
w2_2=1;%weight of em for S2

S1=x+w1_1*n1+w2_1*n2;%MA dom sig
S2=x+w1_2*n1+w2_2*n2;%Em dom sig

xn1=w1_1*n1+w2_1*n2;%MA dom noise
xn2=w1_2*n1+w2_2*n2;%EM dom noise
xn1_p = xn1(1:3600);
xn2_p = xn2(1:3600);

sig1_1 = highpass(S1,0.5,fs);%S1 without BW
sig1_2 = highpass(S2,0.5,fs);%S2 without BW

sig2 = highpass(x,0.5,fs);%clean ecg without baseline

LenS0 = length(sig2);
LenS1 = length(sig1_1);
LenS2 = length(sig1_2);

df0=fs/LenS0;
df1=fs/LenS1;
df2=fs/LenS2;

t0 = 0:1/fs:(LenS1-1)/fs;


fp0 = -fs/2:df0:fs/2-df0;
fp1 = -fs/2:df1:fs/2-df1;
fp2 = -fs/2:df2:fs/2-df2;


power1=fftshift(fft(x));
power2_1=fftshift(fft(xn1));
power2_2=fftshift(fft(xn2));
power_ma=fftshift(fft(S1));
power_em=fftshift(fft(S2));

power1=abs(power1)/LenS1;%power of clean signal
power2_1=abs(power2_1)/LenS1;%power of noise signal MA
power2_2=abs(power2_2)/LenS2;%power of noise signal EM
power_ma=abs(power_ma)/LenS1;
power_em=abs(power_em)/LenS1;

SNR0_1=10*log10(power1/power2_1);%SNR of corrupted signal
SNR0_2=10*log10(power1/power2_2);%SNR of corrupted signal

%preprocessed spectrum
power1pp=fftshift(fft(x(1:2500)));
power1pp=abs(power1pp)/2500;

%Power plot
% figure()
% plot(fp0,power1)
% %title('Power of Clean Signal')
% xlim([-50,50])
% ylim([0,0.05])
% xlabel('Hz')
% ylabel('db')

% figure()
% plot(fp1,power2_1)
% xlim([-50,50])
% ylim([0,0.05])
% xlabel('Hz')
% ylabel('db')
%title('Power of MA')

% figure()
% plot(fp2,power2_2)
% xlim([-50,50])
% ylim([0,0.25])
% xlabel('Hz')
% ylabel('db')
% %title('Power of EM')
%
% figure()
% plot(fp2,power_ma)
% xlim([-50,50])
% ylim([0,0.05])
% xlabel('Hz')
% ylabel('db')
%title('Power of S1')
%
% figure()
% plot(fp2,power_em)
% %title('Power of S2')
% xlim([-50,50])
% ylim([0,0.25])
% xlabel('Hz')
% ylabel('db')

%% MA dom signal
%SPEWT
f1=0.1;%1st frequency boundary for SPEWT
f2=1;%2nd frequency boundary for SPEWT
[Denozsig0_spewt,SNR0_spewt,RMSE0_spewt,SNR0]=SPEWT(fs,x,n1,n2,w1_1,w2_1,f1,f2);
% if length(Denozsig0_spewt)>length(sig2)
%     SNR0spewt=snr(sig2,Denozsig0_spewt(1:length(sig2))'-sig2);
% else
%     SNR0spewt=snr(sig2(1:length(Denozsig0_spewt)),Denozsig0_spewt'-sig2(1:length(Denozsig0_spewt)));
% end

%ADTF
[Denozsig0_ADTF,SNR0_ADTF,RMSE0_ADTF]=ADTFDWT(fs,x,n1,n2,w1_1,w2_1);
% if length(Denozsig0_ADTF)>length(sig2)
%     SNR0ADTF=snr(sig2,Denozsig0_ADTF(1:length(sig2))-sig2);
% else
%     SNR0ADTF=snr(sig2(1:length(Denozsig0_ADTF)),Denozsig0_ADTF-sig2(1:length(Denozsig0_ADTF)));
% end

%QRSDSWT
gr1=0;
%Constants
[LA WA]=size(sig1_1);%size of table A
S1s=1;%remove 1st 2s   
S1e=WA;%remove last 2s
NoD=S1e-S1s+1;%Number of input Raw data
level=9;
wavename='haar';
%SWT Decomposition Constants for OQRSDS-WT
bpt2=1;
bpt3=1;
bpt4=1;
bpt5=1;
bpt6=1;
bpt7=1;
bpt8=1;
bpt9=1;
QRSl=4;%original QRS level from 1-5
QRSs=18;%start point of QRS segment
QRSe=10;%end point of QRS segment
[swd1,swa1,swd_lp1,ecgdata1,NoDs1]=swtdecomp(sig1_1,level,wavename,fs,gr1,NoD);

%OQRSDS-WT Denoised
[Denozsig0_ma]=OQRSDSDWT(sig1_1,swd1,swa1,wavename,ecgdata1,NoDs1,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
if length(Denozsig0_ma)>length(sig2)
    SNR0_ma=snr(sig2,Denozsig0_ma(1:length(sig2))-sig2);
else
    SNR0_ma=snr(sig2(1:length(Denozsig0_ma)),Denozsig0_ma-sig2(1:length(Denozsig0_ma)));
end

if length(Denozsig0_ma)>length(sig2)%RMSE calculation
    n=length(sig2);
    MSE0_ma=1/n*sum((Denozsig0_ma(1:n)-sig2).^2);
    RMSE0_ma=sqrt(MSE0_ma);
else
    n=length(Denozsig0_ma);
    MSE0_ma=1/n*sum((Denozsig0_ma-sig2(1:n)).^2);
    RMSE0_ma=sqrt(MSE0_ma);
end

%MQRSDSDWT+ADLP
QRSl=3;%adaptive QRS level from 1-4 for MA dominated signal S1
QRSs=18;%start point of QRS segment
QRSe=10;%end point of QRS segment
%OQRSDS-WT Denoised
[Denozsig1_ma]=OQRSDSDWT(sig1_1,swd1,swa1,wavename,ecgdata1,NoDs1,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
if length(Denozsig1_ma)>length(sig2)
    SNR1_ma=snr(sig2,Denozsig1_ma(1:length(sig2))-sig2);
else
    SNR1_ma=snr(sig2(1:length(Denozsig1_ma)),Denozsig1_ma-sig2(1:length(Denozsig1_ma)));
end

%MQRSDSDWT+ADLP+LBQD
[Denozsig2_ma]=MQRSDSDWT(sig1_1,swd1,swa1,wavename,ecgdata1,NoDs1,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
if length(Denozsig2_ma)>length(sig2)
    SNR2_ma=snr(sig2,Denozsig2_ma(1:length(sig2))-sig2);
else
    SNR2_ma=snr(sig2(1:length(Denozsig2_ma)),Denozsig2_ma-sig2(1:length(Denozsig2_ma)));
end

%MQRSDSDWT+ADLP+LBQD+beta
% dbpt=0.1;
% bptcount=1;
% for bpt5=0:dbpt:1
%     for bpt6=0:dbpt:1
%         for bpt7=1:dbpt:1
%             for bpt8=0:dbpt:1
%                 for bpt9=0:dbpt:1
%                     [Denozsig3_ma]=MQRSDSDWT(sig1_1,swd1,swa1,wavename,ecgdata1,NoDs1,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
%                     Denozsig3_maM(bptcount,1)=bpt5;
%                     Denozsig3_maM(bptcount,2)=bpt6;
%                     Denozsig3_maM(bptcount,3)=bpt7;
%                     Denozsig3_maM(bptcount,4)=bpt8;
%                     Denozsig3_maM(bptcount,5)=bpt9;
%                     if length(Denozsig3_ma)>length(sig2)
%                         SNR3ma=snr(sig2,Denozsig3_ma(1:length(sig2))-sig2);
%                     else
%                         SNR3ma=snr(sig2(1:length(Denozsig3_ma)),Denozsig3_ma-sig2(1:length(Denozsig3_ma)));
%                     end
%                     Denozsig3_maM(bptcount,6)=SNR3ma;%store denosed signal
%                     bptcount=bptcount+1;
%                 end
%             end
%         end
%     end
% end

% T1=readtable('D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\BETAs\S1_beta_results.csv');
% T1=table2array(T1);
% Denozsig3_maM=T1(2:size(T1,1),:);
% maxM_ma=max(Denozsig3_maM);
% Denozsig3_ma=MQRSDSDWT(sig1_1,swd1,swa1,wavename,ecgdata1,NoDs1,fs,QRSl,...
%                        QRSs,QRSe,bpt2,bpt3,bpt4,maxM_ma(1),maxM_ma(2),...
%                        maxM_ma(3),maxM_ma(4),maxM_ma(5),gr1);%reconstruct the best case with betas
% SNR3_ma=maxM_ma(6);
%% EM dom signal
%SPEWT
f1=0.1;%1st frequency boundary for SPEWT
f2=1;%2nd frequency boundary for SPEWT
[Denozsig1_spewt,SNR1_spewt,RMSE1_spewt,SNR1]=SPEWT(fs,x,n1,n2,w1_2,w2_2,f1,f2);

%ADTF
[Denozsig1_ADTF,SNR1_ADTF,RMSE1_ADTF]=ADTFDWT(fs,x,n1,n2,w1_2,w2_2);

%OQRSDSDWT
%Constants
[LA WA]=size(sig1_2);%size of table A
S1s=1;%remove 1st 2s   
S1e=WA;%remove last 2s
NoD=S1e-S1s+1;%Number of input Raw data
level=9;
wavename='haar';
%SWT Decomposition Constants for OQRSDS-WT
bpt2=1;
bpt3=1;
bpt4=1;
bpt5=1;
bpt6=1;
bpt7=1;
bpt8=1;
bpt9=1;
QRSl=4;%original QRS level from 1-5
QRSs=18;%start point of QRS segment
QRSe=10;%end point of QRS segment
[swd2,swa2,swd_lp2,ecgdata2,NoDs2]=swtdecomp(sig1_2,level,wavename,fs,gr1,NoD);

%OQRSDS-WT
[Denozsig0_em]=OQRSDSDWT(sig1_2,swd2,swa2,wavename,ecgdata2,NoDs2,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
if length(Denozsig0_em)>length(sig2)
    SNR0_em=snr(sig2,Denozsig0_mem(1:length(sig2))-sig2);
else
    SNR0_em=snr(sig2(1:length(Denozsig0_em)),Denozsig0_em-sig2(1:length(Denozsig0_em)));
end

if length(Denozsig0_em)>length(sig2)%RMSE calculation
    n=length(sig2);
    MSE0_em=1/n*sum((Denozsig0_em(1:n)-sig2).^2);
    RMSE0_em=sqrt(MSE0_em);
else
    n=length(Denozsig0_em);
    MSE0_em=1/n*sum((Denozsig0_em-sig2(1:n)).^2);
    RMSE0_em=sqrt(MSE0_em);
end

%MQRSDSDWT+ADLP
QRSl=6;%adaptive QRS level from 1-4 for MA dominated signal S1
QRSs=15;%start point of QRS segment
QRSe=10;%end point of QRS segment
%OQRSDS-WT Denoised
[Denozsig1_em]=OQRSDSDWT(sig1_2,swd2,swa2,wavename,ecgdata2,NoDs2,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
if length(Denozsig1_ma)>length(sig2)
    SNR1_em=snr(sig2,Denozsig1_em(1:length(sig2))-sig2);
else
    SNR1_em=snr(sig2(1:length(Denozsig1_em)),Denozsig1_em-sig2(1:length(Denozsig1_em)));
end

%OQRSDSDWT+ADLP+LBQD
[Denozsig2_em]=MQRSDSDWT(sig1_2,swd2,swa2,wavename,ecgdata2,NoDs2,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
if length(Denozsig2_em)>length(sig2)
    SNR2_em=snr(sig2,Denozsig2_em(1:length(sig2))-sig2);
else
    SNR2_em=snr(sig2(1:length(Denozsig2_em)),Denozsig2_em-sig2(1:length(Denozsig2_em)));
end

%MQRSDSDWT+ADLP+LBQD+beta
% dbpt=0.1;
% bptcount=1;
% for bpt5=0:dbpt:1
%     for bpt6=0:dbpt:1
%         for bpt7=1:dbpt:1
%             for bpt8=0:dbpt:1
%                 for bpt9=0:dbpt:1
%                     [Denozsig3_em]=MQRSDSDWT(sig1_2,swd2,swa2,wavename,ecgdata2,NoDs2,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr1);
%                     Denozsig3_emM(bptcount,1)=bpt5;
%                     Denozsig3_emM(bptcount,2)=bpt6;
%                     Denozsig3_emM(bptcount,3)=bpt7;
%                     Denozsig3_emM(bptcount,4)=bpt8;
%                     Denozsig3_emM(bptcount,5)=bpt9;
%                     if length(Denozsig3_em)>length(sig2)
%                         SNR3em=snr(sig2,Denozsig3_em(1:length(sig2))-sig2);
%                     else
%                         SNR3em=snr(sig2(1:length(Denozsig3_em)),Denozsig3_em-sig2(1:length(Denozsig3_em)));
%                     end
%                     Denozsig3_emM(bptcount,6)=SNR3em;%store denosed signal
%                     bptcount=bptcount+1;
%                 end
%             end
%         end
%     end
% end
% T2=readtable('D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\DATA\BETAs\S2_beta_results.csv');
% T2=table2array(T2);
% Denozsig3_emM=T2(2:size(T2,1),:);
% maxM_em=max(Denozsig3_emM);
% Denozsig3_em=MQRSDSDWT(sig1_2,swd2,swa2,wavename,ecgdata2,NoDs2,fs,QRSl,...
%                        QRSs,QRSe,bpt2,bpt3,bpt4,maxM_em(1),maxM_em(2),...
%                        maxM_em(3),maxM_em(4),maxM_em(5),gr1);%reconstruct the best case with betas
% SNR3_em=maxM_em(6);

%% Save tables for the relationship between betas and output SNR
% T_ma = array2table(Denozsig3_maM);%convert matrix Denozsig3_maM to a table
% T_ma.Properties.VariableNames(1:6) = {'b5','b6','b7','b8','b9','SNR_ma'};%set table headers
% writetable(T_ma,'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\Methods comparison Algorithm\Packed Method Comparison\S1_beta_results.csv')%save table as csv file to a directory
% T_em = array2table(Denozsig3_emM);%convert matrix Denozsig3_maM to a table
% T_em.Properties.VariableNames(1:6) = {'b5','b6','b7','b8','b9','SNR_em'};%set table headers
% writetable(T_em,'D:\UNLV Study MicroSD Backup\Sutdy\CVR Project Denoising\Journal Paper 1\Methods comparison Algorithm\Packed Method Comparison\S2_beta_results.csv')%save table as csv file to a directory

%% Plotting
%Experimental signals comparison
figure('Name','Experimental Signals')
subplot(5,1,1)
plot(sig1_1)
title(['ma dom signal S1 with SNR = ' num2str(SNR0)])
subplot(5,1,2)
plot(xn1_p)
title('Noise Signal MA')
subplot(5,1,3)
plot(sig1_2)
title(['em dom signal S2 with SNR= ' num2str(SNR1)])
subplot(5,1,4)
plot(xn2_p)
title('Noise Signal EM')
subplot(5,1,5)
plot(sig2)
title('Clean signal ECG 115 without BW')

%--------------------------------------------------------------------------%
% %Original Comparison MA
figure('Name','Picked Methods Comparison MA')
subplot(3,1,1)
plot(Denozsig0_spewt)
% title(['(a) SPEWT Result of S1 with SER = ' num2str(SNR0_spewt) 'and RMSE = ' num2str(RMSE0_spewt)]) % 
title(['(a) SPEWT Result of S1 with SER = ' num2str(SNR0_spewt) ' and RMSE = ' num2str(RMSE0_spewt)])
subplot(3,1,2)
plot(Denozsig0_ADTF)
title(['(b) ADTF Result of S1 with SER = ' num2str(SNR0_ADTF) ' and RMSE = ' num2str(RMSE0_ADTF)]) % 
subplot(3,1,3)
plot(Denozsig0_ma)
title(['(c) QRSDS Results of S1 with SER = ' num2str(SNR0_ma) ' and RMSE = ' num2str(RMSE0_ma)])

% %Improved Methods comparison MA
% figure('Name','Improved QRSDS Comparison MA')
% subplot(3,1,1)
% plot(Denozsig1_ma)
% title(['(a) ADMP improved denoised signal with SER = ' num2str(SNR1_ma)]) % ' and RMSE = ' num2str(RMSE1_ma)]
% subplot(3,1,2)
% plot(Denozsig2_ma)
% title(['(b) LBQRS+ADMP improved denoised signal with SER = ' num2str(SNR2_ma)]) % ' and RMSE = ' num2str(RMSE2_ma)]
% subplot(3,1,3)
% plot(Denozsig3_ma)
% title(['(c) LBQRS+ADMP+EPT denoised signal with SER = ' num2str(SNR3_ma)])%' and RMSE = ' num2str(RMSE3_ma)]

%----------------------------------------------------------------------------------------------------------%
%Original Comparison EM
figure('Name','Picked Method Comparison EM')
subplot(3,1,1)
plot(Denozsig1_spewt)
title(['(a) SPEWT Result of S2 with SER = ' num2str(SNR1_spewt) ' and RMSE = ' num2str(RMSE1_spewt)])
subplot(3,1,2)
plot(Denozsig1_ADTF)
title(['(b) ADTF Result of S2 with SER = ' num2str(SNR1_ADTF) ' and RMSE = ' num2str(RMSE1_ADTF)])
subplot(3,1,3)
plot(Denozsig0_em)
title(['(c) QRSDS Results of S2 with SER = ' num2str(SNR0_em) ' and RMSE = ' num2str(RMSE0_em)])


% %Improved Method comparison EM
% figure('Name','Improved QRSDS Comparison EM')
% subplot(3,1,1)
% plot(Denozsig1_em)
% title(['(a) ADMP improved denoised signal with SER = ' num2str(SNR1_em)])%' and RMSE = ' num2str(RMSE1_em)]
% subplot(3,1,2)
% plot(Denozsig2_em)
% title(['(b) LBQRS+ADMP improved denoised signal with SER = ' num2str(SNR2_em)])% ' and RMSE = ' num2str(RMSE2_em)]
% subplot(3,1,3)
% plot(Denozsig3_em)
% title(['(c) LBQRS+ADMP+EPT denoised signal with SER = ' num2str(SNR3_em)])% ' and RMSE = ' num2str(RMSE3_em)]


% %% Useful Plots for Paper
% %MA
% figure()
% subplot(3,1,1)
% plot(sig2(1:2500))
% title('Clean ECG 115')
% ylabel('mV')
% subplot(3,1,2)
% plot(Denozsig0_ma)
% title('OQRSDS-WT MA')
% ylabel('mV')
% subplot(3,1,3)
% plot(Denozsig3_ma)
% title('MQRSDS-WT MA')
% ylabel('mV')
% 
% figure()
% subplot(2,1,1)
% plot(Denozsig0_em)
% title('OQRSDS-WT EM')
% ylabel('mV')
% subplot(2,1,2)
% plot(Denozsig3_em)
% title('MQRSDS-WT EM')
% ylabel('mV')


%% Denoised Spectrum Analysis
% %% Denoisd power spectrum
% power_odma=fftshift(fft(Denozsig0_ma));
% power_odem=fftshift(fft(Denozsig0_em));
%
% power_Idma=fftshift(fft(Denozsig3_ma));
% power_Idem=fftshift(fft(Denozsig3_em));
%
%
% LenAD=length(Denozsig0_ma);
% tad= 0:1/fs:(LenAD-1)/fs;
% dfAD=fs/LenAD;
% fpAD = -fs/2:dfAD:fs/2-dfAD;
%
% power_odma=abs(power_odma)/LenS1;
% power_odem=abs(power_odem)/LenS1;
%
% power_Idma=abs(power_Idma)/LenS1;
% power_Idem=abs(power_Idem)/LenS1;
%
%
% %Original QRSDS MA
% figure()
% plot(fpAD,power_odma)
% ylabel('db')
% xlabel('Hz')
% ylim([0 0.05])
% xlim([-50,50])
% %title('OQRSDS MA')
%
% %Orignal QRSDS EM
% figure()
% plot(fpAD,power_odem)
% ylabel('db')
% xlabel('Hz')
% ylim([0 0.05])
% xlim([-50,50])
% %title('OQRSDS EM')
%
% %Imp QRSDS MA
% figure()
% plot(fpAD,power_Idma)
% ylabel('db')
% xlabel('Hz')
% dim = [.3 .55 .05 .05];
% % annotation('rectangle',dim,'Color','red')
% ylim([0 0.05])
% xlim([-50,50])
% %title('MQRSDS MA')
%
% %Imp QRSDS EM
% figure()
% plot(fpAD,power_Idem)
% ylabel('db')
% xlabel('Hz')
% ylim([0 0.05])
% xlim([-50,50])
% title('MQRSDS MA')

% %MQRSDS S1
% figure()
% [maxY1 X1] = max(power_oIdma);
% maxX1=fpAD(X1);
%
% idmaxmax1 = find(power_oIdma == max(power_oIdma));
%
% plot(fpAD,power_oIdma,'-o','MarkerIndices',[idmaxmax1],...
%     'MarkerFaceColor','blue',...
%     'MarkerSize',8)
% xlabel('Hz')
% ylabel('dB')
% xlim([-50 50])
% ylim([0 0.05])
% title('MQRSDS S1')
%
% %MQRSDS S2
% figure()
% [maxY2 X2] = max(power_Idem);
% maxX2=fpAD(X2);
%
% idemxmax2 = find(power_Idem == max(power_Idem));
%
% plot(fpAD,power_Idem,'-o','MarkerIndices',[idemxmax2],...
%     'MarkerFaceColor','blue',...
%     'MarkerSize',8)
% xlabel('Hz')
% ylabel('dB')
% xlim([-50 50])
% ylim([0 0.05])
% title('MQRSDS EM')
%
% % Cross correlation
% NCCOR_odma=sum(power1pp.*power_odma)/sqrt(sum(power1pp.^2)*sum(power_odma.^2))
% NCCOR_odem=sum(power1pp.*power_odem)/sqrt(sum(power1pp.^2)*sum(power_odem.^2))
% NCCOR_idma=sum(power1pp.*power_Idma)/sqrt(sum(power1pp.^2)*sum(power_Idma.^2))
% NCCOR_idem=sum(power1pp.*power_Idem)/sqrt(sum(power1pp.^2)*sum(power_Idem.^2))