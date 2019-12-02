%Spectrum plot
plot(fp1,power_em,'linewidth',0.2)
xlabel('Hz')
ylabel('db')
ylim([0 0.15])
xlim([-100 100])

%ECG115
plot(sig2(1:2500),'linewidth',0.2)
ylabel('mV')
ylim([-1 2.5])

%Original vs Imp em
subplot(2,1,1)
plot(Denozsig0_em)
title('Orignal QRSDS of S2')
ylabel('mV')
subplot(2,1,2)
plot(Denozsig3_em)
title('Modified QRSDS of S2')
ylabel('mV')

%Original vs Imp ma
subplot(3,1,1)
plot(sig2(1:2500))
title('Clean ECG 115')
ylabel('mV')
subplot(3,1,2)
plot(Denozsig0_ma)
title('Orignal QRSDS of S1')
ylabel('mV')
subplot(3,1,3)
plot(Denozsig3_ma)
title('Modified QRSDS of S1')
ylabel('mV')