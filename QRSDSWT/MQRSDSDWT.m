function [Denozsig1]=MQRSDSDWT(sig1,swd,swa,wavename,ecgdata,NoDs,fs,QRSl,QRSs,QRSe,bpt2,bpt3,bpt4,bpt5,bpt6,bpt7,bpt8,bpt9,gr)
%sig1-input noisy sig1
%swd-wavelet detailed coeff
%swa-wavelet approx coeff
%wavename-wavelet form
%ecgdata-initialized sig1
%NoDs-number of data points
%fs-sampling  rate
%QRSl-how many levels to use QRS detection 
%QRSs-QRS segment start point
%QRSe-QRS segment end point
%bpts-beta coefficients for each level
%gr-plot control 

D1=swd(1,:);
D2=swd(2,:);
D3=swd(3,:);
D4=swd(4,:);
D5=swd(5,:);
D6=swd(6,:);
D7=swd(7,:);
D8=swd(8,:);
D9=swd(9,:);

%QRS detection
if QRSl==1%em
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
elseif QRSl==2
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
elseif QRSl==3
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
    %level 3
    [qrs_amp_raw3,qrs_i_raw3,ecg_h3]=pan_tompkinR(D3,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw3)
        D3(qrs_i_raw3(i)-QRSs:qrs_i_raw3(i)+QRSe)=0;
    end
elseif QRSl==4
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
    %level 3
    [qrs_amp_raw3,qrs_i_raw3,ecg_h3]=pan_tompkinR(D3,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw3)
        D3(qrs_i_raw3(i)-QRSs:qrs_i_raw3(i)+QRSe)=0;
    end
    %level 4
    [qrs_amp_raw4,qrs_i_raw4,ecg_h4]=pan_tompkinR(D4,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw4)
        D4(qrs_i_raw4(i)-QRSs:qrs_i_raw4(i)+QRSe)=0;
    end
elseif QRSl==5
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
    %level 3
    [qrs_amp_raw3,qrs_i_raw3,ecg_h3]=pan_tompkinR(D3,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw3)
        D3(qrs_i_raw3(i)-QRSs:qrs_i_raw3(i)+QRSe)=0;
    end
    %level 4
    [qrs_amp_raw4,qrs_i_raw4,ecg_h4]=pan_tompkinR(D4,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw4)
        D4(qrs_i_raw4(i)-QRSs:qrs_i_raw4(i)+QRSe)=0;
    end
    %level 5
    [qrs_amp_raw5,qrs_i_raw5,ecg_h5]=pan_tompkinR(D5,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw5)
        D5(qrs_i_raw5(i)-QRSs:qrs_i_raw5(i)+QRSe)=0;
    end
    
elseif QRSl==6
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
    %level 3
    [qrs_amp_raw3,qrs_i_raw3,ecg_h3]=pan_tompkinR(D3,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw3)
        D3(qrs_i_raw3(i)-QRSs:qrs_i_raw3(i)+QRSe)=0;
    end
    %level 4
    [qrs_amp_raw4,qrs_i_raw4,ecg_h4]=pan_tompkinR(D4,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw4)
        D4(qrs_i_raw4(i)-QRSs:qrs_i_raw4(i)+QRSe)=0;
    end
    %level 5
    [qrs_amp_raw5,qrs_i_raw5,ecg_h5]=pan_tompkinR(D5,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw5)
        D5(qrs_i_raw5(i)-QRSs:qrs_i_raw5(i)+QRSe)=0;
    end
    %level 6
    [qrs_amp_raw6,qrs_i_raw6,ecg_h6]=pan_tompkinR(D6,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw6)
        D6(qrs_i_raw6(i)-QRSs:qrs_i_raw6(i)+QRSe)=0;
    end
    
elseif QRSl==7
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
    %level 3
    [qrs_amp_raw3,qrs_i_raw3,ecg_h3]=pan_tompkinR(D3,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw3)
        D3(qrs_i_raw3(i)-QRSs:qrs_i_raw3(i)+QRSe)=0;
    end
    %level 4
    [qrs_amp_raw4,qrs_i_raw4,ecg_h4]=pan_tompkinR(D4,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw4)
        D4(qrs_i_raw4(i)-QRSs:qrs_i_raw4(i)+QRSe)=0;
    end
    %level 5
    [qrs_amp_raw5,qrs_i_raw5,ecg_h5]=pan_tompkinR(D5,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw5)
        D5(qrs_i_raw5(i)-QRSs:qrs_i_raw5(i)+QRSe)=0;
    end
    %level 6
    [qrs_amp_raw6,qrs_i_raw6,ecg_h6]=pan_tompkinR(D6,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw6)
        D6(qrs_i_raw6(i)-QRSs:qrs_i_raw6(i)+QRSe)=0;
    end
    %level 7
    [qrs_amp_raw7,qrs_i_raw7,ecg_h7]=pan_tompkinR(D7,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw7)
        D7(qrs_i_raw7(i)-QRSs:qrs_i_raw7(i)+QRSe)=0;
    end
elseif QRSl==8
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
    %level 3
    [qrs_amp_raw3,qrs_i_raw3,ecg_h3]=pan_tompkinR(D3,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw3)
        D3(qrs_i_raw3(i)-QRSs:qrs_i_raw3(i)+QRSe)=0;
    end
    %level 4
    [qrs_amp_raw4,qrs_i_raw4,ecg_h4]=pan_tompkinR(D4,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw4)
        D4(qrs_i_raw4(i)-QRSs:qrs_i_raw4(i)+QRSe)=0;
    end
    %level 5
    [qrs_amp_raw5,qrs_i_raw5,ecg_h5]=pan_tompkinR(D5,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw5)
        D5(qrs_i_raw5(i)-QRSs:qrs_i_raw5(i)+QRSe)=0;
    end
    %level 6
    [qrs_amp_raw6,qrs_i_raw6,ecg_h6]=pan_tompkinR(D6,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw6)
        D6(qrs_i_raw6(i)-QRSs:qrs_i_raw6(i)+QRSe)=0;
    end
    %level 7
    [qrs_amp_raw7,qrs_i_raw7,ecg_h7]=pan_tompkinR(D7,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw7)
        D7(qrs_i_raw7(i)-QRSs:qrs_i_raw7(i)+QRSe)=0;
    end
    %level 8
    [qrs_amp_raw8,qrs_i_raw8,ecg_h8]=pan_tompkinR(D8,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw8)
        D8(qrs_i_raw8(i)-QRSs:qrs_i_raw8(i)+QRSe)=0;
    end
    
elseif QRSl==9
    %level 1
    [qrs_amp_raw1,qrs_i_raw1,ecg_h1]=pan_tompkinR(D1,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw1)
        D1(qrs_i_raw1(i)-QRSs:qrs_i_raw1(i)+QRSe)=0;
    end
    %level 2
    [qrs_amp_raw2,qrs_i_raw2,ecg_h2]=pan_tompkinR(D2,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw2)
        D2(qrs_i_raw2(i)-QRSs:qrs_i_raw2(i)+QRSe)=0;
    end
    %level 3
    [qrs_amp_raw3,qrs_i_raw3,ecg_h3]=pan_tompkinR(D3,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw3)
        D3(qrs_i_raw3(i)-QRSs:qrs_i_raw3(i)+QRSe)=0;
    end
    %level 4
    [qrs_amp_raw4,qrs_i_raw4,ecg_h4]=pan_tompkinR(D4,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw4)
        D4(qrs_i_raw4(i)-QRSs:qrs_i_raw4(i)+QRSe)=0;
    end
    %level 5
    [qrs_amp_raw5,qrs_i_raw5,ecg_h5]=pan_tompkinR(D5,fs,gr);%R peak detection
    for i=2:length(qrs_i_raw5)
        D5(qrs_i_raw5(i)-QRSs:qrs_i_raw5(i)+QRSe)=0;
    end
    %level 6
    [qrs_amp_raw6,qrs_i_raw6,ecg_h6]=pan_tompkinR(D6,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw6)
        D6(qrs_i_raw6(i)-QRSs:qrs_i_raw6(i)+QRSe)=0;
    end
    %level 7
    [qrs_amp_raw7,qrs_i_raw7,ecg_h7]=pan_tompkinR(D7,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw7)
        D7(qrs_i_raw7(i)-QRSs:qrs_i_raw7(i)+QRSe)=0;
    end
    %level 8
    [qrs_amp_raw8,qrs_i_raw8,ecg_h8]=pan_tompkinR(D8,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw8)
        D8(qrs_i_raw8(i)-QRSs:qrs_i_raw8(i)+QRSe)=0;
    end
    %level 9
    [qrs_amp_raw9,qrs_i_raw9,ecg_h9]=pan_tompkinR(D9,fs,gr);%R peak detection
    for i=1:length(qrs_i_raw9)
        D8(qrs_i_raw9(i)-QRSs:qrs_i_raw9(i)+QRSe)=0;
    end
end


%PT wave removal
% nos_PT=512;%number of samples in 1sec window
% temp=[1:fs:length(sig1)];
% dnos=temp(length(temp));
% count2=0;
nos_PT=512;%number of samples in 1sec window
temp=[1:fs:length(D1)];
nose=temp(length(temp));
count2=0;
for i=1:nos_PT:nose-fs
    count2=count2+1;
    if QRSl==1%em
        PT2{count2}=D2(i:i+nos_PT-1);
        Dmax_2loc(count2)=max(PT2{count2});
        Dmin_2loc(count2)=min(PT2{count2});
        Dmax2=median(Dmax_2loc);
        Dmin2=median(Dmin_2loc);
        
        PT3{count2}=D3(i:i+nos_PT-1);
        Dmax_3loc(count2)=max(PT3{count2});
        Dmin_3loc(count2)=min(PT3{count2});
        Dmax3=median(Dmax_3loc);
        Dmin3=median(Dmin_3loc);
        
        PT4{count2}=D4(i:i+nos_PT-1);
        Dmax_4loc(count2)=max(PT4{count2});
        Dmin_4loc(count2)=min(PT4{count2});
        Dmax4=median(Dmax_4loc);
        Dmin4=median(Dmin_4loc);
        
        PT5{count2}=D5(i:i+nos_PT-1);
        Dmax_5loc(count2)=max(PT5{count2});
        Dmin_5loc(count2)=min(PT5{count2});
        Dmax5=median(Dmax_5loc);
        Dmin5=median(Dmin_5loc);
        
        PT6{count2}=D6(i:i+nos_PT-1);
        Dmax_6loc(count2)=max(PT6{count2});
        Dmin_6loc(count2)=min(PT6{count2});
        Dmax6=median(Dmax_6loc);
        Dmin6=median(Dmin_6loc);
        
        PT7{count2}=D7(i:i+nos_PT-1);
        Dmax_7loc(count2)=max(PT7{count2});
        Dmin_7loc(count2)=min(PT7{count2});
        Dmax7=median(Dmax_7loc);
        Dmin7=median(Dmin_7loc);
        
        PT8{count2}=D8(i:i+nos_PT-1);
        Dmax_8loc(count2)=max(PT8{count2});
        Dmin_8loc(count2)=min(PT8{count2});
        Dmax8=median(Dmax_8loc);
        Dmin8=median(Dmin_8loc);
        
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
    elseif QRSl==2
        PT3{count2}=D3(i:i+nos_PT-1);
        Dmax_3loc(count2)=max(PT3{count2});
        Dmin_3loc(count2)=min(PT3{count2});
        Dmax3=median(Dmax_3loc);
        Dmin3=median(Dmin_3loc);
        
        PT4{count2}=D4(i:i+nos_PT-1);
        Dmax_4loc(count2)=max(PT4{count2});
        Dmin_4loc(count2)=min(PT4{count2});
        Dmax4=median(Dmax_4loc);
        Dmin4=median(Dmin_4loc);
        
        PT5{count2}=D5(i:i+nos_PT-1);
        Dmax_5loc(count2)=max(PT5{count2});
        Dmin_5loc(count2)=min(PT5{count2});
        Dmax5=median(Dmax_5loc);
        Dmin5=median(Dmin_5loc);
        
        PT6{count2}=D6(i:i+nos_PT-1);
        Dmax_6loc(count2)=max(PT6{count2});
        Dmin_6loc(count2)=min(PT6{count2});
        Dmax6=median(Dmax_6loc);
        Dmin6=median(Dmin_6loc);
        
        PT7{count2}=D7(i:i+nos_PT-1);
        Dmax_7loc(count2)=max(PT7{count2});
        Dmin_7loc(count2)=min(PT7{count2});
        Dmax7=median(Dmax_7loc);
        Dmin7=median(Dmin_7loc);
        
        PT8{count2}=D8(i:i+nos_PT-1);
        Dmax_8loc(count2)=max(PT8{count2});
        Dmin_8loc(count2)=min(PT8{count2});
        Dmax8=median(Dmax_8loc);
        Dmin8=median(Dmin_8loc);
        
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
    elseif QRSl==3
        PT4{count2}=D4(i:i+nos_PT-1);
        Dmax_4loc(count2)=max(PT4{count2});
        Dmin_4loc(count2)=min(PT4{count2});
        Dmax4=median(Dmax_4loc);
        Dmin4=median(Dmin_4loc);
        
        PT5{count2}=D5(i:i+nos_PT-1);
        Dmax_5loc(count2)=max(PT5{count2});
        Dmin_5loc(count2)=min(PT5{count2});
        Dmax5=median(Dmax_5loc);
        Dmin5=median(Dmin_5loc);
        
        PT6{count2}=D6(i:i+nos_PT-1);
        Dmax_6loc(count2)=max(PT6{count2});
        Dmin_6loc(count2)=min(PT6{count2});
        Dmax6=median(Dmax_6loc);
        Dmin6=median(Dmin_6loc);
        
        PT7{count2}=D7(i:i+nos_PT-1);
        Dmax_7loc(count2)=max(PT7{count2});
        Dmin_7loc(count2)=min(PT7{count2});
        Dmax7=median(Dmax_7loc);
        Dmin7=median(Dmin_7loc);
        
        PT8{count2}=D8(i:i+nos_PT-1);
        Dmax_8loc(count2)=max(PT8{count2});
        Dmin_8loc(count2)=min(PT8{count2});
        Dmax8=median(Dmax_8loc);
        Dmin8=median(Dmin_8loc);
        
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
    elseif QRSl==4
        PT5{count2}=D5(i:i+nos_PT-1);
        Dmax_5loc(count2)=max(PT5{count2});
        Dmin_5loc(count2)=min(PT5{count2});
        Dmax5=median(Dmax_5loc);
        Dmin5=median(Dmin_5loc);
        
        PT6{count2}=D6(i:i+nos_PT-1);
        Dmax_6loc(count2)=max(PT6{count2});
        Dmin_6loc(count2)=min(PT6{count2});
        Dmax6=median(Dmax_6loc);
        Dmin6=median(Dmin_6loc);
        
        PT7{count2}=D7(i:i+nos_PT-1);
        Dmax_7loc(count2)=max(PT7{count2});
        Dmin_7loc(count2)=min(PT7{count2});
        Dmax7=median(Dmax_7loc);
        Dmin7=median(Dmin_7loc);
        
        PT8{count2}=D8(i:i+nos_PT-1);
        Dmax_8loc(count2)=max(PT8{count2});
        Dmin_8loc(count2)=min(PT8{count2});
        Dmax8=median(Dmax_8loc);
        Dmin8=median(Dmin_8loc);
        
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
        
    elseif QRSl==5
        PT6{count2}=D6(i:i+nos_PT-1);
        Dmax_6loc(count2)=max(PT6{count2});
        Dmin_6loc(count2)=min(PT6{count2});
        Dmax6=median(Dmax_6loc);
        Dmin6=median(Dmin_6loc);
        
        PT7{count2}=D7(i:i+nos_PT-1);
        Dmax_7loc(count2)=max(PT7{count2});
        Dmin_7loc(count2)=min(PT7{count2});
        Dmax7=median(Dmax_7loc);
        Dmin7=median(Dmin_7loc);
        
        PT8{count2}=D8(i:i+nos_PT-1);
        Dmax_8loc(count2)=max(PT8{count2});
        Dmin_8loc(count2)=min(PT8{count2});
        Dmax8=median(Dmax_8loc);
        Dmin8=median(Dmin_8loc);
        
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
        
    elseif QRSl==6
        PT7{count2}=D7(i:i+nos_PT-1);
        Dmax_7loc(count2)=max(PT7{count2});
        Dmin_7loc(count2)=min(PT7{count2});
        Dmax7=median(Dmax_7loc);
        Dmin7=median(Dmin_7loc);
        
        PT8{count2}=D8(i:i+nos_PT-1);
        Dmax_8loc(count2)=max(PT8{count2});
        Dmin_8loc(count2)=min(PT8{count2});
        Dmax8=median(Dmax_8loc);
        Dmin8=median(Dmin_8loc);
        
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
        
    elseif QRSl==7
        PT8{count2}=D8(i:i+nos_PT-1);
        Dmax_8loc(count2)=max(PT8{count2});
        Dmin_8loc(count2)=min(PT8{count2});
        Dmax8=median(Dmax_8loc);
        Dmin8=median(Dmin_8loc);
        
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
        
    elseif QRSl==8
        PT9{count2}=D9(i:i+nos_PT-1);
        Dmax_9loc(count2)=max(PT9{count2});
        Dmin_9loc(count2)=min(PT9{count2});
        Dmax9=median(Dmax_9loc);
        Dmin9=median(Dmin_9loc);
        
    elseif QRS==9
        break
    end
end



for ifil=1:length(ecgdata)
    if QRSl==1%EM
        %level 2
        if D2(ifil)>bpt2*Dmin2 && D2(ifil)<bpt2*Dmax2
            D2(ifil)=0;
        else
            D2(ifil)=D2(ifil);
        end
        %level 3
        if D3(ifil)>bpt3*Dmin3 && D3(ifil)<bpt3*Dmax3
            D3(ifil)=0;
        else
            D3(ifil)=D3(ifil);
        end
        %level 4
        if D54(ifil)>bpt4*Dmin4 && D4(ifil)<bpt4*Dmax4
            D4(ifil)=0;
        else
            D4(ifil)=D4(ifil);
        end
        %level 5
        if D5(ifil)>bpt5*Dmin5 && D5(ifil)<bpt5*Dmax5
            D5(ifil)=0;
        else
            D5(ifil)=D5(ifil);
        end
        %level 6
        if D6(ifil)>bpt6*Dmin6 && D6(ifil)<bpt6*Dmax6
            D6(ifil)=0;
        else
            D6(ifil)=D6(ifil);
        end
        %level 7
        if D7(ifil)>bpt7*Dmin7 && D7(ifil)<bpt7*Dmax7
            D7(ifil)=0;
        else
            D7(i)=D7(i);
        end
        %level 8
        if D8(ifil)>bpt8*Dmin8 && D8(ifil)<bpt8*Dmax8
            D8(ifil)=0;
        else
            D8(ifil)=D8(ifil);
        end
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==2
        %level 3
        if D3(ifil)>bpt3*Dmin3 && D3(ifil)<bpt3*Dmax3
            D3(ifil)=0;
        else
            D3(ifil)=D3(ifil);
        end
        %level 4
        if D54(ifil)>bpt4*Dmin4 && D4(ifil)<bpt4*Dmax4
            D4(ifil)=0;
        else
            D4(ifil)=D4(ifil);
        end
        %level 5
        if D5(ifil)>bpt5*Dmin5 && D5(ifil)<bpt5*Dmax5
            D5(ifil)=0;
        else
            D5(ifil)=D5(ifil);
        end
        %level 6
        if D6(ifil)>bpt6*Dmin6 && D6(ifil)<bpt6*Dmax6
            D6(ifil)=0;
        else
            D6(ifil)=D6(ifil);
        end
        %level 7
        if D7(ifil)>bpt7*Dmin7 && D7(ifil)<bpt7*Dmax7
            D7(ifil)=0;
        else
            D7(i)=D7(i);
        end
        %level 8
        if D8(ifil)>bpt8*Dmin8 && D8(ifil)<bpt8*Dmax8
            D8(ifil)=0;
        else
            D8(ifil)=D8(ifil);
        end
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==3
        %level 4
        if D4(ifil)>bpt4*Dmin4 && D4(ifil)<bpt4*Dmax4
            D4(ifil)=0;
        else
            D4(ifil)=D4(ifil);
        end
        %level 5
        if D5(ifil)>bpt5*Dmin5 && D5(ifil)<bpt5*Dmax5
            D5(ifil)=0;
        else
            D5(ifil)=D5(ifil);
        end
        %level 6
        if D6(ifil)>bpt6*Dmin6 && D6(ifil)<bpt6*Dmax6
            D6(ifil)=0;
        else
            D6(ifil)=D6(ifil);
        end
        %level 7
        if D7(ifil)>bpt7*Dmin7 && D7(ifil)<bpt7*Dmax7
            D7(ifil)=0;
        else
            D7(i)=D7(i);
        end
        %level 8
        if D8(ifil)>bpt8*Dmin8 && D8(ifil)<bpt8*Dmax8
            D8(ifil)=0;
        else
            D8(ifil)=D8(ifil);
        end
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==4
        %level 5
        if D5(ifil)>bpt5*Dmin5 && D5(ifil)<bpt5*Dmax5
            D5(ifil)=0;
        else
            D5(ifil)=D5(ifil);
        end
        %level 6
        if D6(ifil)>bpt6*Dmin6 && D6(ifil)<bpt6*Dmax6
            D6(ifil)=0;
        else
            D6(ifil)=D6(ifil);
        end
        %level 7
        if D7(ifil)>bpt7*Dmin7 && D7(ifil)<bpt7*Dmax7
            D7(ifil)=0;
        else
            D7(i)=D7(i);
        end
        %level 8
        if D8(ifil)>bpt8*Dmin8 && D8(ifil)<bpt8*Dmax8
            D8(ifil)=0;
        else
            D8(ifil)=D8(ifil);
        end
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==5
        %level 6
        if D6(ifil)>bpt6*Dmin6 && D6(ifil)<bpt6*Dmax6
            D6(ifil)=0;
        else
            D6(ifil)=D6(ifil);
        end
        %level 7
        if D7(ifil)>bpt7*Dmin7 && D7(ifil)<bpt7*Dmax7
            D7(ifil)=0;
        else
            D7(i)=D7(i);
        end
        %level 8
        if D8(ifil)>bpt8*Dmin8 && D8(ifil)<bpt8*Dmax8
            D8(ifil)=0;
        else
            D8(ifil)=D8(ifil);
        end
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==6
        %level 7
        if D7(ifil)>bpt7*Dmin7 && D7(ifil)<bpt7*Dmax7
            D7(ifil)=0;
        else
            D7(i)=D7(i);
        end
        %level 8
        if D8(ifil)>bpt8*Dmin8 && D8(ifil)<bpt8*Dmax8
            D8(ifil)=0;
        else
            D8(ifil)=D8(ifil);
        end
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==7
        %level 8
        if D8(ifil)>bpt8*Dmin8 && D8(ifil)<bpt8*Dmax8
            D8(ifil)=0;
        else
            D8(ifil)=D8(ifil);
        end
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==8
        %level 9
        if D9(ifil)>bpt9*Dmin9 && D9(ifil)<bpt9*Dmax9
            D9(ifil)=0;
        else
            D9(ifil)=D9(ifil);
        end
    elseif QRSl==9
        break
    end
end


swd_n(1,:)=D1(1,1:NoDs);
swd_n(2,:)=D2(1,1:NoDs);
swd_n(3,:)=D3(1,1:NoDs);
swd_n(4,:)=D4(1,1:NoDs);
swd_n(5,:)=D5(1,1:NoDs);
swd_n(6,:)=D6(1,1:NoDs);
swd_n(7,:)=D7(1,1:NoDs);
swd_n(8,:)=D8(1,1:NoDs);
swd_n(9,:)=D9(1,1:NoDs);

%Reconstruction of Noise and denoise
Noise=iswt(swa,swd_n,wavename);
Denozsig1=ecgdata-Noise;

