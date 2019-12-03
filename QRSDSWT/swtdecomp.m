function [swd,swa,swd_lp,ecgdata,NoDs]=swtdecomp(sig1,level,wavename,fs,gr,NoD)
%gr  - plot switch
%NoD - Number of Data extrated from input data to ecgdata

%Check swt size
cl=2^level;%used to check swt
NoDc=NoD/cl;
NoDs=round(NoDc);
NoDs=NoDs*cl;%Required Number of Data for SWT

%Data Initialization
ecgdata=zeros(1,NoDs);
t=0:1/fs:(NoDs/fs);
if NoDs<NoD
    ecgdata(1:NoDs)=sig1(1:NoDs);
else
    ecgdata(1:NoD)=sig1(1:NoD);
end
% level=9;
% wavename='haar';
LenECG=length(ecgdata);
% t = 0:1/fs:(LenECG-1)/fs;
t = t(1:length(t)-1);

%SWT Decomposed to 9 levels
[swa,swd]=swt(ecgdata,level,wavename);
D1=swd(1,:);
D2=swd(2,:);
D3=swd(3,:);
D4=swd(4,:);
D5=swd(5,:);
D6=swd(6,:);
D7=swd(7,:);
D8=swd(8,:);
D9=swd(9,:);


%from the following figure one can get QRS is about 0.1s
%From detected R peak back 39 and forward 13 points cover the entire
%QRS segment
D1_lp=lowpass(D1,2,fs);
D2_lp=lowpass(D2,2,fs);
D3_lp=lowpass(D3,2,fs);
D4_lp=lowpass(D4,2,fs);
D5_lp=lowpass(D5,2,fs);
D6_lp=lowpass(D6,2,fs);
D7_lp=lowpass(D7,2,fs);
D8_lp=lowpass(D8,2,fs);
D9_lp=lowpass(D9,2,fs);

if gr==1
    figure('Name','SWT Decompositions')
    subplot(9,1,1)
    plot(t,D1)
    title('Detailed Coef D1','position',[-1,0])
    subplot(9,1,2)
    plot(t,D2)
    title('Detailed Coef D2','position',[-1,0])
    subplot(9,1,3)
    plot(t,D3)
    title('Detailed Coef D3','position',[-1,0])
    subplot(9,1,4)
    plot(t,D4)
    title('Detailed Coef D4','position',[-1,0])
    subplot(9,1,5)
    plot(t,D5)
    title('Detailed Coef D5','position',[-1,0])
    subplot(9,1,6)
    plot(t,D6)
    title('Detailed Coef D6','position',[-1,0])
    subplot(9,1,7)
    plot(t,D7)
    title('Detailed Coef D7','position',[-1,0])
    subplot(9,1,8)
    plot(t,D8)
    title('Detailed Coef D8','position',[-1,0])
    subplot(9,1,9)
    plot(t,D9)
    title('Detailed Coef D9','position',[-1,0])
end


if gr==1
    figure('Name','SWT Decompositions After Lowpass')
    subplot(9,1,1)
    plot(t,D1_lp)
    title('Detailed Coef D1','position',[-1,0])
    subplot(9,1,2)
    plot(t,D2_lp)
    title('Detailed Coef D2','position',[-1,0])
    subplot(9,1,3)
    plot(t,D3_lp)
    title('Detailed Coef D3','position',[-1,0])
    subplot(9,1,4)
    plot(t,D4_lp)
    title('Detailed Coef D4','position',[-1,0])
    subplot(9,1,5)
    plot(t,D5_lp)
    title('Detailed Coef D5','position',[-1,0])
    subplot(9,1,6)
    plot(t,D6_lp)
    title('Detailed Coef D6','position',[-1,0])
    subplot(9,1,7)
    plot(t,D7_lp)
    title('Detailed Coef D7','position',[-1,0])
    subplot(9,1,8)
    plot(t,D8_lp)
    title('Detailed Coef D8','position',[-1,0])
    subplot(9,1,9)
    plot(t,D9_lp)
    title('Detailed Coef D9','position',[-1,0])
end

swd_lp(1,:)=D1_lp(1,:);
swd_lp(2,:)=D2_lp(1,:);
swd_lp(3,:)=D3_lp(1,:);
swd_lp(4,:)=D4_lp(1,:);
swd_lp(5,:)=D5_lp(1,:);
swd_lp(6,:)=D6_lp(1,:);
swd_lp(7,:)=D7_lp(1,:);
swd_lp(8,:)=D8_lp(1,:);
swd_lp(9,:)=D9_lp(1,:);




