clc; clear;close all

i=3;
%% teacher picture
%
opt.windowType = 'hann'; % rectwin, hann, or hamming etc.

window_length = 0.02; % Second
overlap = 0.5; % Ratio of overlap and window length


%%
window_gen = str2func( opt.windowType );
% for studnum=1:16
   
for studnum=18:18
 clf;clear sig rms_mag ROR figure(1) z1
 
if studnum>9
    opt1=1;    
else 
    opt1=0;
end
    studnum
    if opt1==0 %學生數字是1~9
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'W0',num2str(i),'-05.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'W',num2str(i),'-05.wav']);

        end
        
    else
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'W0',num2str(i),'-05.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'W',num2str(i),'-05.wav']);

        end
        
        
    end

    
%}
%% student
% y1(:,2)=[];
y2(:,2)=[];
soundsc(y2,sr)
x=y2;
emphasis = 0
if emphasis == 1,
    y_emph = filter([1 -0.98],1,x); 
                %[PARAM] -0.95 may be tuned anywhere from 0.9 to 0.99
else
    y_emph = x;
end
x=y_emph;


window = window_gen( round(window_length*sr) );
noverlap = round(length(window)*overlap);
nfft = power(2, ceil( log2(length(window)) ));
% [S,F,T] =spectrogram( x, window, noverlap, nfft, sr);
%spectrogram( y1, window, noverlap, nfft, sr);
 frame_length = length(window); %%　"M"
    frame_shift = frame_length - noverlap;
nframe = ceil( ( length(x) - frame_length ) / frame_shift );
% numFrames = floor((length(x)-frame_length)/frame_shift)+1;

addmag=zeros(nframe,1);
for i1 = 1:nframe
    
     sig (:,1)=x([(1+( frame_length - noverlap)*(i1-1)):(frame_length+( frame_length - noverlap)*(i1-1))], 1);
    
rms_mag(i1)=10*log10( 1000* sqrt(sum(sig (:,1).^2)/ frame_length) );
end
figure(1)
subplot(4,1,1)
tt=length(x)/sr;
t=linspace(0,tt,length(x));

plot(t,x)
 ylabel('waveform amplitude')
grid on;
subplot(4,1,2)
plot(rms_mag,'-o')
 ylabel('RMS mag')
 
 delta_t=length(frame_length)/sr;
for i2=2:nframe
ROR(i2) = (rms_mag(i2)-rms_mag(i2-1))/ delta_t;
end
subplot(4,1,3)
 plot(ROR)
 ylabel('ROR')
%%%%%%%%%%%%%%%%%%%% 
for i3= 1:nframe
z1(i3)=ZCR(x([(1+( frame_length - noverlap)*(i3-1)):(frame_length+( frame_length - noverlap)*(i3-1))]));
end
z2=ZCR(x);
 subplot(4,1,4)
 plot(z1)
ylabel('ZCR')

saveas(figure(1) , ['C:\Users\jenny\Desktop\克 抗戰鐵路\106-2\paper 實驗 AHG\0306\cat',num2str(studnum),'_5'] , 'bmp' )
pause
end

