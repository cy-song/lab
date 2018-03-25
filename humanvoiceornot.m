%% determine whether there is human voice or not.
%man 85~180 hz  woman 165~255   tone wu set 50~500 
clc; clear ;clf ;
DIR = 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\519 老師重錄\老師重錄 drive-download-20170526T074516Z-001\t01/';
FILENAME = 't01A01-01.wav';


[y,fs1] = audioread([DIR FILENAME]);

soundsc(y,fs1);
y(:,2)=[];
%% new EPD
  [staPT, endPT] = EPD(y,fs1,1)
  y1=y(staPT:endPT);
  pause
  soundsc(y1,fs1)
time6=(endPT-staPT)/fs1;
wavefilea='a.wav';%%%%%save 
 nbits=16;
%  wavwrite(y1,fs1,nbits,wavefilea);
 
%  
% fs = 16000;
% 
% y = resample(y1,fs,fs1);
%% Parameters to play with
framelen = 0.032; % second. [INVESTIGATE]
p = 25; % linear prediction order. [INVESTIGATE]
%p愈少  峰值愈明顯
%%
L = round(framelen*fs1);
overlapratio=0.5;
frame_shift=round(L*(1-overlapratio));
noverlap = round(length(L)*overlapratio);

if L<=p
    disp('Linear prediction requires the num of equations to be greater than the number of variables.');
end

sw.emphasis = 1; % default = 1

numFrames = floor((length(y)-L)/frame_shift)+1;
%ceil((L-winlengh)/shift) 
excitat = zeros(size(y));
e_n = zeros(p+L,1);

LPcoeffs = zeros(p+1,numFrames);
Kcoeffs = zeros(p,numFrames); % reflection coeffs

nfft = power(2, ceil( log2(length(L)) ));

Nfreqs = 2^nextpow2(2*L-1)/2; % Num points for plotting the inverse filter response
df = fs1/2/Nfreqs;
ff = 0:df:fs1/2-df;

if sw.emphasis == 1,
    y_emph = filter([1 -0.95],1,y); 
                %[PARAM] -0.95 may be tuned anywhere from 0.9 to 0.99
else
    y_emph = y;
end
%
%% Linear prediction and estimation of the source e_n
win = hamming(L);%ones(L,1);% hamming(L); % Rectangular window.
ADDff1=[];ADDff2_ff1=[];
ADDff2=[];ADDff3=[];
writerobj = VideoWriter('out.avi')
writerobj.FrameRate =10 ;%numFrames/time6;
open(writerobj);

figg=2;

for kk = 32
    ind = (kk-1)*frame_shift+1:(kk-1)*frame_shift+L;
    ywin = y_emph(ind).*win;
    %A = lpc(ywin,p); %% This is actually the more direct way to obtain the
    % LP coefficients. But, in this script, we used levinson() instead
    % because it gives us the "reflection coefficients".
    %
    % (Below, the fast way to calculate R is
    % copied and modified from MATLAB's lpc() function)
    Y = fft(ywin,2*Nfreqs);
    
    
    
    R = ifft(Y.*conj(Y));
    [A,errvar,K] = levinson(R,p);
    %% [INVESTIGATE] You are encouraged to try implementing the
    % Levinson-Durbin Algorithm, and check your results against levinson().
    % ...
    
    %% Insert your code to calculate formant frequencies from array A, and
    % produce a scattering plot on the f1 vs.(f2-f1) plane.
    % Hint: check subplot(223) and use freqz() to find out the the peak
    % frequencies.
    
    
    
    
    % ...
    
    %% Preparation for data visualization
    if kk == 1,
        e_n(p+1:end) = filter(A,[1],ywin);
    else
        ywin_extended = y((kk-1)*frame_shift+1-p:(kk-1)*frame_shift+L);  %% WORTH TWEAKING
        e_n = filter(A,[1],ywin_extended);
    end
    excitat(ind) = e_n(p+1:end);
    %
    figure(1);
    subplot(211);
    plot(ind/fs1*1000, y(ind));
    xlabel('ms')
    ylabel('orginal signal')
    % set(gca,'xlim',[kk-1 kk]*framelen*1000);
%     subplot(222);
%     plot(ind/fs*1000, e_n(p+1:end));
%     %  set(gca,'xlim',[kk-1 kk]*framelen*1000);
%     xlabel('ms')
%     ylabel('e(z)')
    
    
    subplot(212);
    [H,W] = freqz(1,A,Nfreqs);
    Hmag = 20*log10(abs(H));
    Ymag = 20*log10(abs(Y(1:Nfreqs))); % NEED DEBUG
    Hmax = max(Hmag);
    offset = max(Hmag) - max(Ymag);
%     plot(ff,Hmag);
%     hold on;
    plot(ff,Ymag,'r'); hold off;
%     set(gca,'xlim',[0 fs1/2],'ylim',[Hmax-50, Hmax+5]);
    xlabel('Hz')
    ylabel('orginal signal fourier')
    
    
    
    
%     subplot(224);
%     plot(ff,Ymag+offset-Hmag);
%      set(gca,'xlim',[0 fs/2],'ylim',[-30, 25]);
%     xlabel('Hz');
%     ylabel('e(z) fourier (e(z))')
    drawnow;
    %pause;
    LPcoeffs(:,kk) = A;
    Kcoeffs(:,kk) = K;
    
     if figg==1;
 frame = getframe(gcf);
    writeVideo(writerobj,frame)
    end
    %  K111=1:size(Kcoeffs,2);
    %  scatter(K111,Kcoeffs(:,1))
    %}
    N=size(Hmag(:),1);
    %%  
    
    figure(2)
    % [S,F,T] =spectrogram( x, window, noverlap, [], sr,'yaxis');
spectrogram( y(ind), win, noverlap, 10000, fs1,'yaxis');
% [S,F,T] = spectrogram( y(ind), win, noverlap,nfft, fs1);
% %% Show Spectrogram
% 
% surf(T,F,10*log10(abs(S)),'EdgeColor','none');
% axis xy; axis tight; view([0,90]);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
%  axis([-inf inf 0,5000]);
  if figg==2;
 frame = getframe(gcf);
    writeVideo(writerobj,frame)
    end
end
close(writerobj);











