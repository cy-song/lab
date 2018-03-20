clc ;clear;close all

i=3;
%% end point detection all
for times=1:5

     % {  
    [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\endpoint detection all/t1W0',num2str(i),'-0',num2str(times),'.wav']);
 
%}

%% student
y2(:,2)=[];
y1=y2;
% soundsc(y1,sr)

%% new EPD
  [staPT, endPT] = EPD(y2,sr,1)
  y1=y2(staPT:endPT);
%   pause
  soundsc(y1,sr)
time6=(endPT-staPT)/sr;
% studnum=studnum-1;
wavefilea=['C:\Users\jenny\Desktop\endpoint detection all\t1_cat_',num2str(times),'.wav'];%%%%%save 
 nbits=16;
 audiowrite(wavefilea,y1,sr);






pause
end
%}