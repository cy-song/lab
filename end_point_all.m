clc ;clear;close all

i=1;
%% end point detection all
for studnum=18:34

if studnum>9
    opt1=1;    
else 
    opt1=0;
end
    studnum
    if opt1==0 %�ǥͼƦr�O1~9
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\�J �ܾ��K��\105-2\paper ���� AHG\414 �P��\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'A0',num2str(i),'-01.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\�J �ܾ��K��\105-2\paper ���� AHG\414 �P��\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'A',num2str(i),'-01.wav']);

        end
        
    else
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\�J �ܾ��K��\105-2\paper ���� AHG\414 �P��\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'A0',num2str(i),'-01.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\�J �ܾ��K��\105-2\paper ���� AHG\414 �P��\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'A',num2str(i),'-01.wav']);

        end
        
        
    end
%soundsc(y2,sr)

%% student
y2(:,2)=[];
y1=y2;
soundsc(y1,sr)

%% new EPD
  [staPT, endPT] = EPD(y2,sr,1)
  y1=y2(staPT:endPT);
%   pause
  soundsc(y1,sr)
time6=(endPT-staPT)/sr;
studnum=studnum-1;
wavefilea=['C:\Users\jenny\Desktop\endpoint detection all\a',num2str(studnum),'.wav'];%%%%%save 
 nbits=16;
 audiowrite(wavefilea,y1,sr);






pause
end
%}