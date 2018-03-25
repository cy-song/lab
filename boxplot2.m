

%{
%% boxplot
%teacher 
m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\928\teacer1\A\ADD_f1_A.mat');
Teacher_A_f1=m. ADDff1;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\teacher A\Teacher A f2.mat');
Teacher_A_f2=m. ADDff2;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\teacher A\Teacher A f3.mat');
Teacher_A_f3=m.ADDff3;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\teacher A\Teacher A f2_f1.mat');
Teacher_A_f2_f1=m.ADDff2_ff1;


%student 3 good 
m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s3 good\student03_good_f1.mat');
student03_good_f1=m. ADDff1;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s3 good\student03_good_f2.mat');
student03_good_f2=m. ADDff2;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s3 good\student03_good_f3.mat');
student03_good_f3=m. ADDff3;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s3 good\student03_good_f2-f1.mat');
student03_good_f2_f1=m. ADDff2_ff1;
%student 1 soso
m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\student  soso\student01_soso_f1.mat');
student01_soso_f1=m. ADDff1;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\student  soso\student01_soso_f2.mat');
student01_soso_f2=m. ADDff2;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\student  soso\student01_soso_f3.mat');
student01_soso_f3=m. ADDff3;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\student  soso\student01_soso_f2-f1.mat');
student01_soso_f2_f1=m. ADDff2_ff1;

%student 29 bad
m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s29 bad\student29_bad_f1.mat');
student29_bad_f1=m. ADDff1;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s29 bad\student29_bad_f2.mat');
student29_bad_f2=m. ADDff2;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s29 bad\student29_bad_f3.mat');
student29_bad_f3=m. ADDff3;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s29 bad\student29_bad_f2-f1.mat');
student29_bad_f2_f1=m. ADDff2_ff1;

%student 22 wrong
m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s22 wrong\student22_f1.mat');
student22_f1=m. ADDff1;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s22 wrong\student22_f2.mat');
student22_f2=m. ADDff2;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s22 wrong\student22_f3.mat');
student22_f3=m. ADDff3;

m=load('C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\920\overlap 0.8\s22 wrong\student22_f2-f1.mat');
student22_f2_f1=m. ADDff2_ff1;




%% combineboxplot
figure(1)
combineData = [Teacher_A_f1',student03_good_f1',student01_soso_f1',student29_bad_f1',student22_f1'];        % ?合
group = [ones(1,188)*5,ones(1,131)*1,ones(1,71)*0,ones(1,96)*(-1),ones(1,113)*(-5)];  % ?每一???的值，?定??，?里前20??0，后20??1
boxplot(combineData,group)
grid on;

%%
set(gca,'FontSize',16);
% ％進行二維(x,y)平面描點作圖，線條粗度2
% plot(ｘ,ｙ,'LineWidth',2)　　　　　　　　　　

%座標軸字體加粗大小14
set(gca,'FontWeight','bold','fontsize',14)

% ％繪圖範圍　x(1.52~1.58) y(0~1)
% axis([1.52,1.58,0,1]);　　　　　　　　　　　

% Create title,圖片標題，字體加粗大小16
title('Boxplot','FontWeight','bold','FontSize',16);

% Create xlabel,x座標名稱,字體加粗大小16
xlabel('不同分數比較','FontWeight','bold','FontSize',16);

% Create ylabel,y座標名稱,字體加粗大小16
ylabel('F1','FontWeight','bold','FontSize',16);



%%   F2

figure(2)
combineData = [Teacher_A_f2',student03_good_f2',student01_soso_f2',student29_bad_f2',student22_f2'];        % ?合
group = [ones(1,188)*5,ones(1,131)*1,ones(1,71)*0,ones(1,96)*(-1),ones(1,113)*(-5)];  % ?每一???的值，?定??，?里前20??0，后20??1
boxplot(combineData,group)
grid on;

%%
set(gca,'FontSize',16);
% ％進行二維(x,y)平面描點作圖，線條粗度2
% plot(ｘ,ｙ,'LineWidth',2)　　　　　　　　　　

%座標軸字體加粗大小14
set(gca,'FontWeight','bold','fontsize',14)

% ％繪圖範圍　x(1.52~1.58) y(0~1)
% axis([1.52,1.58,0,1]);　　　　　　　　　　　

% Create title,圖片標題，字體加粗大小16
title('Boxplot','FontWeight','bold','FontSize',16);

% Create xlabel,x座標名稱,字體加粗大小16
xlabel('不同分數比較','FontWeight','bold','FontSize',16);

% Create ylabel,y座標名稱,字體加粗大小16
ylabel('F2','FontWeight','bold','FontSize',16);


%% F3
figure(3)
combineData = [Teacher_A_f3',student03_good_f3',student01_soso_f3',student29_bad_f3',student22_f3'];        % ?合
group = [ones(1,188)*5,ones(1,131)*1,ones(1,71)*0,ones(1,96)*(-1),ones(1,113)*(-5)];  % ?每一???的值，?定??，?里前20??0，后20??1
boxplot(combineData,group)
grid on;

%%
set(gca,'FontSize',16);
% ％進行二維(x,y)平面描點作圖，線條粗度2
% plot(ｘ,ｙ,'LineWidth',2)　　　　　　　　　　

%座標軸字體加粗大小14
set(gca,'FontWeight','bold','fontsize',14)

% ％繪圖範圍　x(1.52~1.58) y(0~1)
% axis([1.52,1.58,0,1]);　　　　　　　　　　　

% Create title,圖片標題，字體加粗大小16
title('Boxplot','FontWeight','bold','FontSize',16);

% Create xlabel,x座標名稱,字體加粗大小16
xlabel('不同分數比較','FontWeight','bold','FontSize',16);

% Create ylabel,y座標名稱,字體加粗大小16
ylabel('F3','FontWeight','bold','FontSize',16);


%%  F2-F1

figure(4)
combineData = [Teacher_A_f2_f1',student03_good_f2_f1',student01_soso_f2_f1',student29_bad_f2_f1',student22_f2_f1'];        % ?合
group = [ones(1,188)*5,ones(1,131)*1,ones(1,71)*0,ones(1,96)*(-1),ones(1,113)*(-5)];  % ?每一???的值，?定??，?里前20??0，后20??1
boxplot(combineData,group)
grid on;

%%
set(gca,'FontSize',16);
% ％進行二維(x,y)平面描點作圖，線條粗度2
% plot(ｘ,ｙ,'LineWidth',2)　　　　　　　　　　

%座標軸字體加粗大小14
set(gca,'FontWeight','bold','fontsize',14)

% ％繪圖範圍　x(1.52~1.58) y(0~1)
% axis([1.52,1.58,0,1]);　　　　　　　　　　　

% Create title,圖片標題，字體加粗大小16
title('Boxplot','FontWeight','bold','FontSize',16);

% Create xlabel,x座標名稱,字體加粗大小16
xlabel('不同分數比較','FontWeight','bold','FontSize',16);

% Create ylabel,y座標名稱,字體加粗大小16
ylabel('F2-F1','FontWeight','bold','FontSize',16);

%% 平均標準差

Teacher_A_f1_mean=mean(Teacher_A_f1)
Teacher_A_f1_std=std(Teacher_A_f1);
Teacher_A_f1_var=var(Teacher_A_f1);
Teacher_A_f1_median=median(Teacher_A_f1);
Teacher_A_f1_mode=mode(Teacher_A_f1);
 skewness(Teacher_A_f1);
  kurtosis(Teacher_A_f1);



student03_good_f1_mean=mean(student03_good_f1)
student03_good_f1_std=std(student03_good_f1);
student03_good_f1_var=var(student03_good_f1);
student03_good_f1_median=median(student03_good_f1);
student03_good_f1_mode=mode(student03_good_f1);
 skewness(student03_good_f1);
  kurtosis(student03_good_f1);

student01_soso_f1_mean=mean(student01_soso_f1)
student01_soso_f1_std=std(student01_soso_f1);
student01_soso_f1_var=var(student01_soso_f1);
Cgood_median=median(student01_soso_f1);
Cgood_mode=mode(student01_soso_f1);
 skewness(student01_soso_f1);
  kurtosis(student01_soso_f1);

student29_bad_f1_mean=mean(student29_bad_f1)
student29_bad_f1_std=std(student29_bad_f1);
student29_bad_f1_var=var(student29_bad_f1);
student29_bad_f1_median=median(student29_bad_f1);
student29_bad_f1_mode=mode(student29_bad_f1);
 skewness(student29_bad_f1);
  kurtosis(student29_bad_f1);

  
student22_f1_mean=mean(student22_f1)
student22_f1_std=std(student22_f1);
student22_f1_var=var(student22_f1);
student22_f1_median=median(student22_f1);
student22_f1_mode=mode(student22_f1);
 skewness(student22_f1);
  kurtosis(student22_f1);
  %}

%{
%% cut the wav file
i=1;
[y,sr] = audioread( ['C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\519 老師重錄\老師重錄 drive-download-20170526T074516Z-001\t02/t02A0',num2str(i),'-01.wav']);
%[y,sr]  = audioread(['C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S22\s22A0',num2str(i),'-01.wav']);

%DIR = 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S03\';
%FILENAME = 't01A01-01.wav';

sound(y,sr)
y(:,2)=[];
   [staPT, endPT] = EPD(y, sr,1)
  y1=y(staPT:endPT);
  pause
  soundsc(y1,sr)
 % y2=y(9529:29801);
   % soundsc(y2,sr)
  
%    fprintf('press any button to save\n') ; pause   
%     wavefilea='a.wav';
%      nbits=8;
%     wavwrite(y1,sr,nbits,wavefilea);
   
   
 %} 
  
 %% 音長boxplotbox
 clear ;clc;clf;close all
 pause1=1;
 pause2=1;
 i=1;
 leteacher=[];
 % teacher 1
 [y,sr] = audioread( ['C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\519 老師重錄\老師重錄 drive-download-20170526T074516Z-001\t01/t01A0',num2str(i),'-01.wav']);
%y(:,2)=[];
soundsc(y,sr)
y(:,2)=[];
 [staPT, endPT] = EPD(y, sr,1);
%  di=endPT-staPT
  y1=y(staPT:endPT);
  teale1=(endPT-staPT)/sr
  leteacher=[leteacher,teale1];
  pause
  soundsc(y1,sr)
  % teacher 2
 [y,sr] = audioread( ['C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\519 老師重錄\老師重錄 drive-download-20170526T074516Z-001\t02/t02A0',num2str(i),'-01.wav']);
%y(:,2)=[];
soundsc(y,sr)
y(:,2)=[];
 [staPT, endPT] = EPD(y, sr,1);
%  di=endPT-staPT
  y1=y(staPT:endPT);
  teale2=(endPT-staPT)/sr
   leteacher=[leteacher,teale2];
  pause
  soundsc(y1,sr)
 
  
  
    % dictionary online
    for i1=1:4
 [y,sr] = audioread( ['C:\Users\jenny\Desktop\克 抗戰鐵路\106-1\paper 實驗 AHG\Vocabulary/1_',num2str(i1),'.wav']);
%y(:,2)=[];
 if pause2==1;
   pause
 end
soundsc(y,sr)
y(:,2)=[];

%normalization
y2ener=sqrt(sum(y.*y));
y=y./y2ener;

 [staPT, endPT] = EPD(y, sr,1);
%  di=endPT-staPT
  y1=y(staPT:endPT);
  teavocabula=(endPT-staPT)/sr;
   leteacher=[leteacher, teavocabula];
 if pause1==1;
   pause
 end
  soundsc(y1,sr)

    end
  
  
  
  
  
  
  
 
 le=[];
 for studnum=1:16

 if studnum>9
    opt1=1;    
else 
    opt1=0;
 end
 
%  q1=input('next student or not? 0 stay the same student')
%  if q1==0
%      studnum=studnum-1;
%  end
 studnum
  
    if opt1==0 %學生數字是1~9
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'A0',num2str(i),'-01.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'A',num2str(i),'-01.wav']);

        end
        
    else
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'A0',num2str(i),'-01.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'A',num2str(i),'-01.wav']);

        end
        
        
    end
    if pause2==1;
        pause
    end
soundsc(y2,sr)

%% student
y2(:,2)=[];

%  [staPT, endPT] = EPD(y2, sr,1);
%  pause
%normalization
y2ener=sqrt(sum(y2.*y2));
y2=y2./y2ener;

 [staPT, endPT] = EPD(y2, sr,1);
 di=(endPT-staPT)/sr
  y1=y2(staPT:endPT);
  le=[le,(endPT-staPT)/sr];
if pause1==1;
   pause
 end
  soundsc(y1,sr)
%pause
end
 
 
 for studnum=18:34

if studnum>9
    opt1=1;    
else 
    opt1=0;
end
    studnum
    if opt1==0 %學生數字是1~9
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'A0',num2str(i),'-01.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S0',num2str(studnum),'/S0',num2str(studnum),'A',num2str(i),'-01.wav']);

        end
        
    else
        if i<10
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'A0',num2str(i),'-01.wav']);
        else
            [y2,sr] =audioread([ 'C:\Users\jenny\Desktop\克 抗戰鐵路\105-2\paper 實驗 AHG\414 同學\handoff_0314-20170331T082045Z-001\handoff_0314\S',num2str(studnum),'/S',num2str(studnum),'A',num2str(i),'-01.wav']);

        end
        
        
    end
    if pause2==1;
        pause
    end
soundsc(y2,sr)

%% student
y2(:,2)=[];
 [staPT, endPT] = EPD(y2, sr,1);
 di=(endPT-staPT)/sr
  y1=y2(staPT:endPT);
  le=[le,(endPT-staPT)/sr];
if pause1==1;
   pause
 end
  soundsc(y1,sr)
%pause
 end
%}
 
legood=[le(3);le(6);le(7);le(8);le(10);le(12);le(13);le(14);le(26-1);le(27-1);le(30-1);le(34-1)];%exceed 17需要-1
lemid=[le(2);le(5);le(16);le(21-1);le(23-1)];
lebad=[le(29-1)];
leerror=[le(1);le(22-1)];
score=ones(12,1);
scatter(score,legood); hold on;
score1=zeros(5,1);
scatter(score1,lemid);
score2=[-1];
scatter(score2,lebad);
score5=[-5,-5];
scatter(score5,leerror);
score0=5*ones(6,1);%teacher 12
scatter(score0,leteacher);
grid on; axis([-6,6,-inf,inf]);
xlabel('老師評分','FontWeight','bold','FontSize',16);
ylabel('音長(sec)','FontWeight','bold','FontSize',16);
title('A ','FontWeight','bold','FontSize',16);

  %% combine
  
  figure(6)
combineData = [leteacher,legood',lemid',lebad',leerror'];        % ?合
group = [ones(1,6)*(5),ones(1,12),zeros(1,5),ones(1,2)*(-1),ones(1,1)*(-5)];  % ?每一???的值，?定??，?里前20??0，后20??1
boxplot(combineData,group)
grid on;
%%
set(gca,'FontSize',16);
% ％進行二維(x,y)平面描點作圖，線條粗度2
% plot(ｘ,ｙ,'LineWidth',2)　　　　　　　　　　

%座標軸字體加粗大小14
set(gca,'FontWeight','bold','fontsize',14)

% ％繪圖範圍　x(1.52~1.58) y(0~1)
% axis([1.52,1.58,0,1]);　　　　　　　　　　　

% Create title,圖片標題，字體加粗大小16
title('A  Boxplot','FontWeight','bold','FontSize',16);

% Create xlabel,x座標名稱,字體加粗大小16
xlabel('老師評分','FontWeight','bold','FontSize',16);

% Create ylabel,y座標名稱,字體加粗大小16
ylabel('音長(sec)','FontWeight','bold','FontSize',16);

  