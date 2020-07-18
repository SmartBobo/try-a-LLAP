clc;
clear;
TEMPERATURE = 30;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;

%prepare the microphone
Microphone1 = audioDeviceReader;
Microphone1.SampleRate = 48000;
Microphone1.SamplesPerFrame = 1920;
Microphone1.Device ='麦克风 (Realtek(R) Audio)';
setup(Microphone1);


%create streaming sine wave and cos wave object
numofFre = 16;
cosWave = cell(numofFre,1);
sinWave = cell(numofFre,1);
delay = 0;
for i = 1:numofFre
    cosWave{i,1} =dsp.SineWave(1,(17150+i*350),pi/2-delay*2*pi*(17150+i*350),'SampleRate',48000,'SamplesPerFrame',1920);
    sinWave{i,1} =dsp.SineWave(1,(17150+i*350),pi-delay*2*pi*(17150+i*350),'SampleRate',48000,'SamplesPerFrame',1920);
end

%set CIC Filter
cicDecim = dsp.CICDecimator(16,17,3);

%set variables
framelength = 1920;
deciSize = framelength/16;
DCvalueI1 = zeros(numofFre,1);
DCvalueR1 = zeros(numofFre,1);
maxR1 = zeros(numofFre,1);
maxI1 = zeros(numofFre,1);
minR1 = zeros(numofFre,1);
minI1 = zeros(numofFre,1);
totDis1 = 400;
I1 = zeros(deciSize,numofFre); 
Q1 = zeros(deciSize,numofFre);
DCvalueI2 = zeros(numofFre,1);
DCvalueR2 = zeros(numofFre,1);
maxR2 = zeros(numofFre,1);
maxI2 = zeros(numofFre,1); 
minR2 = zeros(numofFre,1);
minI2 = zeros(numofFre,1);
totDis2 = 400;
I2 = zeros(deciSize,numofFre); 
Q2 = zeros(deciSize,numofFre);
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,8);
temp2 = zeros(deciSize,8);
data1 = zeros(750,1);
data2 = zeros(750,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
numCount1 = 0;
tempInPhase1 = zeros(framelength*2/5,numofFre);
tempQuard1 = zeros(framelength*2/5,numofFre);
tempInPhase2 = zeros(framelength,numofFre);
tempQuard2 = zeros(framelength,numofFre);
variance1 = zeros(numofFre,1);
variance2 = zeros(numofFre,1);
ampThr = [33000 40000 35000 30000 46000 42000 40000 40000 12000 13000 10000 6000 5000 4000 4000 4000];
L1 = 150;
L2 = 150;

%prepare the audio
dt = 1/48000;
t = 0:dt:35-dt;
fc1 = 17150;
A = 0.2;
y1 = 0; %A*cos(2*pi*fc*t);
for i = 1:8
    y1 = y1 + A*cos(2*pi*(fc1+350*i)*t);
end
y2 = 0; %A*cos(2*pi*fc*t);
fc2 = 19950;
for i = 1:8
    y2 = y2 + A*cos(2*pi*(fc2+350*i)*t);
end
y1 = y1/8;
y2 = y2/8;
player1 = audioplayer(y1,48000,16,5);
player2 = audioplayer(y2,48000,16,4);
play(player1);
play(player2);

% %initiate the figure
% axes1 = subplot(2,1,1);
% pic1 = plot(axes1,0:1/25:30-1/25,data1);
% set(axes1,'xlim', [0 30], 'ylim', [0 800]);
% drawnow;
% axes2 = subplot(2,1,2);
% pic2 = plot(axes2,0:1/25:30-1/25,data2);
% set(axes2,'xlim', [0 30], 'ylim', [0 800]);
% drawnow;

%initiate the figure
axes1 = subplot(1,1,1);
pic = plot(axes1,data1,data2);
set(axes1,'xlim', [0 600], 'ylim', [-150 150]);
drawnow;
pause(5);
disp('Start Recording');
i = 0;
tic;
%可能还要考虑Latency due to sound card's output buffer: 0.046440 seconds
while toc <= 30
   [audioIn1,Overrun1] = step(Microphone1);        % 采样

   if Overrun1 > 0
      warning('  Microphone1 数据溢出 %d 位\n',Overrun1);
   end

   sumxy1 = 0;
   sumy1 = 0;
   numCount1 = 0;
   sumxy2 = 0;
   sumy2 = 0;
   numCount2 = 0;
   distance = 0;
   i = i + 1;
   
%  Calculate the distance of Microphone1
   for freNum = 1:16    
        
        
        %将原始信号转换成Baseband
        [I1,Q1,tempInPhase1(:,freNum),tempQuard1(:,freNum)] = rawToBasebandver2(audioIn1,sinWave{freNum,1}(),cosWave{freNum,1}(),tempInPhase1(:,freNum),tempQuard1(:,freNum),cicDecim);
        
%         if freNum < 9
            %将两路信号用LEVD算法处理并求出DC vector
            [ESCInPhase,ESCQuard,DCvalueR1(freNum),DCvalueI1(freNum),maxR1(freNum),minR1(freNum),maxI1(freNum),minI1(freNum)] = LEVD2(I1,Q1,DCvalueR1(freNum),DCvalueI1(freNum),maxR1(freNum),minR1(freNum),maxI1(freNum),minI1(freNum),ampThr(freNum));
%         else
%             [ESCInPhase,ESCQuard,DCvalueR1(freNum),DCvalueI1(freNum),maxR1(freNum),minR1(freNum),maxI1(freNum),minI1(freNum)] = LEVD1(I1,Q1,DCvalueR1(freNum),DCvalueI1(freNum),maxR1(freNum),minR1(freNum),maxI1(freNum),minI1(freNum));
%         end
            
        %除去信号中的static vector 并将两路信号转变成baseband 
        DVInPhase = I1 - ESCInPhase; 
        DVQuard = Q1 - ESCQuard;
        DVBaseband(:,freNum) = ReImToComp(DVInPhase, DVQuard);
        % 求出相位同时根据threshold继续对DC vector进行处理

        [ph,DCvalueR1(freNum),DCvalueI1(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR1(freNum),minR1(freNum),DCvalueR1(freNum),maxI1(freNum),minI1(freNum),DCvalueI1(freNum),freCount(freNum));

        fre = fc1 + freNum * 350;


        %计算左音响信号方差
        if freCount(freNum) == 1 && freNum < 9
            for a = 1:deciSize
                temp(a,freNum) = ph(a) - ph(1); 
                temp(a,freNum) = temp(a,freNum)*SoundSpeed*1000/(2*pi*fre);
            end
            sumy1 = sum(temp(:,freNum))+sumy1;
            sumxy1 = sum(temp(:,freNum).*t1)+sumxy1;
            numCount1 = numCount1 + 1;
        end
        %计算右音响信号方差
        if freCount(freNum) == 1 && freNum > 8
            for a = 1:deciSize
                temp(a,freNum) = ph(a) - ph(1); 
                temp(a,freNum) = temp(a,freNum)*SoundSpeed*1000/(2*pi*fre);
            end
            sumy2 = sum(temp(:,freNum))+sumy2;
            sumxy2 = sum(temp(:,freNum).*t1)+sumxy2;
            numCount2 = numCount2 + 1;
        end
    end
    
    %find ignore frequency
    if numCount1 == 0 && numCount2 == 0
       data1(i) = data1(i-1);
       data2(i) = data2(i-1);
       continue;
    end

    
    %least square linear regression for speaker1 
    deltax = 8*((deciSize-1)*deciSize*(2*deciSize-1)/6-(deciSize-1)*deciSize*(deciSize-1)/4);
    delta1 = (sumxy1-sumy1*(deciSize-1)/2)/deltax*8/numCount1;
    
    tempvar1 = delta1*t1;
    varsum = 0;
    
    %get variance of each frequency
    for freNum = 1:8
        if freCount(freNum) == 1
            tempvar2 = tempvar1 - temp(:,freNum);
            variance1(freNum) = sum(tempvar2.^2);
            varsum = varsum + variance1(freNum);
        end
    end
    
    for freNum = 1:8
        if variance1(freNum) > varsum/numCount1
            freCount(freNum) = 0;    
        end
    end
    
    sumxy1 = 0;
    sumy1 = 0;
    numCount1 = 0;
    
    %Caluculate the distance based on the frequency after removing the ignore frequency
    for freNum = 1:8
       if freCount(freNum) == 1
           sumy1 = sum(temp(:,freNum))+sumy1;
           sumxy1 = sum(temp(:,freNum).*t1)+sumxy1;
           numCount1 = numCount1 + 1;
       end

    end
    if numCount1 ~= 0
        delta1 = (sumxy1-sumy1*(deciSize-1)/2)/deltax*8/numCount1;
        distance = -delta1*deciSize;
        totDis1 = totDis1 + distance;
    end
    
   %least square linear regression for speaker2
    deltax = 8*((deciSize-1)*deciSize*(2*deciSize-1)/6-(deciSize-1)*deciSize*(deciSize-1)/4);
    delta2 = (sumxy2-sumy2*(deciSize-1)/2)/deltax*8/numCount2;
    
    tempvar1 = delta2*t1;
    varsum = 0;
    
    %get variance of each frequency
    for freNum = 9:16
        if freCount(freNum) == 1
            tempvar2 = tempvar1 - temp(:,freNum);
            variance1(freNum) = sum(tempvar2.^2);
            varsum = varsum + variance1(freNum);
        end
    end
    
    for freNum = 9:16
        if variance1(freNum) > varsum/numCount2
            freCount(freNum) = 0;    
        end
    end
    
    sumxy2 = 0;
    sumy2 = 0;
    numCount2 = 0;
    
    %Caluculate the distance based on the frequency after removing the ignore frequency
    for freNum = 9:16
       if freCount(freNum) == 1
           sumy2 = sum(temp(:,freNum))+sumy2;
           sumxy2 = sum(temp(:,freNum).*t1)+sumxy2;
           numCount2 = numCount2 + 1;
       end

    end
    if numCount2 ~= 0
        delta2 = (sumxy2-sumy2*(deciSize-1)/2)/deltax*8/numCount2;
        distance = -delta2*deciSize;
        totDis2 = totDis2 + distance;
    end
    
 %  calculate the x and y
    data1(i) = sqrt((totDis1^2 - L1^2)*(totDis2^2 - L2^2)*((L1+L2)^2 - (totDis1 - totDis2)^2))/(2*totDis1*L2 + 2*totDis2*L1);
    data2(i) = (totDis2*L1^2 - totDis1*L2^2 - totDis1^2*totDis2+totDis2^2*totDis1)/(2*(totDis1*L2 + totDis2*L1));
%     data1(i) = totDis1;
%     data2(i) = totDis2;
    
    set(pic,'xdata',data1,'ydata',data2);
%     set(pic1,'ydata',data1);                                                                                        
%     drawnow;
%     set(pic2,'ydata',data2);
    drawnow;
    
            
end

