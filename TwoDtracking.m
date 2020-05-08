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

Microphone2 = audioDeviceReader;
Microphone2.SampleRate = 48000;
Microphone2.SamplesPerFrame = 1920;
Microphone2.Device ='麦克风 (2- USB Audio Device)';
setup(Microphone2);

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
% cicDecim = dsp.CICDecimator(16,17,3);

%set variables
framelength = 1920;
deciSize = framelength/16;
DCvalueI1 = zeros(numofFre,1);
DCvalueR1 = zeros(numofFre,1);
maxR1 = zeros(numofFre,1);
maxI1 = zeros(numofFre,1);
minR1 = zeros(numofFre,1);
minI1 = zeros(numofFre,1);
totDis1 = 200;
I1 = zeros(deciSize,numofFre); 
Q1 = zeros(deciSize,numofFre);
DCvalueI2 = zeros(numofFre,1);
DCvalueR2 = zeros(numofFre,1);
maxR2 = zeros(numofFre,1);
maxI2 = zeros(numofFre,1);
minR2 = zeros(numofFre,1);
minI2 = zeros(numofFre,1);
totDis2 = 200;
I2 = zeros(deciSize,numofFre); 
Q2 = zeros(deciSize,numofFre);
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
data1 = zeros(750,1);
data2 = zeros(750,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount1 = zeros(numofFre,1);
freCount2 = zeros(numofFre,1);
numCount1 = 0;
tempInPhase1 = zeros(framelength,numofFre);
tempQuard1 = zeros(framelength,numofFre);
tempInPhase2 = zeros(framelength,numofFre);
tempQuard2 = zeros(framelength,numofFre);
variance1 = zeros(numofFre,1);
variance2 = zeros(numofFre,1);
L1 = 65;
L2 = 65;

%prepare the audio
dt = 1/48000;
t = 0:dt:35-dt;
fc = 17150;
A = 0.2;
y = 0; %A*cos(2*pi*fc*t);
for i = 1:numofFre
    y = y + A*cos(2*pi*(fc+350*i)*t);
end
y = y/numofFre;
sound(y,48000);

%initiate the figure
axes1 = subplot(1,1,1);
pic = plot(axes1,data1,data2);
set(axes1,'xlim', [150 300], 'ylim', [-150 100]);
tic;
drawnow;
pause(5);
disp('Start Recording');
i = 0;

%可能还要考虑Latency due to sound card's output buffer: 0.046440 seconds
while toc <= 35
   [audioIn1,Overrun1] = step(Microphone1);        % 采样
   [audioIn2,Overrun2] = step(Microphone2);
   if Overrun1 > 0
      warning('  Microphone1 数据溢出 %d 位\n',Overrun1);
   end
   if Overrun2 > 0
      warning('  Microphone2 数据溢出 %d 位\n',Overrun2);
   end
   sumxy = 0;
   sumy = 0;
   numCount = 0;
   distance = 0;
   i = i + 1;
   
%  Calculate the distance of Microphone1
   for freNum = 1:4    
        
       
        %将原始信号转换成Baseband
        [I1,Q1,tempInPhase1(:,freNum),tempQuard1(:,freNum)] = rawToBaseband(audioIn1,sinWave{freNum,1}(),cosWave{freNum,1}(),tempInPhase1(:,freNum),tempQuard1(:,freNum),deciSize);
        
        
        %将两路信号用LEVD算法处理并求出DC vector
        [ESCInPhase,ESCQuard,DCvalueR1(freNum),DCvalueI1(freNum),maxR1(freNum),minR1(freNum),maxI1(freNum),minI1(freNum)] = LEVD2(I1,Q1,DCvalueR1(freNum),DCvalueI1(freNum),maxR1(freNum),minR1(freNum),maxI1(freNum),minI1(freNum));


        %除去信号中的static vector 并将两路信号转变成baseband 
        DVInPhase = I1 - ESCInPhase; 
        DVQuard = Q1 - ESCQuard;
        DVBaseband(:,freNum) = ReImToComp(DVInPhase, DVQuard);
        % 求出相位同时根据threshold继续对DC vector进行处理

        [ph,DCvalueR1(freNum),DCvalueI1(freNum),freCount1(freNum)] = DCprocess(DVBaseband(:,freNum),maxR1(freNum),minR1(freNum),DCvalueR1(freNum),maxI1(freNum),minI1(freNum),DCvalueI1(freNum),freCount1(freNum));

        fre = fc + freNum * 350;


        %计算方差
        if freCount1(freNum) == 1
            for a = 1:deciSize
                temp(a,freNum) = ph(a) - ph(1); 
                temp(a,freNum) = temp(a,freNum)*SoundSpeed*1000/(2*pi*fre);
            end
            sumy = sum(temp(:,freNum))+sumy;
            sumxy = sum(temp(:,freNum).*t1)+sumxy;
            numCount = numCount + 1;
        end
    end
    
    %find ignore frequency
    if numCount == 0
       data1(i) = totDis1;
       data2(i) = totDis2;
        continue;
    end
    
    %least square linear regression
    deltax = numofFre*((deciSize-1)*deciSize*(2*deciSize-1)/6-(deciSize-1)*deciSize*(deciSize-1)/4);
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    
    temp1 = delta*t1;
    varsum = 0;
    
    %get variance of each frequency
    for freNum = 1:numofFre
        if freCount1(freNum) == 1
            temp2 = temp1 - temp(:,freNum);
            variance1(freNum) = sum(temp2.^2);
            varsum = varsum + variance1(freNum);
        end
    end
    
    for freNum = 1:numofFre
        if variance1(freNum) > varsum/numCount
            freCount1(freNum) = 0;    
        end
    end
    
    sumxy = 0;
    sumy = 0;
    numCount = 0;
    
    %Caluculate the distance based on the frequency after removing the ignore frequency
    for freNum = 1:numofFre
       if freCount1(freNum) == 1
           sumy = sum(temp(:,freNum))+sumy;
           sumxy = sum(temp(:,freNum).*t1)+sumxy;
           numCount = numCount + 1;
       end

    end
    if numCount == 0
       data1(i) = totDis1;
       data2(i) = totDis2;
       continue;
    end
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    distance = -delta*deciSize/2;
    totDis1 = totDis1 + distance;
    
    
    %  Calculate the distance of Microphone2
    sumxy = 0;
    sumy = 0;
    numCount = 0;
    distance = 0;
    for freNum = 1:4    
        
       
        %将原始信号转换成Baseband
        [I2,Q2,tempInPhase2(:,freNum),tempQuard2(:,freNum)] = rawToBaseband(audioIn2,sinWave{freNum,1}(),cosWave{freNum,1}(),tempInPhase2(:,freNum),tempQuard2(:,freNum),deciSize);
        
        
        %将两路信号用LEVD算法处理并求出DC vector
        [ESCInPhase,ESCQuard,DCvalueR2(freNum),DCvalueI2(freNum),maxR2(freNum),minR2(freNum),maxI2(freNum),minI2(freNum)] = LEVD2(I2,Q2,DCvalueR2(freNum),DCvalueI2(freNum),maxR2(freNum),minR2(freNum),maxI2(freNum),minI2(freNum));


        %除去信号中的static vector 并将两路信号转变成baseband 
        DVInPhase = I2 - ESCInPhase; 
        DVQuard = Q2 - ESCQuard;
        DVBaseband(:,freNum) = ReImToComp(DVInPhase, DVQuard);
        % 求出相位同时根据threshold继续对DC vector进行处理

        [ph,DCvalueR2(freNum),DCvalueI2(freNum),freCount2(freNum)] = DCprocess(DVBaseband(:,freNum),maxR2(freNum),minR2(freNum),DCvalueR2(freNum),maxI2(freNum),minI2(freNum),DCvalueI2(freNum),freCount2(freNum));

        fre = fc + freNum * 350;


        %计算方差
        if freCount2(freNum) == 1
            for a = 1:deciSize
                temp(a,freNum) = ph(a) - ph(1); 
                temp(a,freNum) = temp(a,freNum)*SoundSpeed*1000/(2*pi*fre);
            end
            sumy = sum(temp(:,freNum))+sumy;
            sumxy = sum(temp(:,freNum).*t1)+sumxy;
            numCount = numCount + 1;
        end
    end
    
    %find ignore frequency
    if numCount == 0
       data1(i) = totDis1;
       data2(i) = totDis2;
       continue;
    end
    
    %least square linear regression
    deltax = numofFre*((deciSize-1)*deciSize*(2*deciSize-1)/6-(deciSize-1)*deciSize*(deciSize-1)/4);
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    
    temp1 = delta*t1;
    varsum = 0;
    
    %get variance of each frequency
    for freNum = 1:numofFre
        if freCount2(freNum) == 1
            temp2 = temp1 - temp(:,freNum);
            variance2(freNum) = sum(temp2.^2);
            varsum = varsum + variance2(freNum);
        end
    end
    
    for freNum = 1:numofFre
        if variance2(freNum) > varsum/numCount
            freCount2(freNum) = 0;    
        end
    end
    
    sumxy = 0;
    sumy = 0;
    numCount = 0;
    
    %Caluculate the distance based on the frequency after removing the ignore frequency
    for freNum = 1:numofFre
       if freCount2(freNum) == 1
           sumy = sum(temp(:,freNum))+sumy;
           sumxy = sum(temp(:,freNum).*t1)+sumxy;
           numCount = numCount + 1;
       end

    end
    if numCount == 0
       data1(i) = totDis1;
       data2(i) = totDis2;
       continue;
    end
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    distance = -delta*deciSize/2;
    totDis2 = totDis2 + distance;
    
    %calculate the x and y
    data1(i) = sqrt((totDis1^2 - L1^2)*(totDis2^2 - L2^2)*((L1+L2)^2 - (totDis1 - totDis2)^2))/(2*totDis1*L2 + 2*totDis2*L1);
    data2(i) = (L2*totDis1^2 - L1*totDis2^2 - totDis1^2*totDis2+totDis2^2*totDis1)/2*(totDis1*L2 + totDis2*L1);

    
    set(pic,'xdata',data1,'ydata',data2);
    drawnow;
    
            
end

