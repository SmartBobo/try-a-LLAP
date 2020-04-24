clc;
clear;
TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;

%prepare the microphone
Microphone = audioDeviceReader;
Microphone.SampleRate = 48000;
Microphone.SamplesPerFrame = 1920;
setup(Microphone);

%create streaming sine wave and cos wave object
numofFre = 16;
cosWave = cell(numofFre,1);
sinWave = cell(numofFre,1);
for i = 1:numofFre
    cosWave{i,1} =dsp.SineWave(1,(17150+i*350),pi/2,'SampleRate',48000,'SamplesPerFrame',1920);
    sinWave{i,1} =dsp.SineWave(1,(17150+i*350),pi,'SampleRate',48000,'SamplesPerFrame',1920);
end

%set CIC Filter
% cicDecim = dsp.CICDecimator(16,17,3);

%set variables
framelength = 1920;
deciSize = framelength/16;
DCvalueI = zeros(numofFre,1);
DCvalueR = zeros(numofFre,1);
maxR = zeros(numofFre,1);
maxI = zeros(numofFre,1);
minR = zeros(numofFre,1);
minI = zeros(numofFre,1);
totDis = 200;
I = zeros(deciSize,numofFre);
Q = zeros(deciSize,numofFre);
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
data1 = zeros(750,1);
data2 = zeros(750,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
variance = zeros(numofFre,1);
numCount = 0;
tempInPhase = zeros(framelength,numofFre);
tempQuard = zeros(framelength,numofFre);
%prepare the audio
dt = 1/48000;
t = 0:dt:30-dt;
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
pic = plot(axes1,0:1/25:30-1/25,data2);
set(axes1,'xlim', [0 30], 'ylim', [180 400]);
tic;
drawnow;

i = 0;
%可能还要考虑Latency due to sound card's output buffer: 0.046440 seconds
while toc <= 30
   [audioIn,Overrun] = step(Microphone);        % 采样
   if Overrun > 0
      warning('  数据溢出 %d 位\n',Overrun);
   end
   sumxy = 0;
   sumy = 0;
   numCount = 0;
   distance = 0;
   i = i + 1;
   disp(i);
   for freNum = 1:8    
        
        %将原始信号转换成Baseband
        [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBaseband(audioIn,sinWave{freNum,1}(),cosWave{freNum,1}(),tempInPhase(:,freNum),tempQuard(:,freNum),deciSize);
        
        
        %将两路信号用LEVD算法处理并求出DC vector
        [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(I,Q,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));


        %除去信号中的static vector 并将两路信号转变成baseband 
        DVInPhase = I - ESCInPhase; 
        DVQuard = Q - ESCQuard;
        DVBaseband(:,freNum) = ReImToComp(DVInPhase, DVQuard);
        % 求出相位同时根据threshold继续对DC vector进行处理

        [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));

        fre = fc + freNum * 350;


        %计算方差
        if freCount(freNum) == 1
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
        data1(i) = 0;
        data2(i) = totDis;
        continue;
    end
    
    %least square linear regression
    deltax = numofFre*((deciSize-1)*deciSize*(2*deciSize-1)/6-(deciSize-1)*deciSize*(deciSize-1)/4);
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    
    temp1 = delta*t1;
    varsum = 0;
    
    %get variance of each frequency
    for freNum = 1:numofFre
        if freCount(freNum) == 1
            temp2 = temp1 - temp(:,freNum);
            variance(freNum) = sum(temp2.^2);
            varsum = varsum + variance(freNum);
        end
    end
    
    for freNum = 1:numofFre
        if variance(freNum) > varsum/numCount
            freCount(freNum) = 0;    
        end
    end
    
    sumxy = 0;
    sumy = 0;
    numCount = 0;
    
    %Caluculate the distance based on the frequency after removing the ignore frequency
    for freNum = 1:numofFre
       if freCount(freNum) == 1
           sumy = sum(temp(:,freNum))+sumy;
           sumxy = sum(temp(:,freNum).*t1)+sumxy;
           numCount = numCount + 1;
       end

    end
    if numCount == 0
       data1(i) = 0;
       data2(i) = totDis;
       continue;
    end
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
    distance = -delta*deciSize/2;
    totDis = totDis + distance;
    data2(i) = totDis;
    data1(i) = distance;
    
    set(pic,'ydata',data2);
    drawnow;
    
            
end

