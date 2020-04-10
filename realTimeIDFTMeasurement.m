clc;
clear;
close all;
load('16Fre来回移动5次.mat');
t = reshape(0:1/48000:10-1/48000,480000,1);
numofFre = 16;
framelength = 1920;
deciSize = framelength/16;
DCvalueI = zeros(numofFre,1);
DCvalueR = zeros(numofFre,1);
maxR = zeros(numofFre,1);
maxI = zeros(numofFre,1);
minR = zeros(numofFre,1);
minI = zeros(numofFre,1);
totDis = 0;
baseband = zeros(16,1);
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);

data1 = zeros(250,1); 
data2 = zeros(250,1);
data3 = zeros(250,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
var = zeros(numofFre,1);
numCount = 0;
baseband1 = zeros(250,numofFre);
baseband2 = zeros(250,numofFre);
baseband3 = zeros(250,numofFre);
DVInPhase = zeros(deciSize,numofFre);
DVQuard = zeros(deciSize,numofFre);

TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;
f = 17150;
receivedSignal = R;
cosSignal = zeros(480000,numofFre);
sinSignal = zeros(480000,numofFre);
tempInPhase = zeros(480000,numofFre);
tempQuard = zeros(480000,numofFre);

%16个频率下的signal
for freNum = 1:numofFre
    fre = f + freNum * 350;
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000,1);
end

for i = 1:250
    framesignal=receivedSignal((i-1)*framelength+1:i*framelength);

    for freNum = 1:numofFre
        cos1 = cosSignal(((i-1)*framelength+1:i*framelength),freNum);
        sin1 = sinSignal(((i-1)*framelength+1:i*framelength),freNum);
        
        [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBaseband(framesignal,sin1,cos1,tempInPhase(:,freNum),tempQuard(:,freNum),i,deciSize,framelength);
       
        [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(I,Q,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));
        DVInPhase(:,freNum) = I - ESCInPhase;
        DVQuard(:,freNum) = Q - ESCQuard;
        DVBaseband(:,freNum) = ReImToComp(DVInPhase(:,freNum), DVQuard(:,freNum));
        
        [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));
       
        baseband1(i,freNum) = DVBaseband(120,freNum);
        baseband2(i,freNum) = DVBaseband(120,freNum)-DVBaseband(1,freNum);
        baseband3(i,freNum) = sum(DVBaseband(:,freNum));
    

    end
    
    data1(i) = IDFTDistanceMeasure(baseband1(i,:),numofFre,SoundSpeed);
    data2(i) = IDFTDistanceMeasure(baseband1(i,:),numofFre,SoundSpeed);
    data3(i) = IDFTDistanceMeasure(baseband1(i,:),numofFre,SoundSpeed);

end
x = 1/25:1/25:10;
figure;
plot(x,data1,'r');
title('只用每一帧的最后一个点计算');
xlabel('time , s');
ylabel('distance , m');
legend('IDFT测量距离');
figure;
plot(x,data2,'r');
title('用每一帧的最后一个点减去第一个点的差值计算');
xlabel('time , s');
ylabel('distance , m');
legend('IDFT测量距离');
figure;
plot(x,data3,'r');
title('用每一帧所有点的和去计算');
xlabel('time , s');
ylabel('distance , m');
legend('IDFT测量距离');

load train;
sound(y,Fs);