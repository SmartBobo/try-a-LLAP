clc;
clear;
load('17500移动3次.mat');
t = reshape(0:1/48000:10-1/48000,480000,1);
numofFre = 1;

TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;
f = 17150;
cosSignal = zeros(480000,numofFre);
sinSignal = zeros(480000,numofFre);

%8个频率下的signal
for freNum = 1:1
    fre = f + freNum * 350;
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000,1);
end
I = zeros(30000,1);
Q = zeros(30000,1);
Iprocess = zeros(30000,1);
Qprocess = zeros(30000,1);
receivedSignal = R;
framelength = 1920;
deciSize = framelength/16;
DCvalueI = zeros(numofFre,1);
DCvalueR = zeros(numofFre,1);
maxR = zeros(numofFre,1);
maxI = zeros(numofFre,1);
minR = zeros(numofFre,1);
minI = zeros(numofFre,1);
totDis = 250;
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
cicDecim = dsp.CICDecimator(16,17,3);
data1 = zeros(250,1); 
data2 = zeros(250,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
variance = zeros(numofFre,1);
numCount = 0;



tempInPhase = zeros(480000,numofFre);
tempQuard = zeros(480000,numofFre);

InPhase = receivedSignal.*cosSignal(:,1);
InPhaseFi = cicDecim(InPhase);
In = double(InPhaseFi);

for i = 1:250
    framesignal=receivedSignal((i-1)*framelength+1:i*framelength);
    sumxy = 0;
    sumy = 0;
    freCount = 0;
    freNum = 1;
    cos1 = cosSignal(((i-1)*framelength+1:i*framelength),freNum);
    sin1 = sinSignal(((i-1)*framelength+1:i*framelength),freNum);

    [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBaseband(framesignal,sin1(:,freNum),cos1(:,freNum),tempInPhase(:,freNum),tempQuard(:,freNum),i,deciSize,framelength);
        

    
    %将两路信号用LEVD算法处理并求出DC vector
    [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(I,Q,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));

    %除去信号中的static vector 并将两路信号转变成baseband 
    DVInPhase = I - ESCInPhase; 
    DVQuard = Q - ESCQuard;
    DVBaseband(:,freNum) = ReImToComp(DVInPhase, DVQuard);
    
    Baseband(((i-1)*deciSize+1:i*deciSize),freNum) = DVBaseband(:,freNum);
    

    % 求出相位同时根据threshold继续对DC vector进行处理

    [ph,DCvalueR(freNum),DCvalueI(freNum),freCount] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount);

    fre = f + freNum * 350;
   

    %linear regression

    for a = 1:deciSize
        temp(a,freNum) = ph(a) - ph(1); 
        temp(a,freNum) = temp(a,freNum)*SoundSpeed*1000/(2*pi*fre);
    end
    sumy = sum(temp(:,freNum))+sumy;
    sumxy = sum(temp(:,freNum).*t1)+sumxy;
    
    if freCount(freNum) == 0
        data1(i) = 0;
        data2(i) = totDis;
        continue;
    end

    deltax = ((deciSize-1)*deciSize*(2*deciSize-1)/6-(deciSize-1)*deciSize*(deciSize-1)/4);
    delta = (sumxy-sumy*(deciSize-1)/2)/deltax;

    

   
    distance = -delta*deciSize/2;    
    data1(i) = distance; 
    totDis = totDis + distance;
    data2(i) = totDis;

end

x = 1/3000:1/3000:10;
figure;
y = 1/25:1/25:10;
plot(y,data2,'r');
xlabel('time , s');
ylabel('distance, mm');
legend('moving distance');
figure;
plot(x,abs(Baseband));