clc;
clear;
close all;
load('voice-20200401-201127.mat');
t = reshape(1/48000:1/48000:10,480000,1);
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
Hm = mfilt.cicdecim(16,17,3);
data1 = zeros(250,1); 
data2 = zeros(250,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
var = zeros(numofFre,1);
numCount = 0;
baseband1 = zeros(250,numofFre);


TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;
f = 17150;
receivedSignal = R;
cosSignal = zeros(480000,numofFre);
sinSignal = zeros(480000,numofFre);

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
        
        [DVBaseband(:,freNum),DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = rawToBaseband(framesignal,sin1,cos1,Hm,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));
        
        
        [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));
       
        baseband1(i,freNum) = sum(DVBaseband(120,freNum));
    

    end
    
    data1(i) = IDFTDistanceMeasure(baseband1(i,:),numofFre,SoundSpeed);

end
x = 1/25:1/25:10;

plot(x,data1,'r');
xlabel('time , s');
ylabel('distance , m');
legend('预期距离' , 'IDFT测量距离');

load train;
sound(y,Fs);