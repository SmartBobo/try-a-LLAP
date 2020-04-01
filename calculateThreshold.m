load('staticSignal17150-22750.mat');
t = reshape(1/48000:1/48000:10,480000,1);
numofFre = 8;

TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;
f = 17150;
cosSignal = zeros(480000,numofFre);
sinSignal = zeros(480000,numofFre);

%8个频率下的signal
for freNum = 1:numofFre
    fre = f + freNum * 350;
    
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000,1);
end
receivedSignal = R;
framelength = 1920;
deciSize = framelength/16;
DCvalueI = zeros(numofFre,1);
DCvalueR = zeros(numofFre,1);
maxR = zeros(numofFre,1);
maxI = zeros(numofFre,1);
minR = zeros(numofFre,1);
minI = zeros(numofFre,1);
totDis = 0;
InPhasePro = zeros(30000,numofFre);
QuardPro = zeros(30000,numofFre);
DVBaseband = zeros(30000,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
Hm = mfilt.cicdecim(16,17,3);
data1 = zeros(250,8); 
data2 = zeros(250,8);
data3 = zeros(250,8);




sumxy = 0;
sumy = 0;
freCount = zeros(numofFre,1);
temp = zeros(deciSize,numofFre);
numCount = 0;
framesignal = receivedSignal;
for freNum = 1:numofFre
    cos1 = cosSignal(:,freNum);
    sin1 = sinSignal(:,freNum);
    InPhase = framesignal.*sin1;
    Quard = framesignal.*cos1;

    %将两路信号用CICfilter处理
    InPhaseFi = filter(Hm,InPhase);
    QuardFi = filter(Hm,Quard);

    InPhasePro(:,freNum) = double(InPhaseFi);
    QuardPro(:,freNum) = double(QuardFi);
    
    DVBaseband(:,freNum) = ReImToComp(InPhasePro(:,freNum),QuardPro(:,freNum));

end
for i = 1:250
    InPhase1 = InPhasePro((i-1)*120+1:i*120,:);
    Quard1 = QuardPro((i-1)*120+1:i*120,:);
    Baseband1 = DVBaseband((i-1)*120+1:i*120,:);
    data1(i,:) = reshape(std(InPhase1),8,1);
    data2(i,:) = reshape(std(Quard1),8,1);
    data3(i,:) = reshape(var(Baseband1),8,1);
end

disp(['The standard deviation for In Phase signal is ' num2str(min(std(InPhasePro)))]);
disp(['The standard deviation for Quard signal is ' num2str(min(std(QuardPro)))]);
disp(['The standard deviation for Baseband signal is ' num2str(min(std(DVBaseband)))]);
