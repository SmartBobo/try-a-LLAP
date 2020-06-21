clc;
clear;
TEMPERATURE = 30;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;

%prepare the microphone
Microphone1 = audioDeviceReader;
Microphone1.SampleRate = 48000;
Microphone1.SamplesPerFrame = 1920;
Microphone1.Device ='麦克风 (Realtek(R) Audio)';
%Microphone1.Device ='麦克风 (Sound Blaster R3/A6U)';

% Microphone2 = audioDeviceReader;
% Microphone2.SampleRate = 48000;
% Microphone2.SamplesPerFrame = 1920;
% Microphone2.Device ='麦克风 (2- USB Audio Device)';


% %create streaming sine wave and cos wave object
% numofFre = 8;
% cosWave = cell(numofFre,1);
% sinWave = cell(numofFre,1);
% for i = 1:numofFre
%     cosWave{i,1} =dsp.SineWave(1,(17150+i*350),pi/2,'SampleRate',48000,'SamplesPerFrame',1920);
%     sinWave{i,1} =dsp.SineWave(1,(17150+i*350),pi,'SampleRate',48000,'SamplesPerFrame',1920);
% end

%另一种 streaming方法
numofFre = 10;
t = reshape(0:1/48000:10*4-1/48000,480000*4,1);
cosSignal = zeros(480000*4,numofFre);
sinSignal = zeros(480000*4,numofFre);
f = 15650;
for freNum = 1:10
    fre = f + freNum * 350;
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000*4,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000*4,1);
end

%set CIC Filter
cicDecim = dsp.CICDecimator(16,17,3);

%set variables
framelength = 1920;
deciSize = framelength/16;
DCvalueI = zeros(numofFre,1);
DCvalueR = zeros(numofFre,1);
maxR = zeros(numofFre,1);
maxI = zeros(numofFre,1);
minR = zeros(numofFre,1);
minI = zeros(numofFre,1);
I = zeros(deciSize,numofFre); 
Q = zeros(deciSize,numofFre);
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
IDFTdata1 = zeros(750,1);
data2 = zeros(750,1);
tempdistance = zeros(50,1);
idftTemp = zeros(750,numofFre);
cali = zeros(4,1);
baseband1 = zeros(250,numofFre);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
variance = zeros(numofFre,1);
numCount = 0;
tempInPhase = zeros(2*framelength/5,numofFre);
tempQuard = zeros(2*framelength/5,numofFre);

%prepare the audio
dt = 1/48000;
t = 0:dt:35-dt;
fc = 15650;
A = 0.2;
y = 0; %A*cos(2*pi*fc*t);
for i = 1:numofFre
    y = y + A*cos(2*pi*(fc+350*i)*t);
end
y = y/numofFre;
player2 = audioplayer(y,48000,16,5);
play(player2);

%initiate the figure
axes1 = subplot(1,1,1);
pic = plot(axes1,0:1/25:30-1/25,data2);
set(axes1,'xlim', [0 30], 'ylim', [0 700]);
tic;
drawnow;
idft = 1;
disp('Start Recording');
i = 0;
k = 0;
flag = 0;
loopcount = 0;
j = 1;
idftdistance = 0;
%可能还要考虑Latency due to sound card's output buffer: 0.046440 seconds
while toc <= 30
   [audioIn1,Overrun1] = step(Microphone1);        % 采样
%    [audioIn2,Overrun2] = step(Microphone2);
   if Overrun1 > 0
      warning('  Microphone1 数据溢出 %d 位\n',Overrun1);
   end
%    if Overrun2 > 0
%       warning('  Microphone2 数据溢出 %d 位\n',Overrun2);
%    end
     sumxy = 0;
   sumy = 0;
   numCount = 0;
   distance = 0;
   i = i + 1;
   disp(i);
   if idft == 1
        for freNum = 1:10    

            %将原始信号转换成Baseband
            [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBasebandver2(audioIn1,sinSignal((i-1)*framelength+1-k:i*framelength-k,freNum),cosSignal((i-1)*framelength+1-k:i*framelength-k,freNum),tempInPhase(:,freNum),tempQuard(:,freNum),cicDecim);
 
            DVBaseband(:,freNum) = ReImToComp(I,Q);
            
            
        end
%         idft = 0;
        if flag ~= 1
            k = k + 1;
            x = 1;
            for a = 1:39:120 
                index = IDFTDistanceMeasure(DVBaseband(a,:),numofFre,SoundSpeed);
                if index == 2
                    cali(x) = 1;
                    x = x + 1;
                end
            end
        end
        
        if mean(cali(1:4)) == 1
            flag = 1;
            idft = 0;
            disp('Calibration finished and place your hand');
            IDFTdata1 = zeros(750,1);
            continue;
        else
            cali = zeros(4,1);
        end
        

        IDFTdata1(i) = (index-1)*SoundSpeed/((10-1)*350)*1000;   


        
   end
    
    if idft == 0
        for freNum = 1:10    

            %将原始信号转换成Baseband
            [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBasebandver2(audioIn1,sinSignal((i-1)*framelength+1-k:i*framelength-k,freNum),cosSignal((i-1)*framelength+1-k:i*framelength-k,freNum),tempInPhase(:,freNum),tempQuard(:,freNum),cicDecim);

            %将两路信号用LEVD算法处理并求出DC vector
            [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(I,Q,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));


            %除去信号中的static vector 并将两路信号转变成baseband 
            DVInPhase = I - ESCInPhase; 
            DVQuard = Q - ESCQuard;
            DVBaseband(:,freNum) = ReImToComp(DVInPhase,DVQuard);
            % 求出相位同时根据threshold继续对DC vector进行处理

            [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));

            baseband1(j,freNum) = DVBaseband(120,freNum);
            
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
        index = IDFTDistanceMeasure(baseband1(j,:),numofFre,SoundSpeed);
        
        %j在0到44的时候给手一个放上去的时间,44到49多个采样点确保初始位置的准确性
        %当j等于50的时候 进行对初始位置的记录 后续每用idft求出一个新的点就与前面比较
        if j > 44 && j <50
            idftdistance = idftdistance+(index-1)*SoundSpeed/((10-1)*350)*1000;
        elseif j == 50
            idftdistance = idftdistance/5;
            totDis = idftdistance;
        elseif j > 50
            IDFTdata1(j) = (index-1)*SoundSpeed/((10-1)*350)*1000;
%             if abs(IDFTdata1(j)-idftdistance) < 2*SoundSpeed/((10-1)*350)*1000 
%                 
%                 idftdistance = IDFTdata1(j);
%             else
%                 IDFTdata1(j) = IDFTdata1(j-1);
%             end
        end
        
        if j > 50
            %find ignore frequency
            if numCount == 0
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
               data2(i) = totDis;
               continue;
            end
            delta = (sumxy-sumy*(deciSize-1)/2)/deltax*numofFre/numCount;
            distance = -delta*deciSize/2;
            
            %按周期以及计算出来的手的移动速度为标准判断是否导入idft测算出来的距离
            if (mod(j,50) == 0 && IDFTdata1(j) == IDFTdata1(j-1)) || (distance > 3 && loopcount == 0)
                totDis = IDFTdata1(j);
                loopcount = loopcount + 1;
                if loopcount == 10
                    loopcount = 0;
                end
            else
                totDis = totDis + 2*distance;
            end
            data2(j) = totDis;
        end
        j = j + 1;
    end
    set(pic,'ydata',data2);
    drawnow;
            
end

