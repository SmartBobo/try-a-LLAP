clc;
clear;
TEMPERATURE = 30;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;

%prepare the microphone
Microphone = audioDeviceReader;
Microphone.SampleRate = 48000;
Microphone.SamplesPerFrame = 1920;
Microphone.Device ='麦克风 (Realtek(R) Audio)';


% % create streaming sine wave and cos wave object
% numofFre = 10;
% cosWave = cell(numofFre,1);
% sinWave = cell(numofFre,1);
% for i = 1:numofFre
%     cosWave{i,1} =dsp.SineWave(1,(15650+i*350),pi/2,'SampleRate',48000,'SamplesPerFrame',1920);
%     sinWave{i,1} =dsp.SineWave(1,(15650+i*350),pi,'SampleRate',48000,'SamplesPerFrame',1920);
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
% set CIC Filter
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
totDis = 200;
I = zeros(deciSize,numofFre); 
Q = zeros(deciSize,numofFre);
DVBaseband = zeros(deciSize,numofFre);
temp = zeros(deciSize,numofFre);
temp1 = zeros(deciSize,1);
temp2 = zeros(deciSize,1);
calidata = zeros(10,1);
data1 = zeros(750,1);
data2 = zeros(750,1);
t1 = reshape(0:deciSize-1,deciSize,1);
freCount = zeros(numofFre,1);
variance = zeros(numofFre,1);
numCount = 0;
cali = zeros(4,1);
tempdistance = zeros(50,1);
tempInPhase = zeros(framelength*2/5,numofFre);
tempQuard = zeros(framelength*2/5,numofFre);
baseband1 = zeros(250,numofFre);
cicDecim = dsp.CICDecimator(16,17,3);
%prepare the audio
dt = 1/48000;
t = 0:dt:40-dt;
fc = 15650;
A = 0.2;
y = 0; %A*cos(2*pi*fc*t);
for i = 1:numofFre
    y = y + A*cos(2*pi*(fc+350*i)*t);
end
y = y/numofFre;
player2 = audioplayer(y,48000,16,5);
tic;
play(player2);
delay = toc;
setup(Microphone);

%initiate the figure
axes1 = subplot(1,1,1);
pic = plot(axes1,0:1/25:30-1/25,data1);
set(axes1,'xlim', [0 30], 'ylim', [0 1]);
tic;
drawnow;
idft = 1;
disp('start recording');
i = 0;
k = 0;
flag = 0;
j = 1;
%可能还要考虑Latency due to sound card's output buffer: 0.046440 seconds
while toc<30
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
   if idft == 1
        for freNum = 1:10    

            %将原始信号转换成Baseband
            [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBasebandver2(audioIn,sinSignal((i-1)*framelength+1-k:i*framelength-k,freNum),cosSignal((i-1)*framelength+1-k:i*framelength-k,freNum),tempInPhase(:,freNum),tempQuard(:,freNum),cicDecim);
 
            DVBaseband(:,freNum) = ReImToComp(I,Q);
            
            
        end
%         idft = 0;
        if flag ~= 1
            k = k + 1;
            x = 1;
            for a = 1:39:120 
                index = IDFTDistanceMeasure(DVBaseband(a,:),numofFre,SoundSpeed);
                if index == 3
                    cali(x) = 1;
                    x = x + 1;
                end
            end
        end
        
        if mean(cali(1:4)) == 1
            flag = 1;
            idft = 0;
            disp('Calibration finished and place your hand');
            data1 = zeros(750,1);
            continue;
        else
            cali = zeros(4,1);
        end
        

        data1(i) = (index-1)*SoundSpeed/((10-1)*350);   


        
   end
    
    if idft == 0
        for freNum = 1:10    

            %将原始信号转换成Baseband
            [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBasebandver2(audioIn,sinSignal((i-1)*framelength+1-k:i*framelength-k,freNum),cosSignal((i-1)*framelength+1-k:i*framelength-k,freNum),tempInPhase(:,freNum),tempQuard(:,freNum),cicDecim);

            %将两路信号用LEVD算法处理并求出DC vector
            [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(I,Q,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));


            %除去信号中的static vector 并将两路信号转变成baseband 
            DVInPhase = I - ESCInPhase; 
            DVQuard = Q - ESCQuard;
            DVBaseband(:,freNum) = ReImToComp(DVInPhase,DVQuard);
            % 求出相位同时根据threshold继续对DC vector进行处理

            [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));

            baseband1(j,freNum) = DVBaseband(120,freNum);

        end
        index = IDFTDistanceMeasure(baseband1(j,:),numofFre,SoundSpeed);
        %当j小于等于50的时候 进行对初始位置的判断 后续每用idft求出一个新的点就与前面比较
        if j == 50
            idftdistance = (index-1)*SoundSpeed/((10-1)*350);
%             idftdistance = mode(tempdistance);
        elseif j > 50
            data1(j) = (index-1)*SoundSpeed/((10-1)*350);
            if abs(data1(j)-idftdistance) < 3*SoundSpeed/((10-1)*350) 
                
                idftdistance = data1(j);
            else
                data1(j) = data1(j-1);
            end
        end
        

        
        j = j + 1;
   end
    
   set(pic,'ydata',data1);
   drawnow;
            
end
release(Microphone);


