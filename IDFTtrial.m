clc;
clear;
TEMPERATURE = 16;
SoundSpeed = 331.3 + 0.606 * TEMPERATURE;

%prepare the microphone
Microphone = audioDeviceReader;
Microphone.SampleRate = 48000;
Microphone.SamplesPerFrame = 1920;
Microphone.Device ='麦克风 (Realtek(R) Audio)';
setup(Microphone);

%create streaming sine wave and cos wave object
% numofFre = 16;
% cosWave = cell(numofFre,1);
% sinWave = cell(numofFre,1);
% for i = 1:numofFre
%     cosWave{i,1} =dsp.SineWave(1,(17150+i*350),pi/2,'SampleRate',48000,'SamplesPerFrame',1920);
%     sinWave{i,1} =dsp.SineWave(1,(17150+i*350),pi,'SampleRate',48000,'SamplesPerFrame',1920);
% end

%另一种 streaming方法
numofFre = 16;
t = reshape(0:1/48000:10*4-1/48000,480000*4,1);
cosSignal = zeros(480000*4,numofFre);
sinSignal = zeros(480000*4,numofFre);
f = 17150;
for freNum = 1:16
    fre = f + freNum * 350;
    cosSignal(:,freNum) = reshape(cos(2*pi*fre*t),480000*4,1);
    sinSignal(:,freNum) = reshape(-sin(2*pi*fre*t),480000*4,1);
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
tempInPhase = zeros(framelength*2/5,numofFre);
tempQuard = zeros(framelength*2/5,numofFre);
baseband1 = zeros(250,numofFre);

%prepare the audio
dt = 1/48000;
t = 0:dt:40-dt;
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
pic = plot(axes1,0:1/25:30-1/25,data1);
set(axes1,'xlim', [0 30], 'ylim', [0 1]);
tic;
drawnow;
idft = 0;
disp('start recording');
i = 0;
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
        for freNum = 1:16    

            %将原始信号转换成Baseband
            [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBasebandver2(audioIn,sinSignal((i-1)*framelength+1:i*framelength,freNum),cosSignal((i-1)*framelength+1:i*framelength,freNum),tempInPhase(:,freNum),tempQuard(:,freNum),deciSize,framelength);


            %将两路信号用LEVD算法处理并求出DC vector
            [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(I,Q,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));


     
            DVBaseband(:,freNum) = ReImToComp(I,Q);
            
            baseband1(i,freNum) = DVBaseband(120,freNum);

        end
        idft = 0;
        [idftdistance,index] = IDFTDistanceMeasure(baseband1(i,:),numofFre,SoundSpeed);
        delay = (1 - index)/(15*350);
%         for i = 1:numofFre
%             cosWave{i,1} =dsp.SineWave(1,(17150+i*350),pi/2-delay*2*pi*(17150+i*350),'SampleRate',48000,'SamplesPerFrame',1920);
%             sinWave{i,1} =dsp.SineWave(1,(17150+i*350),pi-delay*2*pi*(17150+i*350),'SampleRate',48000,'SamplesPerFrame',1920);
%         end
        for freNum = 1:16
            fre = f + freNum * 350;
            cosSignal(:,freNum) = reshape(cos(2*pi*fre*(t-delay)),480000*4,1);
            sinSignal(:,freNum) = reshape(-sin(2*pi*fre*(t-delay)),480000*4,1);
        end
        data1(i) = idftdistance;
   end
    
    if idft == 0
        for freNum = 1:16    

            %将原始信号转换成Baseband
            [I,Q,tempInPhase(:,freNum),tempQuard(:,freNum)] = rawToBasebandver2(audioIn,sinSignal((i-1)*framelength+1:i*framelength,freNum),cosSignal((i-1)*framelength+1:i*framelength,freNum),tempInPhase(:,freNum),tempQuard(:,freNum),deciSize,framelength);


            %将两路信号用LEVD算法处理并求出DC vector
            [ESCInPhase,ESCQuard,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum)] = LEVD2(I,Q,DCvalueR(freNum),DCvalueI(freNum),maxR(freNum),minR(freNum),maxI(freNum),minI(freNum));


            %除去信号中的static vector 并将两路信号转变成baseband 
            DVInPhase = I - ESCInPhase; 
            DVQuard = Q - ESCQuard;
            DVBaseband(:,freNum) = ReImToComp(DVInPhase,DVQuard);
            % 求出相位同时根据threshold继续对DC vector进行处理

            [ph,DCvalueR(freNum),DCvalueI(freNum),freCount(freNum)] = DCprocess(DVBaseband(:,freNum),maxR(freNum),minR(freNum),DCvalueR(freNum),maxI(freNum),minI(freNum),DCvalueI(freNum),freCount(freNum));

            baseband1(i,freNum) = DVBaseband(120,freNum);

        end
        [idftdistance,index] = IDFTDistanceMeasure(baseband1(i,:),numofFre,SoundSpeed);
        data1(i) = idftdistance;
   end
    
   set(pic,'ydata',data1);
   drawnow;
            
end
% for i = 1:numofFre
%     cosWave{i,1} =dsp.SineWave(1,(17150+i*350),pi/2-delay*2*pi*(17150+i*350),'SampleRate',48000,'SamplesPerFrame',1920);
%     sinWave{i,1} =dsp.SineWave(1,(17150+i*350),pi-delay*2*pi*(17150+i*350),'SampleRate',48000,'SamplesPerFrame',1920);
% end

