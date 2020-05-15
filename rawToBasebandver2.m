function [I,Q,tempI,tempQ] = rawToBasebandver2(inputSignal,sin1,cos1,preI,preQ,deciSize,framelength)
    
    cicDecim = dsp.CICDecimator(16,17,3);
    InPhase = inputSignal.*cos1;
    Quard = inputSignal.*sin1;
    
    tempI = InPhase(3/5*framelength+1:framelength);
    tempQ = Quard(3/5*framelength+1:framelength);
    
    %上一帧拼接这一帧 去做cicfilter
    InPhaseFi = cicDecim([preI;InPhase]);
    QuardFi = cicDecim([preQ;Quard]);
    I = double(InPhaseFi(2*deciSize/5+1:7*deciSize/5));
    Q = double(QuardFi(2*deciSize/5+1:7*deciSize/5));
   
end

