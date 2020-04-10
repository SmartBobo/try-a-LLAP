function [I,Q,tempI,tempQ] = rawToBaseband(inputSignal,sin1,cos1,tempI,tempQ,i,deciSize,framelength)
    
    cicDecim = dsp.CICDecimator(16,17,3);
    InPhase = inputSignal.*cos1;
    Quard = inputSignal.*sin1;
    
    tempI((i-1)*framelength+1:i*framelength) = InPhase;
    tempQ((i-1)*framelength+1:i*framelength) = Quard;
    
    %将两路信号用CICfilter处理
    InPhaseFi = cicDecim(tempI(1:i*framelength));
    QuardFi = cicDecim(tempQ(1:i*framelength));

    I = double(InPhaseFi((i-1)*deciSize+1:i*deciSize));
    Q = double(QuardFi((i-1)*deciSize+1:i*deciSize));

    
end

