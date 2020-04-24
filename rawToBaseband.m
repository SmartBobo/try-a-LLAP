function [I,Q,tempI,tempQ] = rawToBaseband(inputSignal,sin1,cos1,preI,preQ,deciSize)
    
    cicDecim = dsp.CICDecimator(16,17,3);
    InPhase = inputSignal.*cos1;
    Quard = inputSignal.*sin1;
    
    tempI = InPhase;
    tempQ = Quard;

    InPhaseFi = cicDecim([preI;InPhase]);
    QuardFi = cicDecim([preQ;Quard]);
    I = double(InPhaseFi((2-1)*deciSize+1:2*deciSize));
    Q = double(QuardFi((2-1)*deciSize+1:2*deciSize));
   
end

