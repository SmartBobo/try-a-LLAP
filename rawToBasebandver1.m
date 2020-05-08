
function [I,Q,tempI,tempQ] = rawToBasebandver1(inputSignal,sin1,cos1,preI,preQ,deciSize,frameLength,i)
    tempI = zeros(480000,1);
    tempQ = zeros(480000,1);
    cicDecim = dsp.CICDecimator(16,17,3);
    InPhase = inputSignal.*cos1;
    Quard = inputSignal.*sin1;
    if i == 1
        tempI(1:i*frameLength) = InPhase;
        tempQ(1:i*frameLength) = Quard;
    else
        tempI(1:i*frameLength) = [preI(1:(i-1)*frameLength);InPhase];
        tempQ(1:i*frameLength) = [preQ(1:(i-1)*frameLength);Quard];
    end
%     disp(tempI);
    InPhaseFi = cicDecim(tempI(1:i*frameLength));
    QuardFi = cicDecim(tempQ(1:i*frameLength));
    I = double(InPhaseFi((i-1)*deciSize+1:i*deciSize));
    Q = double(QuardFi((i-1)*deciSize+1:i*deciSize));
   
end



