function [DVBaseband,DCvalueR,DCvalueI,maxR,minR,maxI,minI] = rawToBaseband(inputSignal,sin1,cos1,Hm,DCvalueR,DCvalueI,maxR,minR,maxI,minI)
    InPhase = inputSignal.*sin1;
    Quard = inputSignal.*cos1;

    %将两路信号用CICfilter处理
    InPhaseFi = filter(Hm,InPhase);
    QuardFi = filter(Hm,Quard);

    InPhasePro = double(InPhaseFi);
    QuardPro = double(QuardFi);

    %将两路信号用LEVD算法处理并求出DC vector
    [ESCInPhase,ESCQuard,DCvalueR,DCvalueI,maxR,minR,maxI,minI] = LEVD2(InPhasePro,QuardPro,DCvalueR,DCvalueI,maxR,minR,maxI,minI);

    %除去信号中的static vector 并将两路信号转变成baseband 
    DVInPhase = InPhasePro - ESCInPhase; 
    DVQuard = QuardPro - ESCQuard;
    DVBaseband = ReImToComp(DVInPhase, DVQuard);
end

