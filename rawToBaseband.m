function [DVBaseband,DCvalueR,DCvalueI,maxR,minR,maxI,minI] = rawToBaseband(inputSignal,sin1,cos1,Hm,DCvalueR,DCvalueI,maxR,minR,maxI,minI)
    InPhase = inputSignal.*sin1;
    Quard = inputSignal.*cos1;

    %����·�ź���CICfilter����
    InPhaseFi = filter(Hm,InPhase);
    QuardFi = filter(Hm,Quard);

    InPhasePro = double(InPhaseFi);
    QuardPro = double(QuardFi);

    %����·�ź���LEVD�㷨�������DC vector
    [ESCInPhase,ESCQuard,DCvalueR,DCvalueI,maxR,minR,maxI,minI] = LEVD2(InPhasePro,QuardPro,DCvalueR,DCvalueI,maxR,minR,maxI,minI);

    %��ȥ�ź��е�static vector ������·�ź�ת���baseband 
    DVInPhase = InPhasePro - ESCInPhase; 
    DVQuard = QuardPro - ESCQuard;
    DVBaseband = ReImToComp(DVInPhase, DVQuard);
end

