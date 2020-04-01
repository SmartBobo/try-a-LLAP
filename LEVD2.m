function [outputR,outputI,newDCR,newDCI,LmaxR,LminR,LmaxI,LminI] = LEVD2(inputR,inputI,preDCR,preDCI,premaxR,preminR,premaxI,preminI)
%对输入的矩阵进行LEVD处理
size1 = length(inputR);
maxR = max(inputR);
minR = min(inputR);
maxI = max(inputI);
minI = min(inputI);
LmaxR = premaxR;
LmaxI = premaxI;
LminR = preminR;
LminI = preminI;
ampThr = 30000;%设为静态环境下接收的信号的标准差的三倍 16Fre是15000 8Fre是30000 
powerThr = 15000;
%get real part variance
temp = inputR - ones(size1,1)*inputR(1);
dsum = abs(sum(temp)/size1);
vsum = abs(sum(temp.^2)/size1);
%get imaginary part variance
temp = inputI - ones(size1,1)*inputI(1);
dsum = dsum + abs(sum(temp)/size1);
vsum = vsum + abs(sum(temp.^2)/size1);

newDCR = preDCR;
newDCI = preDCI;

if vsum + dsum * dsum >  powerThr %15000 is power threshold
    if (maxR > LmaxR || (maxR > LminR + ampThr && (LmaxR-LminR) > ampThr * 4)) %amplitude threshold
        LmaxR = maxR;

    end
    
    if (minR < LminR || (minR < LmaxR - ampThr && (LmaxR-LminR) > ampThr * 4))
        LminR = minR;

    end
    
    if (maxI > LmaxI || (maxI > LminI + ampThr && (LmaxI-LminI) > ampThr * 4)) %amplitude threshold
        LmaxI = maxI;

    end
    
    if (minI < LminI || (minI < LmaxI - ampThr && (LmaxI-LminI) > ampThr * 4))
        LminI = minI;

    end
    
    if ((LmaxR - LminR) > ampThr && (LmaxI - LminI) > ampThr)
        newDCR = (1-0.25)*preDCR + (LmaxR + LminR)/2*0.25;
        newDCI = (1-0.25)*preDCI + (LmaxI + LminI)/2*0.25;
    end
    

end

outputR = ones(size1,1)*newDCR;
outputI = ones(size1,1)*newDCI;

end