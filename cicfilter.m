function [output] = cicfilter(input,deciSize,dec,delay)
    cicBuffer = zeros(deciSize+delay,1);
    cicBuffer1 = zeros(deciSize+delay,1);
    cicBuffer2 = zeros(deciSize+delay,1);
    cicBuffer3 = zeros(deciSize+delay,1);
    output = zeros(deciSize,1);
    for i = 1:deciSize
        cicBuffer(delay+i) = sum(input((i-1)*dec+1:i*dec,1)); 
    end
    for i = 1:deciSize
        cicBuffer1(delay+i) = sum(cicBuffer(i:i+delay-1,1));
    end
    for i = 1:deciSize
        cicBuffer2(delay+i) = sum(cicBuffer1(i:i+delay-1,1));
    end
    for i = 1:deciSize
        cicBuffer3(delay+i) = sum(cicBuffer2(i:i+delay-1,1));
    end
    for i = 1:deciSize
        output(i) = sum(cicBuffer3(i:i+delay-1,1));
    end
    
    
end

