function [distance] = IDFTDistanceMeasure(baseband1,numofFre,SoundSpeed)

    baseband = reshape(baseband1(1,:),numofFre,1);
    ift = dsp.IFFT;    
    y = ift(baseband);

    [~, indy] = max(abs(y));

    distance = (indy-1)*(SoundSpeed/((numofFre-1)*350));
end

