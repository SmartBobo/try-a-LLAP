function Comp = ReImToComp(Re, Im)
size3 = length(Re);
Comp = zeros(size3,1);
for a = 1:size3
    Comp(a) = Re(a) + Im(a)*1i;
end
end

