function y = gen_hvec(Heff,O,Nb,Mr)
oNbMr = O*Nb*Mr;
nbMr = Nb*Mr;
y = zeros(oNbMr,1);
for oo = 1 : O
    hTmp = Heff(:,(oo-1)*Nb+(1:Nb)).';
    y((oo-1)*nbMr+(1:nbMr)) = hTmp(:);
end

end