function risPhsSet = gen_ris_phaseSet(Mris, b, errStd)
Nb = 2^b;
errStd = errStd/180*pi;
risPhsSet = zeros(Mris,Nb);
stdPsSet = (0:Nb-1)'/Nb*2*pi;
for mm = 1 : Mris
    errStdSet = [0;(2*rand(Nb-1,1)-1)*errStd];
    geniePs = stdPsSet+errStdSet;
    risPhsSet(mm,:) = exp(1j*geniePs);
end
end