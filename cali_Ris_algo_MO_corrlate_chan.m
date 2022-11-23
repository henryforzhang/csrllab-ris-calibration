function risEstPhsSet = cali_Ris_algo_MO_corrlate_chan(Heff,O,L,Mr,Mris,stdPhsSet,permtrSet)
% ***************************************
%  estimate RIS phases
%  author - Wei Zhang
%  input: 
%            Heff: estimated channels
%            O: number of measurements per gear 
%            L: number of gears 
%            Mr: receiving antennas 
%            Mris: number of RIS elements
%            stdPhsSet: nominal phases
%            permtrSet: permutation matrice
%  output: 
%            risEstPhsSet:estimated phase shifts
%copyright - CSRL@Fudan,2022/11/23
%  ************************************
Yvec = gen_hvec(Heff,O,L,Mr);
%% calibration
risStdPhsSet = ones(Mris,1)*stdPhsSet;
risEstPhsSet = risStdPhsSet;
if  Mris == 64
    iterNum = 5e3;
else 
    iterNum = 1e2;
end
testData = zeros(1,iterNum);
for ii = 1 : iterNum
    risEstEmplyed = gen_multi_d(risEstPhsSet,permtrSet);
    %% channel estimation
    HEst = wls_chan_est_fast_version(risEstEmplyed,Heff);
    E = gen_E(HEst,permtrSet);
    %% phase estimation
    F = risEstPhsSet.';
    f = F(:);
    testData(ii) = norm(Yvec-E*f,'fro')^2/(O*L*Mr);
    [f, cost] = my_manopt_method_correlate_chan(Yvec,E,f,Mr,L,O);
    F = reshape(f,L,Mris);
    risEstPhsSet = F.';
    if ii > 5 && abs(testData(ii)-testData(ii-5))/abs(testData(ii)) < 1e-2
        break;
    end
end
% figure;semilogy(testData)
end