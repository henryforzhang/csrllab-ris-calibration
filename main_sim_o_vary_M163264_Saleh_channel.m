%% full duplex communication system
%
%
clear;
% basic setting
%% RIS setting
MrisSet = [16 32 64];
MrisSetLen = length(MrisSet);
MrisMax = max(MrisSet);
%% Bs setting
Mr = 8;
Mt = 1;
% pathloss
pilotLen = 100;
pilotGain = pow2db(pilotLen);
gain = pilotGain;

snrDb = 20;
noisePow = db2pow(-(snrDb+gain));


b = 4;
errStd = 20;
L = 2^b;
OSet = (9:16);
OSetMax = max(OSet);
OSetMaxNb = OSetMax*L;
OSetLen = length(OSet);

stdPsSet = (0:L-1)/L*2*pi;
stdPhsSet = exp(1j*stdPsSet);

risGeniePhsSetMax = gen_ris_phaseSet(MrisMax, b, errStd);
permtrSetMax = gen_permtr_mtr(MrisMax,L,OSetMax);
risEmplyedMax = gen_multi_d(risGeniePhsSetMax,permtrSetMax);

%% channel
Nc = 3; % number of cluster
Nray = 5; % number of rays per cluster
[HbrRSet,hrBtSet] = save_channel(MrisSet,Mr,Nc,Nray,8);
%% crb
crbSet = zeros(MrisSetLen,OSetLen);
for mm = 1 : MrisSetLen
    Mris = MrisSet(mm);
    HrBt = hrBtSet{mm};
    HbrR = HbrRSet{mm};
    Hp = diag(HrBt);
    Hcas = HbrR*Hp;
    risGeniePhsSet = risGeniePhsSetMax(1:Mris,:); 
    for oo = 1 : OSetLen
        O = OSet(oo);
        oNb = O*L;
        permtrSet = permtrSetMax(1:O,1:Mris,:,:);
        fisher = cal_CRB_of_RIS_cali_correlate_chan(permtrSet,risGeniePhsSet,Hcas);
        fisherInv = fisher^(-1);
        fisherInvDiag = diag(fisherInv);
        crbSet(mm,oo) = noisePow/2*mean(fisherInvDiag(1:(Mris*(L-1))));
    end
end

%% rmse
mseSet = zeros(MrisSetLen,OSetLen);
nnMonte = 1e2;
for nn = 1 : nnMonte
    offUnitNoiseTl = sqrt(1/2)*(randn(Mr,1)+1j*randn(Mr,1));
    onUnitNoiseTl = sqrt(1/2)*(randn(Mr,OSetMaxNb)+1j*randn(Mr,OSetMaxNb));
    for mm = 1 : MrisSetLen
        Mris = MrisSet(mm);
        if nn == 1 || mod(nn,5) == 0
            fprintf(['\nThe Iteration Number: %d, the number of RIS: %d, at ', datestr(now,"HH:MM"),'\n'], nn, Mris);
        end
        HrBt = hrBtSet{mm};
        HbrR = HbrRSet{mm};
        Hp = diag(HrBt);
        Hcas = HbrR*Hp;
        risGeniePhsSet = risGeniePhsSetMax(1:Mris,:);
        
        for oo = 1 : OSetLen
            O = OSet(oo);
            oNb = O*L;
            permtrSet = permtrSetMax(1:O,1:Mris,:,:);
            onUnitNoise = onUnitNoiseTl(:,1:oNb);
            offUnitNoise = offUnitNoiseTl*ones(1,oNb);
            risEmplyed = risEmplyedMax(1:Mris,1:oNb);
            HeffDeNoise = Hcas*risEmplyed;
            Heff = HeffDeNoise+sqrt(noisePow)*(onUnitNoise-offUnitNoise);
            risEstPhsSet = cali_Ris_algo_MO_corrlate_chan(Heff,O,L,Mr,Mris,stdPhsSet,permtrSet);
            risEstPhsSetRmAm = remove_ambiguity(risEstPhsSet);
            mse = norm(angle(risEstPhsSetRmAm./risGeniePhsSet),'fro')^2/(Mris*(L-1));
            mseSet(mm,oo) = mseSet(mm,oo)+ mse;
        end
    end
end
mseSetAve = mseSet/nn;
rmseSet = sqrt(mseSetAve)*180/pi;
rcrbSet = sqrt(crbSet)*180/pi;
%% plot
plt1 = {'ro','kx','b^'};
plt2 = {'r-','k--','b:'};
figure;
for mm = 1 : MrisSetLen
    semilogy(OSet,rmseSet(mm,:),plt1{mm},'linewidth',2,'markersize',8);
    hold on
    semilogy(OSet,rcrbSet(mm,:),plt2{mm},'linewidth',2,'markersize',8);
end
grid on
xlabel('Number of measurements per gear');
ylabel('RMSE (^\circ)');
title([num2str(b),'-bit RIS, ', 'Phase Deviations [-',num2str(errStd),'\circ',', ',num2str(errStd),'\circ',']', ', Mr = ', num2str(Mr)])
legend('M_{ris} = 16, Estimate','M_{ris} = 16, CRB',...
    'M_{ris} = 32, Estimate','M_{ris} = 32, CRB',...
     'M_{ris} = 64, Estimate','M_{ris} = 64, CRB')










