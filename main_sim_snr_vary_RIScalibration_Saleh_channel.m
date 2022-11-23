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
MrSet = [8];
MrSetLen = length(MrSet);
% pathloss
pilotLen = 100;
pilotGain = pow2db(pilotLen);
gain = pilotGain;

snrDbSet = [10,15,20,25,30,35,40];
snrDbSetLen = length(snrDbSet);

b = 4;
errStd = 20;
L = 2^b;
O = 15;
oNb = O*L;
stdPsSet = (0:L-1)/L*2*pi;
stdPhsSet = exp(1j*stdPsSet);

risGeniePhsSetMax = gen_ris_phaseSet(MrisMax, b, errStd);
permtrSetMax = gen_permtr_mtr(MrisMax,L,O);
risEmplyedMax = gen_multi_d(risGeniePhsSetMax,permtrSetMax);

%% channel
Nc = 3; % number of cluster
Nray = 5; % number of rays per cluster
[HbrRSet,hrBtSet] = save_channel(MrisSet,MrSet,Nc,Nray,8);

mseSet = zeros(MrisSetLen,MrSetLen,snrDbSetLen);
crbSet = zeros(MrisSetLen,MrSetLen,snrDbSetLen);
%% crb
for ll = 1 : MrSetLen
    Mr = MrSet(ll);
    for mm = 1 : MrisSetLen
        Mris = MrisSet(mm);
        hrBt = hrBtSet{mm};
        HbrR = HbrRSet{mm,ll};
        Hp = diag(hrBt);
        Hbreve = HbrR*Hp;
        risGeniePhsSet = risGeniePhsSetMax(1:Mris,:);
        permtrSet = permtrSetMax(:,1:Mris,:,:);
        for ss = 1 : snrDbSetLen
            snrDb = snrDbSet(ss);
            noisePow = db2pow(-(snrDb+gain));
            fisher = cal_CRB_of_RIS_cali_correlate_chan(permtrSet,risGeniePhsSet,Hbreve);
            fisherInv = fisher^(-1);
            fisherInvDiag = diag(fisherInv);
            crbSet(mm,ll,ss) = noisePow/2*mean(fisherInvDiag(1:(Mris*(L-1))));
        end
    end
end
%% rmse
nnMonte = 1e2;
for nn = 1 : nnMonte
    for ll = 1 : MrSetLen
        Mr = MrSet(ll);
        onUnitNoise = sqrt(1/2)*(randn(Mr,oNb)+1j*randn(Mr,oNb));
        offUnitNoise = sqrt(1/2)*(randn(Mr,1)+1j*randn(Mr,1))*ones(1,oNb);
        for mm = 1 : MrisSetLen
            Mris = MrisSet(mm);
            if nn == 1 || mod(nn,5) == 0
                fprintf(['\nThe Iteration Number: %d, the number of RIS: %d, at ', datestr(now,"HH:MM"),'\n'], nn, Mris);
            end
            hrBt = hrBtSet{mm};
            HbrR = HbrRSet{mm,ll};
            Hp = diag(hrBt);
            Hbreve = HbrR*Hp;
            risGeniePhsSet = risGeniePhsSetMax(1:Mris,:);
            permtrSet = permtrSetMax(:,1:Mris,:,:);
            risEmplyed = risEmplyedMax(1:Mris,:);
            HeffDeNoise = Hbreve*risEmplyed;
            for ss = 1 : snrDbSetLen
                snrDb = snrDbSet(ss);
                noisePow = db2pow(-(snrDb+gain));
                Heff = HeffDeNoise+sqrt(noisePow)*(onUnitNoise-offUnitNoise);
                risEstPhsSet = cali_Ris_algo_MO_corrlate_chan(Heff,O,L,Mr,Mris,stdPhsSet,permtrSet);
                risEstPhsSetRmAm = remove_ambiguity(risEstPhsSet);
                mse = norm(angle(risEstPhsSetRmAm./risGeniePhsSet),'fro')^2/(Mris*(L-1));
                mseSet(mm,ll,ss) = mseSet(mm,ll,ss)+ mse;
            end
        end
    end
end
mseSetAve = mseSet/nn;
rmseSet = sqrt(mseSetAve)*180/pi;
rcrbSet = sqrt(crbSet)*180/pi;
%% plot
plt1 = {[1,0,0],[0,0,0],[0,0,1];[1 0 1],[0.5,0,0],[0,1,1]};
plt2 = {'o','x','^'};
plt3 = {'-','-.',':'};
figure;
for ll = 1 : MrSetLen
    for mm = 1 : MrisSetLen
        semilogy(snrDbSet,squeeze(rmseSet(mm,ll,:)),plt2{mm},'color',plt1{ll,mm},'linewidth',2,'markersize',8);
        hold on
        semilogy(snrDbSet,squeeze(rcrbSet(mm,ll,:)),plt3{mm},'color',plt1{ll,mm},'linewidth',2,'markersize',8);
    end
end
grid on
xlabel('SNR (dB)');
ylabel('RMSE (^\circ)');
title([num2str(b),'-bit RIS, ', 'Phase Deviations [-',num2str(errStd),'\circ',', ',num2str(errStd),'\circ',']'])
legend('Mr=8,M_{ris}=16,Estimate','Mr=8,M_{ris}=16,CRB',...
    'Mr=8,M_{ris}=32,Estimate','Mr=8,M_{ris}=32,CRB',...,
    'Mr=8,M_{ris}= 64,Estimate','Mr=8,M_{ris}=64,CRB')








