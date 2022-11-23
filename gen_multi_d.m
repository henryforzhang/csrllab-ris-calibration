function risEmplyed = gen_multi_d(risPhsSet,permtrSet)
[O,Mris,Nb,~] = size(permtrSet);
oNb=O*Nb;
risEmplyed = zeros(Mris,oNb);
for oo = 1:O
    permtr = squeeze(permtrSet(oo,:,:,:));
    tmp = zeros(Mris,Nb);
    for mm = 1:Mris
       tmp(mm,:) =  squeeze(permtr(mm,:,:))*risPhsSet(mm,:).';
    end
    risEmplyed(:,(oo-1)*Nb+(1:Nb)) = tmp;
end
end