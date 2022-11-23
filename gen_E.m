function E = gen_E(HbreveEst,permtrSet)
Mr = size(HbreveEst,1);
[O,Mris,Nb,~] = size(permtrSet);
mrNb = Mr*Nb;
E = zeros(Mr*O*Nb,Mris);
for oo = 1 : O
    for mm = 1 : Mris
        E((oo-1)*mrNb+(1:mrNb),(mm-1)*Nb+(1:Nb)) = kron(HbreveEst(:,mm),squeeze(permtrSet(oo,mm,:,:)));
    end
end
end