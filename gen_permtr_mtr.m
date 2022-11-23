function permtrSet = gen_permtr_mtr(Mris,Nb,O)
permtrSet = zeros(O,Mris,Nb,Nb);
for oo = 1 : O
    for mm = 1 : Mris
        permtrSet(oo,mm,:,:) = random_permut_mtr(Nb);
    end
end

end