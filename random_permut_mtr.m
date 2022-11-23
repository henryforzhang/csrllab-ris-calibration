function permtr = random_permut_mtr(Nb)
permtr = zeros(Nb);
leftSet = 1:Nb;
for ii = 1 : Nb
    randSeed = randi([1 Nb-ii+1]);
    idxTmp = leftSet(randSeed);
    leftSet = setdiff(leftSet,idxTmp);
    permtr(ii,idxTmp) = 1;
end
end