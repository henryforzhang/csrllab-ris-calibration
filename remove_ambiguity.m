function risEstPhsSetRmAm = remove_ambiguity(risEstPhsSet)
%% remove ambiguity
[Mris,Nb] = size(risEstPhsSet);
risEstPhsSetRmAm = zeros(Mris,Nb);
for mm = 1 : Mris
    risEstPhsSetRmAm(mm,:) = risEstPhsSet(mm,:)/risEstPhsSet(mm,1);
end
end