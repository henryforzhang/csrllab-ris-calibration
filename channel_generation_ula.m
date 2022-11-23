function [H]=channel_generation_ula(NrSet,Nc,Nray,p0)
% ***************************************
%  generate URA channel
%  author - Wei Zhang
%  input:
%            NrSet: receiver antenna number
%            Nc: cluster number
%            Nray: ray number
%            p0: line number for transmitter
%  output-H: 
%             URA Channel
%copyright - CSRL@Fudan,2021/01/12
%  ************************************
angle_sigma = 20/180*pi; %standard deviation of the angles in azimuth and elevation both of Rx and Tx
Np = Nc*Nray;
gamma = sqrt(1/Np); %normalization factor
sigma = 1; %according to the normalization condition of the H
AoD = zeros(2,Np);

for c = 1:Nc
    AoD_m = unifrnd(0,2*pi,1,2);
    
    AoD(1,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoD_m(1),angle_sigma);
    AoD(2,(c-1)*Nray+1:Nray*c) = laprnd(1,Nray,AoD_m(2),angle_sigma);
end

NrSetLen = length(NrSet);
H = cell(NrSetLen,1);
for nn = 1 : NrSetLen
    Nr = NrSet(nn);
    Htmp = zeros(Nr,1);
    for j = 1:Np
        Ar = array_response_UPA(AoD(1,j),AoD(2,j),Nr,p0); %UPA array response
        alpha = normrnd(0,sqrt(sigma/2)) + normrnd(0,sqrt(sigma/2))*sqrt(-1);
        Htmp = Htmp + alpha * Ar;
    end
    H{nn} = gamma * Htmp;
end