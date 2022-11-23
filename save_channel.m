function [HbrRSet,hrBtSet] = save_channel(MrisSet,MrSet,Nc,Nray,p)
HbrRSet=channel_generation_ura(MrisSet,MrSet,Nc,Nray,p);
hrBtSet=channel_generation_ula(MrisSet,Nc,Nray,p);
end