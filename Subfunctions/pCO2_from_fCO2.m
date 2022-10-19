function [pCO2] = pCO2_from_fCO2(fCO2, T)
% taken from CO2SYS / Dickson 2007
% 2018_10_24 SMB

TempK = T + 273.15;

RGasConstant = 83.1451;  % ml bar-1 K-1 mol-1, DOEv2
RT  = RGasConstant.*TempK;


Delta = (57.7 - 0.118.*TempK);
b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;

P1atm = 1.01325; % in bar
FugFac = exp((b + 2.*Delta).*P1atm./RT);


pCO2 = fCO2./FugFac;
end