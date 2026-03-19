function [qa, qb] = computeSSLSTALoadingBinary(P,T,y,qsNPA,qsLPA,b01NPA,b01LPA,delU1NPA,delU1LPA,qsNPB,qsLPB,b01NPB,b01LPB,delU1NPB,delU1LPB,kgate,cgate,sval)
PA = y.*P;
PB = (1-y).*P;
yval = ((1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))).^qsLPA) ...
     ./((1+(b01NPA.*PA.*exp(delU1NPA./(8.314.*T))) + (b01NPB.*PB.*exp(delU1NPB./(8.314.*T)))).^qsNPA).* ...
    exp(-(kgate-T.*cgate)./(8.314.*T));

sigmaval = yval.^sval./(1+yval.^sval);

qa = (1-sigmaval).*qsNPA.*(b01NPA.*PA.*exp(delU1NPA./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))) + ...
    sigmaval.*qsLPA.*(b01LPA.*PA.*exp(delU1LPA./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T))));
qb = (1-sigmaval).*qsNPB.*(b01NPB.*PB.*exp(delU1NPB./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))) + ...
    sigmaval.*qsLPB.*(b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T))));

end