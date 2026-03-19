function [QA,QB] = computeSSLSTAHeatBinaryBT(P, y, T, parameters)


isoparams1 = parameters(1,:)';
isoparams2 = parameters(2,:)';


qsNPA   =   isoparams1(1,1);
b01NPA    = isoparams1(2,1);
delU1NPA  = isoparams1(3,1);
qsLPA   =   isoparams1(4,1);
b01LPA    = isoparams1(5,1);
delU1LPA  = isoparams1(6,1);
kgate  = isoparams1(7,1);
cgate  = isoparams1(8,1);
sval  = isoparams1(9,1);


qsNPB   =   isoparams2(1,1);
b01NPB    = isoparams2(2,1);
delU1NPB  = isoparams2(3,1);
qsLPB   =   isoparams2(4,1);
b01LPB    = isoparams2(5,1);
delU1LPB  = isoparams2(6,1);

PA = y.*P;
PB = (1-y).*P;
yval = ((1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))).^qsLPA) ...
     ./((1+(b01NPA.*PA.*exp(delU1NPA./(8.314.*T))) + (b01NPB.*PB.*exp(delU1NPB./(8.314.*T)))).^qsNPA).* ...
    exp(-(kgate-T.*cgate)./(8.314.*T));

sigmaval = yval.^sval./(1+yval.^sval);
qNPA = qsNPA.*(b01NPA.*PA.*exp(delU1NPA./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))) ;
qLPA = qsLPA.*(b01LPA.*PA.*exp(delU1LPA./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T))));
qNPB = qsNPB.*(b01NPB.*PB.*exp(delU1NPB./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))) ;
qLPB = qsLPB.*(b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T))) + (b01LPB.*PB.*exp(delU1LPB./(8.314.*T))));
% qa = (1-sigmaval).*qNPA+ ...
%     sigmaval.*qLPA;
% 
% qb = (1-sigmaval).*qNPB+ ...
%     sigmaval.*qLPB;
% 


% yval = ((1+(b01LP.*P.*exp(delU1LP./(8.314.*T)))).^qsLP)./((1+(b01NP.*P.*exp(delU1NP./(8.314.*T)))).^qsNP).* ...
%     exp(-(kgate-T.*cgate)./(8.314.*T));
% 
% sigmaval = yval.^sval./(1+yval.^sval);
% 
% qNP = qsNP.*(b01NP.*P.*exp(delU1NP./(8.314.*T)))./(1+(b01NP.*P.*exp(delU1NP./(8.314.*T))));
% qLP = qsLP.*(b01LP.*P.*exp(delU1LP./(8.314.*T)))./(1+(b01LP.*P.*exp(delU1LP./(8.314.*T))));

% qa = (1-sigmaval).*qNP +  sigmaval.*qLP;

QnetA = (qLPA.*delU1LPA - qNPA.*delU1NPA - kgate)./(qLPA-qNPA);
QnetB = (qLPB.*delU1LPB - qNPB.*delU1NPB - kgate)./(qLPB-qNPB);


if qLPB-qNPB == 0
    QnetB = 0;
end
if (qLPA-qNPA) == 0
    QnetA = 0;
end

QnumA = (1-sigmaval).*qNPA./(1+(b01NPA.*PA.*exp(delU1NPA./(8.314.*T)))).*delU1NPA ...
       +  (sigmaval).*qLPA./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T)))).*delU1LPA ...
       + sval.*(1-sigmaval).*sigmaval.*(qLPA-qNPA).^2.*QnetA;
QnumB = (1-sigmaval).*qNPB./(1+(b01NPB.*PB.*exp(delU1NPB./(8.314.*T)))).*delU1NPB ...
       +  (sigmaval).*qLPB./(1+(b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))).*delU1LPB ...
       + sval.*(1-sigmaval).*sigmaval.*(qLPB-qNPB).^2.*QnetB;
QdenA = (1-sigmaval).*qNPA./(1+(b01NPA.*PA.*exp(delU1NPA./(8.314.*T)))) ...
       + (sigmaval).*qLPA./(1+(b01LPA.*PA.*exp(delU1LPA./(8.314.*T)))) ...
       + sval.*(1-sigmaval).*sigmaval.*(qLPA-qNPA).^2;
QdenB = (1-sigmaval).*qNPB./(1+(b01NPB.*PB.*exp(delU1NPB./(8.314.*T)))) ...
       + (sigmaval).*qLPB./(1+(b01LPB.*PB.*exp(delU1LPB./(8.314.*T)))) ...
       + sval.*(1-sigmaval).*sigmaval.*(qLPB-qNPB).^2;

QA = QnumA./QdenA;
QB = QnumB./QdenB;

end