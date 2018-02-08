function [r_spon,r_memo,FPStab] = CompBifDiagram(ModPar,ConPar,LrnPar,r0,rx,g,D)

warning('off','all')

r_spon = zeros(length(g),4);
r_memo = zeros(length(g),4);
FPStab = zeros(length(g),2);

% spontaneous states
for i=1:length(g)
    LrnPar.m = g(i); 
    r_spon(i,:) = CompRate_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r0,rx,0);
    FPStab(i,1) = CompEigVal_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r_spon(i,:),rx,0);
end

% memory states
e = 1.0;
for i=1:length(g)
    LrnPar.m = g(i); 
    % decide for initial condition to be tested
    if (i==1 || e<0.5)
        r0 = r_spon(i,:) + [0,D,0,0];
    else
        r0 = r_memo(i-1,:) + [0,10.0,0,0];
    end
    % compute potential memory state and error
    r_memo(i,:) = CompRate_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r0,rx,0);
    FPStab(i,2) = CompEigVal_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r_memo(i,:),rx,0);
    if FPStab(i,2)>=0
        r_memo(i,:) = CompRate_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r0,rx,1);
        FPStab(i,2) = CompEigVal_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r_memo(i,:),rx,0);
    end
    e = abs(r_memo(i,2)-r_memo(i,3));    
end