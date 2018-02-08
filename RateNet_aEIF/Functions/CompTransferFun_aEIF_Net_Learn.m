function rates = CompTransferFun_aEIF_Net_Learn(r,rx,ModPar,ConPar,LrnPar)
% This function computes the firing rates of a network of aEIF neurons
% with given (predicted/initial) rates after learning
%
% r = (N+2)-vector of rates (+, sel, 0 and all remaining populations)
% rx = N-vector of external rates
% ModPar = matrix of neuron model parameters (rows = N neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% LrnPar = structure containing all parameters related to learning
% rates = (N+2)-vector of firing rates

tm = ModPar(:,1)/1000.0;

f = LrnPar.f;
p = LrnPar.p;
m = LrnPar.m; % more (g_+)
l = 1 - f*(m-1)/(1-f); % less (g_-)
SDJ = LrnPar.SDJ;  
L = (1 + SDJ^2);

eps = ConPar.R;
N = ConPar.NN;
J = ConPar.J; 
NTypes = length(N); % plus external excitatory pool
rates = zeros(1,NTypes+1); 
KEE = tm(1)*eps(1,1)*N(1)*J(1,1);

% excitatory pools after learning
for i=1:3 % exc pools (+, sel, 0)
    
    % influence of external pool
    mu = eps(1,NTypes)*N(NTypes)*J(1,NTypes)*rx(1)*tm(1);
    sig2 = J(1,NTypes)*mu*L;
    
    % influence of remaining populations
    for j=2:NTypes-1
        mu = mu + tm(1)*eps(1,j)*N(j)*r(2+j)*J(1,j);
        sig2 = sig2 + L*tm(1)*eps(1,j)*N(j)*r(2+j)*J(1,j).^2;
    end
    
    % influence of exc pools
    if i==1 % +
        mu = mu + (((p-2)*l + m)*f*r(1) + f*l*r(2) + (1-p*f)*l*r(3))*KEE;
        sig2 = sig2 + L*(((p-2)*l^2 + m^2)*f*r(1) + f*l^2*r(2) + (1-p*f)*l^2*r(3))*KEE*J(1,1);
    elseif i==2 % sel
        mu = mu + ((p-1)*f*l*r(1) + f*m*r(2) + (1-p*f)*l*r(3))*KEE;
        sig2 = sig2 + L*((p-1)*f*l^2*r(1) + f*m^2*r(2) + (1-p*f)*l^2*r(3))*KEE*J(1,1);
    elseif i==3 % 0
        mu = mu + ((p-1)*f*r(1) + f*r(2) + (1-p*f)*r(3))*KEE;
        sig2 = sig2 + L*((p-1)*f*r(1) + f*r(2) + (1-p*f)*r(3))*KEE*J(1,1);
    end
     
    % compute transfer function  
    CModPar = num2cell(ModPar(1,:));
    [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
    lb = min(EL,Vr) - 20.0;

    S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(i)/gL).*x + sf^2*exp((x-VT)./sf));
    q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
    rates(i) = 1000.0*sig2/(2*tau*q);

end
rE = f*(r(2) + (p-1)*r(1)) + (1-p*f)*r(3);

% remaining populations
for i=2:NTypes-1
    
    % influence of external pool
    mu = eps(i,NTypes)*N(NTypes)*J(i,NTypes)*rx(i)*tm(i);
    sig2 = J(i,NTypes)*mu*L;
    
    % influence of excitatory pools
    mu = mu + tm(i)*eps(i,1)*N(1)*rE*J(i,1);
    sig2 = sig2 + L*tm(i)*eps(i,1)*N(1)*rE*J(i,1).^2;
    
    % influence of remaining population
    for j=2:NTypes-1
        mu = mu + tm(i)*eps(i,j)*N(j)*r(2+j)*J(i,j);
        sig2 = sig2 + L*tm(i)*eps(i,j)*N(j)*r(2+j)*J(i,j).^2;
    end 
    
    % compute transfer function  
    CModPar = num2cell(ModPar(i,:));
    [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
    lb = min(EL,Vr) - 20.0;

    S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(2+i)/gL).*x + sf^2*exp((x-VT)./sf));
    q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
    rates(2+i) = 1000.0*sig2/(2*tau*q);
end

end

