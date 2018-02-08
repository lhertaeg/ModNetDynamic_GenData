function [Vm,Vv] = CompMemPotStat_aEIF(r,ModPar,InpPar)
% This functions computes the mean and the variance of the membrane
% potential of an aEIF neuron in the presence of Gaussian white noise
%
% r = mean firing rate
% ModPar = vector of neuron model parameters
% InpPar = vector of input mean and standard deviation
% Vm = mean membrane potential
% Vv = variance of membrane potential

% rename properly
CModPar = num2cell(ModPar);
CInpPar = num2cell(InpPar);
[tm,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
[mu,sig] = CInpPar{:};
lb = min(EL,Vr) - 20.0;
%ATol=1e-10;
%RTol=1e-10;

% solve/compute integral
S = @(x) (2.0/sig^2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r/gL).*x + sf^2*exp((x-VT)./sf));
Vm = 2.*(tm./1000.0).*r.*integral2(@(V,u) V.*exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup)./sig^2;
V2 = 2.*(tm./1000.0).*r.*integral2(@(V,u) V.^2.*exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup)./sig^2;
Vv = V2 - Vm.^2;

% if rate is very small, numerical issue occur, to avoid this apply
% rectification
if Vv<0
    Vv=0;
end

end

