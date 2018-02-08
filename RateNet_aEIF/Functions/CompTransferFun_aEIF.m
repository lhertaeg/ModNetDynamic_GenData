function rate = CompTransferFun_aEIF(r,ModPar,InpPar)
% This function computes the firing rate of an aEIF neuron in the presence
% of Gaussian white noise with a given (predicted/initial) rate for the
% averaged adaptation
%
% r = (estimated) rate to compute averaged adaptation is computed
% ModPar = vector of neuron model parameters
% InpPar = vector of input mean and standard deviation
% rate = firing rate of the neuron

% rename properly
CModPar = num2cell(ModPar);
CInpPar = num2cell(InpPar);
[tm,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
[mu,sig] = CInpPar{:};
lb = min(EL,Vr) - 20.0;

% solve/compute integral
S = @(x) (2.0/sig^2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r/gL).*x + sf^2*exp((x-VT)./sf));
q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
rate = 1000.0*sig^2/(2*tm*q);

end

