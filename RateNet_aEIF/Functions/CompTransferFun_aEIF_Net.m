function rates = CompTransferFun_aEIF_Net(r,rx,ModPar,ConPar)
% This function computes the firing rates of an a network of aEIF neurons
% with a given (predicted/initial) rate for the averaged adaptation
%
% r = (estimated) rate to compute averaged adaptation and input parameters
% rx = vector of external rates
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% rates = vector of firing rates

tm = ModPar(:,1)/1000.0;

eps = ConPar.R;
N = ConPar.NN;
J = ConPar.J; 
NTypes = length(N); % plus external excitatory pool
rates = zeros(1,NTypes-1);

for i=1:NTypes-1
    
    % influence of external pool
    Iex = eps(i,NTypes)*N(NTypes)*J(i,NTypes)*rx(i)*tm(i);
    mu = Iex;
    sig2 = J(i,NTypes)*Iex;
    
    for j=1:NTypes-1
        mu = mu + tm(i)*eps(i,j)*N(j)*r(j)*J(i,j);
        sig2 = sig2 + tm(i)*eps(i,j)*N(j)*r(j)*J(i,j).^2;
    end 
    
    % compute transfer function  
    CModPar = num2cell(ModPar(i,:));
    [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
    lb = min(EL,Vr) - 20.0;

    S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(i)/gL).*x + sf^2*exp((x-VT)./sf));
    q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
    rates(i) = 1000.0*sig2/(2*tau*q);
end

end

