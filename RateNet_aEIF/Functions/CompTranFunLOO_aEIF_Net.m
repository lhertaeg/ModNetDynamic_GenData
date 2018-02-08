function rate = CompTranFunLOO_aEIF_Net(r0,r_fix,rx,ModPar,ConPar,idx)
% This function computes the firing rate(s) of an a network of aEIF neurons
% with a given (predicted/initial) rate for the averaged adaptation,
% whereby one neuron population activity is fixed
%
% r0 = (estimated) rate (vector) that needs to be optimized
% r_fix = fixed rate (the one that is held constant)
% rx = vector of external rates
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% idx = index of neuron population that has fixed activity (idx-nullcline)
% rate = new firing rate

tm = ModPar(:,1)/1000.0;

eps = ConPar.R;
N = ConPar.NN;
J = ConPar.J; 
NTypes = length(N); % plus external excitatory pool

ind = 1:NTypes-1;
ind(idx) = [];

r = zeros(NTypes-1);
r(ind) = r0;
r(idx) = r_fix;

for i=idx
    
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
    rate = 1000.0*sig2/(2*tau*q);
end

end