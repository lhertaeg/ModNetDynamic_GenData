function E = CompEnergy_aEIF_Net(rates,ModPar,ConPar,r0,rx,id)
% This function computes the enery landscape of an aEIF network (along one
% dimension, whereby all others are kept sonstant - at FP)
%
% rates = vector of rates at which energy landscape is computed
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% r0 = vector of steady-state activity
% rx = vector of external rates
% id = cell type or population number that is varied
% E = "energy" landscape computed at values given in rates
    
tm = ModPar(:,1)/1000.0;

eps = ConPar.R;
N = ConPar.NN;
J = ConPar.J; 
NTypes = length(N); % plus external excitatory pool
r = zeros(1,length(rates));
E = zeros(1,length(rates));
list = 1:NTypes-1;
list(id) =[];

% influence of external pool
Iex = eps(id,NTypes)*N(NTypes)*J(id,NTypes)*rx(id)*tm(id);
m = Iex;
s2 = J(id,NTypes)*Iex;

% influence of all other populations except id
for j=list
    m = m + tm(id)*eps(id,j)*N(j)*r0(j)*J(id,j);
    s2 = s2 + tm(id)*eps(id,j)*N(j)*r0(j)*J(id,j).^2;
end

% compute potential landscape
for i=1:length(rates)
    
    mu = m + tm(id)*eps(id,id)*N(id)*rates(i)*J(id,id);
    sig2 = s2 + tm(id)*eps(id,id)*N(id)*rates(i)*J(id,id).^2;
    
    CModPar = num2cell(ModPar(id,:));
    [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
    lb = min(EL,Vr) - 20.0;
     
    S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*rates(i)/gL).*x + sf^2*exp((x-VT)./sf));
    q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
    r(i) = 1000.0*sig2/(2*tau*q); 
    
    if i>1
        F = @(z) interp1(rates(1:i),r(1:i),z,'spline');
        E(i) = -integral(@(x) (F(x)-x),rates(1),rates(i))/tm(id);
    else
        E(1) = 0;
    end
    
end

end
