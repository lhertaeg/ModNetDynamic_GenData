function [u,v] = CompVecField_aEIF_Net(ra,rb,rx,ModPar,ConPar)
% This function computes the vector field of an aEIF network around a FP
%
% x = 
% rx = vector of external rates
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)

tm = ModPar(:,1);
u = zeros(length(rb),length(ra));
v = zeros(length(rb),length(ra));

for j=1:length(rb)
    for i=1:length(ra)
        r = [ra(i),rb(j)];
        rates = CompTransferFun_aEIF_Net(r,rx,ModPar,ConPar);
        u(j,i) = (rates(1) - ra(i))/tm(1);
        v(j,i) = (rates(2) - rb(j))/tm(2);     
    end
end
