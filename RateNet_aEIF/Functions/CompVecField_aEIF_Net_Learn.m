function [u,v] = CompVecField_aEIF_Net_Learn(ra,rb,rx,FP,a,b,ModPar,ConPar,LrnPar)
% This function computes the vector field of an aEIF network around a FP 
% after learning (that is, several excitatory pools have been generated and
% respond to specific stimuli) 
%
% ra/rb = vectors that contain points at which the vector field is computed
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% LrnPar = structure containing all parameters related to learning
% rx = vector of external rates

tm = ModPar(:,1);
u = zeros(length(rb),length(ra));
v = zeros(length(rb),length(ra));

for j=1:length(rb)
    for i=1:length(ra)
        r = FP;
        r(a) = ra(i);
        r(b) = rb(j);
        rates = CompTransferFun_aEIF_Net_Learn(r,rx,ModPar,ConPar,LrnPar);
        u(j,i) = (rates(a) - ra(i))/tm(1);
        v(j,i) = (rates(b) - rb(j))/tm(2);     
    end
end
