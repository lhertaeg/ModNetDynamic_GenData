function NullArray = CompNullclines_aEIF_Net(ModPar,ConPar,r0,rx,Values)
% This function computes the nullclines of an aEIF network 
%
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% r0 = vector of steady-state activity
% rx = vector of external rates
% Values = matrix with fixed population activities (row = consecutive
%          fixed values for population activity, column = populations)
% NullArray = 3D array containing nullclines
%             first index = consecutive values for activity
%             second index = populations
%             third index = index of population that is held constant

N = size(ModPar,1);
M = size(Values,2);
NullArray = zeros(M,N,N);
r_FP = r0;

for i=1:N % i-th nullcline
    r0 = r_FP;
    ind = 1:N;
    ind(i) = [];
    for j=1:M % j-th data point
        options=optimset('TolFun',1e-5,'TolX',1e-5,'Display','off');
        MSE = @(x) sum((CompTranFunLOO_aEIF_Net(x,Values(i,j),rx,ModPar,ConPar,i)-Values(i,j)).^2);
        [rates,err] = fminsearch(MSE,r0(ind),options);
        NullArray(j,ind,i) = rates;
        NullArray(j,i,i) = Values(i,j);
        if err>N*1e-8
            warning('Optimization error > Ne-10')
        end
    end
end

end

