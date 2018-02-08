%% Computing neuron (population) properties in the presence of Gaussian white noise
%
%
%% 
% The functions needed to compute the properties are in
addpath('Functions\')

%% Neuron model and parameters
%
% The neurons are modeled by the adaptive exponential I&F model (aEIF) that
% obeys the following differential equations:
%
% $$ C\cdot\frac{dV}{dt}=-g_L \cdot (V-E_L) + g_L \cdot\Delta_T\cdot e^{(\frac{V-V_T}{\Delta_T})}+I-w $$
%
% $$ \tau_w\cdot\frac{dw}{dt}=-w $$
% 
% When $V>V_\mathrm{up}$, $V\rightarrow V_r$ and $w\rightarrow w_r=w+b$.
%
% Please note that we only allow for spike-triggered adaptation (denoted by
% _b_). In the following, we compute the mean firing rate, the variance of
% the interspike intervals & the mean and variance of the membrane
% potential in the presence of Gaussian white noise.
% We set the neuron parameters to be the following:
%
ModPar = [10.0,5.0,-70.0,2.0,-50.0,-70.0,-20.0,5.0,100.0]; % [tm,gL,EL,sf,VT,Vr,Vup,b,tw]

%% Mean firing rate and variance of interspike intervals
%
% In the presence of Gaussian white noise (input $\in N(\mu,\sigma)$), the firing rate $\nu_0$ is given/approximated by
%
% $$ \nu_0=\left(\frac{2 \tau_m}{\sigma^2}\int\limits^{V_\mathrm{up}}_{-\infty}dV\int\limits^{V_\mathrm{up}}_{max(V,V_r)} e^{-\frac{2}{\sigma^2} \int\limits^u_V \left(f(x) + \mu - b\nu_0\tau_w/g_L\right) dx}du\right)^{-1}. $$
%
% with $f(V)=-(V-E_L)+\Delta_T\cdot e^{\frac{V-V_T}{\Delta_T}}$. Furthermore, the variance of the interspike-intervals is given by
%
% $$ \sigma^2_{ISI}=\frac{8\tau^2}{\sigma_0^4}\int\limits_{V_r}^{\infty}dV \
% e^{-\frac{2}{\sigma_0^2}F(V)}\int\limits_{-\infty}^{V}du \
% e^{-\frac{2}{\sigma_0^2}F(u)}\left[\int\limits_{-\infty}^{u}dy \
% e^{\frac{2}{\sigma_0^2}F(y)}\  \right]^2 $$
%
% where $F(x)=E\cdot x - \frac{x^2}{2} + \Delta_T e^{(x-V_T)/\Delta_T}$ and
% $E=E_L+\mu-\frac{b\tau_w\nu}{gL}$.
%
% We characterize the firing rate and the variance of the interspike
% intervals as a function of the external rate. Therefore, we define the
% external stimuli to be tested:
%
rext = linspace(4.0,10.0,40);
%%%
% The mean and the standard deviation of the input are a function of the
% external stimulation defined by the external rate. We define the number
% of connections, N, and the postsynaptic potentials by
N = 100;
J = 2;
%%%
% The firing rate and the var(ISI) as a function of the external rate:
%
tm = ModPar(1);
rate = zeros(1,length(rext));
VarISI = zeros(1,length(rext));

for i=1:length(rext)
    mu = (tm/1000.0)*N*rext(i)*J;
    sig = sqrt((tm/1000.0)*N*rext(i)*J.^2);
    InpPar = [mu,sig];
    
    rate(i) = CompRate_aEIF(ModPar,InpPar,r,1);
    VarISI(i) = CompVarISI_aEIF(rate(i),ModPar,InpPar);
end

subplot(1,2,1); plot(rext,rate,'LineWidth',2); 
xlabel('r_{ext} (Hz)'); ylabel('rate (Hz)')
subplot(1,2,2); plot(rext,sqrt(VarISI),'LineWidth',2); 
set(gca, 'YScale', 'log')%ylim([0,1000])
xlabel('r_{ext} (Hz)'); ylabel('\sigma_{ISI} (s)')
set(gcf,'units','centimeters','position',[1,1,18,7])
snapnow; close

%% Mean and variance of the membrane potential
%
% The probability distribution of the membrane potential is given by
%
% $$ P(V)=\frac{2\nu_0 \tau_m}{\sigma^2}\int\limits^{V_\mathrm{up}}_{max(V,V_r)} e^{-\frac{2}{\sigma^2} \int\limits^u_V \left(f(x) + \mu - b\nu_0\tau_w/g_L\right) dx}du. $$
%
% Consequently, first and second moment of the membrane potential are given by
%
% $$ \langle V \rangle  = \frac{2\nu_0 \tau_m}{\sigma^2}\int\limits^{V_\mathrm{up}}_{-\infty} V dV\int\limits^{V_\mathrm{up}}_{max(V,V_r)}
% e^{-\frac{2}{\sigma^2} \int\limits^u_V \left(f(x) + \mu -
% b\nu_0\tau_w/g_L\right) dx}du. $$
%
% $$ \langle V^2 \rangle  = \frac{2\nu_0 \tau_m}{\sigma^2}\int\limits^{V_\mathrm{up}}_{-\infty} V^2 dV\int\limits^{V_\mathrm{up}}_{max(V,V_r)}
% e^{-\frac{2}{\sigma^2} \int\limits^u_V \left(f(x) + \mu -
% b\nu_0\tau_w/g_L\right) dx}du. $$
%
% (Please note that spikes should be discarded for the computation.
% Therefore the upper bound should be corrected - if necessary)
%
%
Vm = zeros(1,length(rext));
Vv = zeros(1,length(rext));

for i=1:length(rext)
    mu = (tm/1000.0)*N*rext(i)*J;
    sig = sqrt((tm/1000.0)*N*rext(i)*J.^2);
    InpPar = [mu,sig];
    
    [Vm(i),Vv(i)] = CompMemPotStat_aEIF(rate(i),ModPar,InpPar);
end

subplot(1,2,1); plot(rext,Vm,'LineWidth',2); 
xlabel('r_{ext} (Hz)'); ylabel('mean membrane potential (mV)')
subplot(1,2,2); plot(rext,sqrt(Vv),'LineWidth',2); 
xlabel('r_{ext} (Hz)'); ylabel('SD of membrane potential (mV)')
set(gcf,'units','centimeters','position',[1,1,18,7])
snapnow; close

%%
% Note: Check if we need to set the upper threshold for calculating mean and SD of V to V'<Vup to cut out spikes ...
%