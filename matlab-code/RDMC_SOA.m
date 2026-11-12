function X = RDMC_SOA(A, k, v_a, v_c,c)
sigma=0.1;
% Boundary separation
dt=0.001; tmax=1.5; t=linspace(dt,tmax,tmax/dt);
% Weighting factor for the automatic drift
% Weighting factor for the controlled drift
w_a=A.*exp(-k.*t);
w_c=1-w_a;
% Simulate time-depedent Wiener process X(t) for an automatic prcessing
dX_a = w_a.*v_a.*dt;
% cumulative sum corresponds to X_a(t)
X_a = cumsum(dX_a);
% Simulate time-depedent Wiener process X(t) for a controlled prcessing
dX_c = w_c.*v_c.*dt;
% cumulative sum corresponds to X_c(t)
X_c = cumsum(dX_c);
% The total drift for a trial
v = w_a.*v_a + w_c.*v_c;
% simulate time-depedent Wiener process X(t) for a single trial 
dX = v.*dt+sigma*sqrt(dt)*randn(1, length(t));
% cumulative sum corresponds to X(t)
X = cumsum(dX)+c;
end