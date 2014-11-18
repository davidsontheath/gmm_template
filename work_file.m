clear all; close all; clc

%simulate data and run regular OLS
T = 100;
x = randn(T,1);
eps_raw = trnd(10,T-1,1);
y = 2*x + [eps_raw;0] + exp(x).*[0;eps_raw];
figure;plot(x,y, '.'); title('scatterplot Y vs X')
OLS_results = regstats(y,x);

%just evaluate our J function at some parameters,
%make sure it's working
testparams = [0;1];
Jstat_test = OLS_J(testparams, [y,x], eye(2))


%define function call that only takes params
%using the identity weighting matrix
OLS_J1 = @(param_vec) OLS_J(param_vec, [y,x], eye(2));
Jstat_test2 = OLS_J1([0;1])


%okay, time to search for best params
%initial grid search
bestgridJ = 1e10; bestgridtheta = [NaN; NaN];
for gridalpha = -10:1:10
    for gridbeta = -10:1:10
        Jstat = OLS_J1([gridalpha;gridbeta]);
        if Jstat < bestgridJ
            bestgridJ = Jstat;
            bestgridtheta = [gridalpha;gridbeta];
        end
    end %i
end %j

%do a proper search starting from best grid point
theta_hat = fminunc(OLS_J1, bestgridtheta);
disp('OLS coeffs   1st stage coeffs')
disp([OLS_results.beta, theta_hat]);


%%%%%%%%%%%%%%
%second stage
%%%%%%%%%%%%%%
%Evaluate model at first stage thetahat
[Jstat, g_t, g_T] = OLS_J1(theta_hat);

%look at the errors
figure; plot(g_t); title('Time series of the moments g_t')
figure; plot(g_t(1:end-1,1),g_t(2:end,1),'+'); title('g1_t vs g1_{t-1}')
figure; plot(g_t(1:end-1,2),g_t(2:end,2),'o'); title('g2_t vs g2_{t-1}')

%estimate Acovg
Acovg = g_t.'*g_t/T;
num_lags = 1;
for n = 1:num_lags
    NWweight = 1 - n/(num_lags+1);
    lag_cov = g_t(1+n:end,:).'*g_t(1:end-n,:)/T;
    Acovg = Acovg + NWweight*(lag_cov+lag_cov');
end

%optimal weight matrix, and function call that uses that
W2 = inv(Acovg);
OLS_J2 = @(param_vec) OLS_J(param_vec, [y,x], W2);

%2nd stage estimate, efficient under heteroskedasticity and autocorrelation
theta_hat2 = fminunc(OLS_J2, theta_hat);

disp('OLS coeffs   2nd stage coeffs')
disp([OLS_results.beta, theta_hat2]);

%standard errors, robust to heteroskedasticity and autocorrelatoin
%for this we need gprime the change in g_T w.r.t. theta
%delta method, or "standard errors for anything"...

stepsize = 1e-10;
[ans, ans, g_T] = OLS_J2(theta_hat2);
for i = 1:2
    theta_hat2_fd = theta_hat2;
    theta_hat2_fd(i) = theta_hat2(i)+stepsize;
    [ans, ans, g_T_fd] = OLS_J2(theta_hat2_fd);
    dgT(:,i) = (g_T_fd - g_T)'/stepsize;
end

thetahat2_SE = sqrt(diag(inv(dgT'*W2*dgT))/T);

%compare GMM vs OLS standard errors
disp('OLS SEs   GMM SEs')
disp([OLS_results.tstat.se thetahat2_SE]);


