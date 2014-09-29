function [Jstat, g_t, g_T] = OLS_J(param_vec, data, W)

alpha = param_vec(1); beta = param_vec(2);
y = data(:,1); x = data(:,2);
T = size(y,1);

u_t = nan(T,1);
g_t = nan(T,2);
for t = 1:T 
      u_t(t) = y(t) - alpha - beta*x(t);
      g_t(t,:) = kron(u_t(t),[1 x(t)]);
end

g_T = mean(g_t);

Jstat = g_T*W*g_T';