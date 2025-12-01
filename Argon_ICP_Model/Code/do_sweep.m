function Y = do_sweep(PRFvec, y0, params, tspan, opts)
% Run a PRF sweep with warm starts; returns steady states row-wise.
N = numel(PRFvec);
Y = nan(N, numel(y0));
y = y0;
for kk = 1:N
    params.PRF = PRFvec(kk);                %#ok<NASGU>  (captured by handle)
    sol = ode15s(@(t,z) rhs_global(t,z,params,[]), tspan, y, opts);
    y   = sol.y(:,end);
    Y(kk,:) = y.';
end
end
