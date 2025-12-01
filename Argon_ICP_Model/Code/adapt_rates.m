function R = adapt_rates(Rin, ps)
% Optional gentle scalers. Pass only fields you want to change.
% Supported ps fields:
%   scale_exc_gm  -> k6
%   scale_exc_gr  -> k7
%   scale_exc_gp  -> k8
%   scale_xfer_to_p -> k13,k14
%   scale_iz_g    -> k2
%   scale_iz_m    -> k3
%   scale_iz_r    -> k4
%   scale_iz_p    -> k5

if nargin < 2 || ~isstruct(ps), ps = struct(); end
R = Rin;

s_gm = getp(ps,'scale_exc_gm',   1.00);
s_gr = getp(ps,'scale_exc_gr',   1.00);
s_gp = getp(ps,'scale_exc_gp',   1.00);
s_xp = getp(ps,'scale_xfer_to_p',1.00);

s_iz_g = getp(ps,'scale_iz_g',   1.00);
s_iz_m = getp(ps,'scale_iz_m',   1.00);
s_iz_r = getp(ps,'scale_iz_r',   1.00);
s_iz_p = getp(ps,'scale_iz_p',   1.00);

if isfield(R,'k6'),  R.k6  = R.k6  * s_gm; end
if isfield(R,'k7'),  R.k7  = R.k7  * s_gr; end
if isfield(R,'k8'),  R.k8  = R.k8  * s_gp; end
if isfield(R,'k13'), R.k13 = R.k13 * s_xp; end
if isfield(R,'k14'), R.k14 = R.k14 * s_xp; end

if isfield(R,'k2'),  R.k2  = R.k2  * s_iz_g; end
if isfield(R,'k3'),  R.k3  = R.k3  * s_iz_m; end
if isfield(R,'k4'),  R.k4  = R.k4  * s_iz_r; end
if isfield(R,'k5'),  R.k5  = R.k5  * s_iz_p; end
end

function val = getp(S, name, def)
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    val = S.(name);
else
    val = def;
end
end