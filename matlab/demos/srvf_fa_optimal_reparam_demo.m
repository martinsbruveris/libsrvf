load -ASCII ages.mat;
load -ASCII boys_rates.mat;

T1 = ages;
T2 = ages;

F1 = boys_rates(1,:);
F2 = boys_rates(2,:);

% srvf_fa_optimal_reparam assumes that the functions are 
% scaled to unit arclength, and have constant-speed parametrizations
F1 = F1 ./ plf_arclength(F1, T1);
F2 = F2 ./ plf_arclength(F2, T2);

[Gn1 TGn1] = plf_constant_speed_reparam(F1, T1);
[Gn2 TGn2] = plf_constant_speed_reparam(F2, T2);

[F1n TF1n] = plf_compose(F1, T1, Gn1, TGn1);
[F2n TF2n] = plf_compose(F2, T2, Gn2, TGn2);

Q1 = plf_to_srvf(F1n, TF1n);
Q2 = plf_to_srvf(F2n, TF2n);

% srvf_fa_optimal_reparam also assumes that the SRVF values
% alternate between 1 and -1 on adjacent intervals
[Q1n, TQ1n] = srvf_make_alternating(Q1, TF1n);
[Q2n, TQ2n] = srvf_make_alternating(Q2, TF2n);

[G1, TG1, G2, TG2] = srvf_fa_optimal_reparam(Q1n, TQ1n, Q2n, TQ2n);

[F1r, TF1r] = plf_compose(F1n, TF1n, G1, TG1);
[F2r, TF2r] = plf_compose(F2n, TF2n, G2, TG2);

figure();
clf();
plot(TF1r, F1r, 'b', TF2r, F2r, 'r');
