load -ASCII ages.mat;
load -ASCII boys_rates.mat;

T1 = ages;
T2 = ages;

F1 = boys_rates(1,:);
F2 = boys_rates(3,:);

%% Calculate distance
% Distance before alignment
Q1 = plf_to_srvf(F1, T1);
Q2 = plf_to_srvf(F2, T2);
dist1 = srvf_l2distance(Q1, T1, Q2, T2);

% Distance after alignment
dist2 = plf_reparam_distance(F1, T1, F2, T2, 'klr');

fprintf('The distance between F1 and F2 before and after alignment.\n');
fprintf('Before alignment, dist = %.4f\n', dist1);
fprintf('After alignment, dist = %.4f\n\n', dist2); 

%% Find optimal matching
[G1, TG1, G2, TG2, segPts] = srvf_klr_optimal_reparam(Q1, T1, Q2, T2);

[F1r, TF1r] = plf_compose(F1, T1, G1, TG1);
[F2r, TF2r] = plf_compose(F2, T2, G2, TG2);

figure;
plot(T1, F1, '-bo', T2, F2, '-ro');
title('Functions before alignment.');
 
figure;
plot(TF1r, F1r, '-bo', TF2r, F2r, '-ro');
title('Functions after alignment.');

figure;
plot_gamma(G1, T1, G2, T2, '-ro');
hold on;
plot(segPts(1,:), segPts(2,:), 'ko');
hold off;
title('Optimal reparametrizations.');

%% Show symmetry of distance and matching
% Compute distance from F2 to F1
dist3 = plf_reparam_distance(F2, T2, F1, T1, 'klr');

fprintf('Symmetry of the distance.\n');
fprintf('dist(F1, F2) = %.8f\n', dist2);
fprintf('dist(F2, F1) = %.8f\n', dist3);
fprintf('Difference = %.2e\n', abs(dist2 - dist3));