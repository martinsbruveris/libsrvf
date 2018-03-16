load -ASCII ages.mat;
load -ASCII boys_rates.mat;

nfuncs = size(boys_rates, 1);

% Original piecewise-linear functions.  Values in Fs{i}, parameters in Ts{i}.
Fs = cell(1, nfuncs);
Ts = cell(1, nfuncs);

% PLFs scaled to unit arclength and reparametrized to constant speed
Fsn = cell(1, nfuncs);
Tsn = cell(1, nfuncs);

% The reparametrizations that were applied to the scaled-down Fs to yield Fsn
Gsn = cell(1, nfuncs);
TGsn = cell(1, nfuncs);


% Compute the normalized PLFs
for i=1:nfuncs
  Ts{i} = ages;
  Fs{i} = boys_rates(i,:) ./ plf_arclength(boys_rates(i,:), ages);
  [Gsn{i} TGsn{i}] = plf_constant_speed_reparam(Fs{i}, Ts{i});
  [Fsn{i} Tsn{i}] = plf_compose(Fs{i}, Ts{i}, Gsn{i}, TGsn{i});
end


% Fm is the function to which all of the others are aligned
Fm = Fs{1};
Tm = Ts{1};

% The normalized Fm and corresponding reparametrization
Fmn = Fsn{1};
Tmn = Tsn{1};
Gmn = Gsn{1};
TGmn = TGsn{1};


% The SRVFs for the normalized PLFs
Qs = cell(1, nfuncs);
TQs = cell(1, nfuncs);
for i=1:nfuncs
  Qs{i} = plf_to_srvf(Fsn{i}, Tsn{i});
  TQs{i} = Tsn{i};
end
Qm = plf_to_srvf(Fmn, Tmn);
TQm = Tmn;


% All SRVFs are assumed to alternate between 1 and -1 on adjacent intervals
Qsn = cell(1, nfuncs);
TQsn = cell(1, nfuncs);
for i=1:nfuncs
  [Qsn{i}, TQsn{i}] = srvf_make_alternating(Qs{i}, TQs{i});
end
[Qmn, TQmn] = srvf_make_alternating(Qm, TQm);


% Compute the reparametrizations to align all SRVFs to Qmn
[Gmr, TGmr, Gsr, TGsr] = srvf_fa_groupwise_reparam(Qmn, TQmn, Qsn, TQsn);


% Apply the reparametrizations
for i=1:nfuncs
  [Fsr{i}, TFsr{i}] = plf_compose(Fsn{i}, Tsn{i}, Gsr{i}, TGsr{i});
end
[Fmr, TFmr] = plf_compose(Fmn, Tmn, Gmr, TGmr);


% Plot the normalized, groupwise-aligned functions
figure();
clf();
hold on;
for i=1:nfuncs
  plot(TFsr{i}, Fsr{i});
end
