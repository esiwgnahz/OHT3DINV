%% Cleanup

clear all; close all; clc;

%% Describe the setup for the forward models

% Bounds of the domain, and number of points in the discretization (last)
x_disc = [-22.5 22.5 20];
y_disc = [-22.5 22.5 20];
z_disc = [0 15 15];
xmin = x_disc(1); xmax = x_disc(2);
ymin = y_disc(1); ymax = y_disc(2);
zmin = z_disc(1); zmax = z_disc(2);


%Specify boundary types and boundary values (x / y constant head
%boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals);

% Locations of all wells (pumping and observation)
well_locs = [...
    0 0 3; ...
    0 0 6; ...
    0 0 9; ...
    0 0 12; ...
    0 10 3; ...
    0 10 6; ...
    0 10 9; ...
    0 10 12; ...
    10 0 3; ...
    10 0 6; ...
    10 0 9; ...
    10 0 12; ...
    0 -10 3; ...
    0 -10 6; ...
    0 -10 9; ...
    0 -10 12; ...
    -10 0 3; ...
    -10 0 6; ...
    -10 0 9; ...
    -10 0 12; ...
    ];


V = 0.005;
P = [10 50 100];
test_list = [];
for i = 1:1:numel(P)
    for s = 17:1:20
        for r = 1:1:16
            test_list = [test_list; ...
                (2*pi)/P(i) s V*pi/P(i) r; ...
                ];
        end
    end
end

% Periods (in time units) of the various pumping tests
periods = (2*pi)./test_list(:,1);

%% Create the "official" input files needed by the steady periodic model.

%This step should automatically create the input files needed to run all of
%the steady periodic models. 

[domain] = equigrid_setup(x_disc,y_disc,z_disc);

[experiment] = OHT_create_inputs(well_locs,test_list,domain);


%% Calculate counts (used often for plotting and other functionality)

% Calculates number of observations, size of grid, and sets up a
% meshgridded set of points to do plotting of results.

num_omegas = size(experiment,1);
numobs = size(test_list,1);
num_x = numel(domain.x) - 1;
num_y = numel(domain.y) - 1;
num_z = numel(domain.z) - 1;
num_cells = num_x*num_y*num_z;

%% Setup grid of cell centers for plotting and geostatistical setups
[coords, numgrid] = plaid_cellcenter_coord(domain);

%% For synthetic problem - true parameter field statistics

mean_lnK = -9.2;
var_lnK = 2;

mean_lnSs = -11.5;
var_lnSs = 0.01;

corr_x = 5;
corr_y = 10;
corr_z = 2;

distmat_row = dimdist(coords(1,:),coords);
corr_row = exp(-(...
    (distmat_row(:,:,1)./corr_x).^2 + ...
    (distmat_row(:,:,2)./corr_y).^2 + ...
    (distmat_row(:,:,3)./corr_z).^2).^.5);

randn('state',0)
[corr_relz] = toepmat_vector_math(corr_row,'r',[],3,[num_y num_x num_z]);

lnK_true = corr_relz(:,1).*var_lnK.^.5 + mean_lnK;
lnSs_true = corr_relz(:,2).*var_lnSs.^.5 + mean_lnSs;
params_true = [lnK_true; lnSs_true];

lnK_true_grid = reshape(lnK_true,num_y,num_x,num_z);
lnSs_true_grid = reshape(lnSs_true,num_y,num_x,num_z);

K_range = [min(lnK_true) max(lnK_true)];
Ss_range = [min(lnSs_true) max(lnSs_true)];

figure(1)
subplot(1,2,1)
sl1 = slice(numgrid{1},numgrid{2},numgrid{3},lnK_true_grid,0,0,5);
set(sl1,'LineStyle','none','FaceColor','interp')
colorbar
axis equal
axis([xmin xmax ymin ymax zmin zmax]);
title('True ln(K)');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
colormap jet
set(gca,'FontSize',16);

caxis(K_range)
subplot(1,2,2)
sl2 = slice(numgrid{1},numgrid{2},numgrid{3},lnSs_true_grid,0,0,5);
set(sl2,'LineStyle','none','FaceColor','interp')
colorbar
axis equal
axis([xmin xmax ymin ymax zmin zmax]);
title('True ln(Ss)');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
colormap jet
set(gca,'FontSize',16);

%% Perform all model runs to generate data

% Steady_periodic_multimodel_run performs the model runs to generate all
% observations. 
% IMPORTANT:
% The first output is simulated observations, output as
% follows: [A_obs(1); A_obs(2); A_obs(3) ... A_obs(numobs); B_obs(1) ...
% B_obs(numobs)]
%
% NOTE: B in this case is the imaginary part of the phasor. Note that this
% is not the same as the coefficient in front of the sine (which would be
% negative B)
%
% The second output is the full phasor field for each pumping test.
% TODO: May change this so that it is the given pumping test for each
% observation. Also, need to add description of structure.
%
% H is the sensitivity matrix. Each row represents a different observation,
% organized the same way as sim_obs. Each column represents a different
% parameter sensitivity - the first set of columns is with respect to each
% K value, and the second set of columns is with respect to each Ss value.
% K and Ss grid cells are stepped through in sequential order by y, then x,
% then z, as would be done by meshgrid.

h_function = @(params) OHT_run_distribKSs(params, ...
    domain,bdrys,experiment,1);

h_homog_function = @(KSs) h_function([KSs(1)*ones(num_cells,1); KSs(2)*ones(num_cells,1)]);

Phi_function = @(params) OHT_run_distribKSs(params, ...
    domain,bdrys,experiment,2);

H_tilde_func = @(params) OHT_run_distribKSs(params, ...
    domain,bdrys,experiment,3);

tic
sim_obs = h_function(params_true);
toc

tic
Phi_true = Phi_function(params_true);
toc

tic

params_init = [-9*ones(num_cells,1); -9*ones(num_cells,1)];

tic
H_adj = H_tilde_func(params_init);
toc

tic
Phi_init = Phi_function(params_init);
toc

%% Cross-checks - visualize results and sensitivities

%These cross-checks are used simply to make sure that the model appears to
%be running correctly. We calculate the amplitude and phase first (which
%are more intuitive than A and B), and then look at how these vary with
%space.
%
%The second set of plots looks at the sensitivity structure of the signal
%amplitude to K and Ss, which should look like "blobs" centered around the
%pumping and observation well for each observation.

A_obs = sim_obs(1:numobs);
B_obs = sim_obs((numobs+1):(2*numobs));
amp = (A_obs.^2 + B_obs.^2).^.5;
phase = atan2(-B_obs,A_obs);

A_full = Phi_true(1:num_cells,:);
B_full = Phi_true((num_cells+1):(2*num_cells),:);
amp_full = (A_full.^2 + B_full.^2).^.5;
phase_full = atan2(-B_full,A_full);

H_AK = H_adj((1:numobs),(1:num_cells));
H_BK = H_adj((numobs+1):(2*numobs),(1:num_cells));

H_ASs = H_adj((1:numobs),((num_cells+1):(2*num_cells)));
H_BSs = H_adj((numobs+1):(2*numobs),((num_cells+1):(2*num_cells)));

%Check all phasor fields produced through numerical model running
num_totalstims = size(Phi_true,2);
figure(2)
for i = 1:1:num_totalstims
    amp_field = reshape(log(amp_full(:,i)),num_y,num_x,num_z);
    phase_field = reshape(phase_full(:,i),num_y,num_x,num_z);
    subplot(1,2,1);
    sl6a = slice(numgrid{1},numgrid{2},numgrid{3},amp_field,[0],[0],[5]);
    set(sl6a,'FaceColor','interp','LineStyle','none')
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    title('ln(Amplitude)')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    colormap jet
    subplot(1,2,2);
    sl6b = slice(numgrid{1},numgrid{2},numgrid{3},phase_field,[0],[0],[5]);
    set(sl6b,'FaceColor','interp','LineStyle','none')
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    title('Phase')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    colormap jet
    subplot(1,2,2)
    pause
end

figure(3);
for i = 1:1:numobs
    ampK_sens = (A_obs(i)./amp(i)).*H_AK(i,:) + (B_obs(i)./amp(i)).*H_BK(i,:);
    ampK_sensmap = reshape(ampK_sens,num_y,num_x,num_z);
    subplot(1,2,1)
    sl3 = slice(numgrid{1},numgrid{2},numgrid{3},ampK_sensmap,[0],[0],[5]);
    set(sl3,'FaceColor','interp','LineStyle','none')
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    title('Amplitude sensitivity to ln(K)')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    colormap jet

    ampSs_sens = (A_obs(i)./amp(i)).*H_ASs(i,:) + (B_obs(i)./amp(i)).*H_BSs(i,:);
    ampSs_sensmap = reshape(ampSs_sens,num_y,num_x,num_z);
    subplot(1,2,2)
    sl4 = slice(numgrid{1},numgrid{2},numgrid{3},ampSs_sensmap,[0],[0],[5]);
    set(sl4,'FaceColor','interp','LineStyle','none')
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    title('Amplitude sensitivity to ln(Ss)')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    colormap jet

    pause
end

%% Prior setup

var_lnK_est = 2;
var_lnSs_est = 0.01;

%For now assume that the Q matrices are known exactly
QK_row = var_lnK_est.*corr_row;
QSs_row = var_lnSs_est.*corr_row;

Qprod_func = @(invec) covar_product_K_Ss(QK_row,QSs_row,invec,num_x,num_y,num_z);

%% Inversion

lnK_homog_guess = -9;
lnSs_homog_guess = -11;
data_errorcov_guess = 1e-9;

params_init = [lnK_homog_guess*ones(num_cells,1); lnSs_homog_guess*ones(num_cells,1)];
beta_init = [lnK_homog_guess; lnSs_homog_guess];
X = [ones(num_cells,1) zeros(num_cells,1); zeros(num_cells,1) ones(num_cells,1)];
R = data_errorcov_guess.*eye(2*numobs);
y = sim_obs;

tic;
[params_best, beta_best, H_local, negloglike] = ql_geostat_inv(y,params_init,beta_init,X,R,Qprod_func,h_function,H_tilde_func);
inv_time = toc;

negloglike_func = @(s,beta) negloglike_eval(y,X,s,beta,Q,R,h_function);

%% Plotting of true vs. results

lnK_best_grid = reshape(params_best(1:num_cells),num_y,num_x,num_z);
lnSs_best_grid = reshape(params_best((num_cells+1):(2*num_cells)),num_y,num_x,num_z);

marksize = 8;
faceal = 0.7;

figure(4)
subplot(1,2,1)
pl = plot3(well_locs(:,1),well_locs(:,2),well_locs(:,3),'ok');
set(pl,'MarkerSize',marksize,'MarkerFaceColor','k')
hold on
sl1 = slice(numgrid{1},numgrid{2},numgrid{3},lnK_true_grid,0,0,5);
hold off
set(sl1,'LineStyle','none','FaceColor','interp','FaceAlpha',faceal)
colorbar
colormap('jet')
caxis(K_range)
axis equal
axis([xmin xmax ymin ymax zmin zmax])
title('True ln(K[m/s])')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
set(gca,'FontSize',16)
grid on

subplot(1,2,2)
pl = plot3(well_locs(:,1),well_locs(:,2),well_locs(:,3),'ok');
set(pl,'MarkerSize',marksize,'MarkerFaceColor','k')
hold on
sl2 = slice(numgrid{1},numgrid{2},numgrid{3},lnSs_true_grid,0,0,5);
hold off
set(sl2,'LineStyle','none','FaceColor','interp','FaceAlpha',faceal)
colorbar
colormap('jet')
caxis(Ss_range);
axis equal
axis([xmin xmax ymin ymax zmin zmax])
title('True ln(Ss[1/m])')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
set(gca,'FontSize',16)
grid on

figure(5)
subplot(1,2,1)
pl = plot3(well_locs(:,1),well_locs(:,2),well_locs(:,3),'ok');
set(pl,'MarkerSize',marksize,'MarkerFaceColor','k')
hold on
sl1 = slice(numgrid{1},numgrid{2},numgrid{3},lnK_best_grid,0,0,5);
hold off
set(sl1,'LineStyle','none','FaceColor','interp','FaceAlpha',faceal)
colorbar
colormap('jet')
axis equal
axis([xmin xmax ymin ymax zmin zmax])
title('Estimated ln(K[m/s])')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
set(gca,'FontSize',16)
grid on

subplot(1,2,2)
pl = plot3(well_locs(:,1),well_locs(:,2),well_locs(:,3),'ok');
set(pl,'MarkerSize',marksize,'MarkerFaceColor','k')
hold on
sl2 = slice(numgrid{1},numgrid{2},numgrid{3},lnSs_best_grid,0,0,5);
hold off
set(sl2,'LineStyle','none','FaceColor','interp','FaceAlpha',faceal)
colorbar
colormap('jet')
axis equal
axis([xmin xmax ymin ymax zmin zmax])
title('Estimated ln(Ss[1/m])')
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
set(gca,'FontSize',16)
grid on

v = [27 42];
figure(4); subplot(1,2,1); view(v); caxis(K_range); subplot(1,2,2); view(v); caxis(Ss_range)
figure(5); subplot(1,2,1); view(v); caxis(K_range); subplot(1,2,2); view(v); caxis(Ss_range)

set(4,'Position',[100 100 2000 600])
set(5,'Position',[100 100 2000 600])


%% Save output

save('inverse_test_3D_out.mat')
