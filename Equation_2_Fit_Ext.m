%% EQUATION (2) - FRICTION MODEL IDENTIFICATION
% MATLAB 2025
% Nonlinear least squares using lsqnonlin
% Dataset: 406 points per stroke

clearvars -except F_Fr_Ext
clc;

%% Load experimental data

% Get the directory containing this script
repoDir = fileparts(mfilename('fullpath'));

% Construct the path relative to the repository root
dataFile = fullfile(repoDir, 'data', 'F_Fr_Ext.mat');

% Load data
load(dataFile, 'F_Fr_Ext');


%% 1. EXPERIMENTAL DATA

% Cylinder bore diameters [mm]
Dc_all = [8 10 12 16 20 25 32 40 50 63];

% Rod diameters [mm]
Dh_all = [4 4 6 6 8 10 12 16 20 20];

% Supply pressures [bar]
Ps_all = [4 6 8];

% Velocity points [m/s]
vel_all = linspace(50,700,14)/1000;

% Construct data matrices
[Ps,Dc,Vel] = ndgrid(Ps_all,Dc_all,vel_all);

% Calculate the sum of bore and rod diameters [mm]
Dsum = Dc + repmat(reshape(Dh_all,1,[],1),3,1,14);

% Flatten the experimental friction data
F = F_Fr_Ext(:);

% Construct input matrix:
% X(:,1) = Dc + Dh [mm]
% X(:,2) = Ps [bar]
% X(:,3) = velocity [m/s]
X = [Dsum(:), Ps(:), Vel(:)];

% Remove invalid and zero-valued observations
valid = isfinite(F) & all(isfinite(X),2) & (F ~= 0);

F = F(valid);
X = X(valid,:);

fprintf('Valid data points: %d\n',numel(F));
fprintf('Removed data points: %d\n',numel(F_Fr_Ext)-numel(F));

%% 2. INITIAL PARAMETER ESTIMATE

theta0 = [0.54 0.0719 -0.074 4.94 -0.029 ...
         -0.73 1.219 1.46];



%% 3. NONLINEAR LEAST-SQUARES FIT

% MATLAB lsqnonlin default options
options = optimoptions('lsqnonlin');

% Estimate parameters
[theta_fit,resnorm,residual,exitflag,output] = ...
    lsqnonlin(@(theta) modelFun(theta,X)-F, ...
              theta0,[],[],options);

%% 4. FITTED PARAMETERS

fprintf('\nFitted coefficients:\n');

for k = 1:numel(theta_fit)
    fprintf('theta(%d) = %.10g\n',k,theta_fit(k));
end

fprintf('\nResidual sum of squares = %.8g\n',resnorm);
fprintf('Exit flag = %d\n',exitflag);

disp(output);

%% 5. MODEL PERFORMANCE

F_pred = modelFun(theta_fit,X);

error = F_pred-F;

MAE  = mean(abs(error));
RMSE = sqrt(mean(error.^2));
Bias = mean(error);

SSres = sum(error.^2);
SStot = sum((F-mean(F)).^2);

R2 = 1-SSres/SStot;

% MAPE is undefined for zero measured values
idx = abs(F) > eps;

MAPE = mean(abs(error(idx)./F(idx)))*100;

fprintf('\nModel performance:\n');
fprintf('MAE  = %.4f N\n',MAE);
fprintf('RMSE = %.4f N\n',RMSE);
fprintf('Bias = %.4f N\n',Bias);
fprintf('MAPE = %.2f %%\n',MAPE);
fprintf('R^2  = %.6f\n',R2);

function yhat = modelFun(theta,X)

    Dc  = X(:,1);
    ps  = X(:,2);
    vel = X(:,3);

    yhat = exp(-theta(1)./( ...
        vel.^theta(7) + Dc.^theta(5) + theta(6))) ...
        .* exp(Dc.*theta(2)) ...
        .* (Dc.*theta(3) + ...
        theta(4).*log(ps + 1.013 + theta(8)));

end