clear;
clc;

%% ============================================================
% 1. MODEL INPUTS AND EXPERIMENTAL CONDITIONS
% ============================================================

% Operating conditions
vel = 0.05:0.05:0.7;          % Velocity [m/s]
Supply_Pressure = [4 6 8];    % Supply pressure [bar]

% Cylinder dimensions
Dc = [8 10 12 16 20 25 32 40 50 63];  % Bore diameter [mm]
Dh = [4 4 6 6 8 10 12 16 20 20];      % Rod diameter [mm]

% Effective area ratios
AbAb = [0.75 0.84 0.75 0.8594 0.84 ...
        0.84 0.8594 0.84 0.84 0.8992];

% Model constants
b  = 0.3;
Po = 1.013e5;                 % Atmospheric pressure [Pa]

% Number of conditions
nCyl      = numel(Dc);
nVel      = numel(vel);
nPressure = numel(Supply_Pressure);

%% ============================================================
% 2. FRICTION-FORCE MODEL COEFFICIENTS
% ============================================================

% Bfr = [B1 B2 B3 B4 B5 B6 B7 B8]
Bfr = [0.419536900707924 ...
    0.0729716605961434 ...
    -0.0661216837044507 ...
    4.09689940676814 ...
    -0.191343440791252 ...
    -0.277014346294269 ...
    1.33234982210445 ...
    1.35799377866848]; 


%% ============================================================
% 3. POOLED MODEL DEFINITION — EQUATION (8)
% ============================================================

% Input columns:
% X(:,1) = velocity [m/s]
% X(:,2) = supply pressure [bar]
% X(:,3) = optimal friction force [N]
%
% Output:
% Ud = friction coefficient

yfcn = @(B,X) ...
    X(:,1).^B(7) .* ( ...
    B(1) .* (log(X(:,3)) + B(4).*X(:,2)) ...
    + B(2) .* (log(X(:,3)) + B(5).*X(:,2)).^2 ...
    + B(3) .* (log(X(:,3)) + B(6).*X(:,2)).^3 ...
    );

% Initial parameter estimates
b0 = [0.26 -0.10 0.015 0.76 0.22 0.014 0.54]';


% Optimization options
options = optimoptions( ...
    'lsqcurvefit', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', 1e5, ...
    'MaxIterations', 1e4);

%% ============================================================
% 4. INITIALIZE DATASET
% ============================================================

% Columns:
% 1 = cylinder ID
% 2 = velocity [m/s]
% 3 = supply pressure [bar]
% 4 = optimal friction force [N]
% 5 = reference friction coefficient Ud

nDataExpected = nCyl*nVel*nPressure;

Data_set = NaN(nDataExpected,5);

cont = 1;

%% ============================================================
% 5. GENERATE DATASET
% ============================================================

for p = 1:nPressure

    % Supply pressure [Pa]
    Ps = Supply_Pressure(p)*1e5;

    for k = 1:nVel

        % Operating velocity [m/s]
        v = vel(k);

        for j = 1:nCyl

            %% ------------------------------------------------
            % 5.1 Cylinder geometry
            %% ------------------------------------------------

            Aa = pi*(Dc(j)/1000)^2/4;

            Dsum = Dc(j) + Dh(j);

            Ra = AbAb(j);

            Ab = Aa*Ra;
            
            
            %% ------------------------------------------------
            % Friction model
            %
            % thetaFold_Ext(j,:) corresponds to the LOCO
            % coefficients for the friction-force model.
            %% ------------------------------------------------

            Fat = ...
                exp(-Bfr(1)/( ...
                vel(k)^Bfr(7) + ...
                (Dc(j)+Dh(j))^Bfr(5) + ...
                Bfr(6))) ...
                .* exp(Bfr(2)*(Dc(j)+Dh(j))) ...
                .* ( ...
                Bfr(3)*(Dc(j)+Dh(j)) ...
                + Bfr(4)*log(Ps/1e5 + 1.013 + Bfr(8)) );


            %% ------------------------------------------------
            % Friction coefficient
            %% ------------------------------------------------

            Ud_ref = Fat/(Ps*Aa);


            %% ------------------------------------------------
            % Geometry ratio
            %% ------------------------------------------------

            Ra = AbAb(j);


%% ------------------------------------------------
% Minimum pressure ratio (retraction)
%% ------------------------------------------------

popasmin = real( ...
    ( ...
    b*(1/Ra) + ...
    sqrt( ...
    b^2*(1/Ra)^2 ...
    + 5*b^2 ...
    + 1 ...
    - 4*b ...
    - 2*b^3 ...
    + (1/Ra)^2 ...
    - 2*b*(1/Ra)^2 ...
    ) ...
    )*(1/Ra) ...
    /(1 - 2*b + b^2 + (1/Ra)^2) ...
    * 1.00001 );


%% ------------------------------------------------
% Calculate Fc_Opt (retraction)
%% ------------------------------------------------

pbps = [];
popa = [];
Vc_Fc = [];

ii = 1;

Fr_coef = Ud_ref;

for pb = Ps*popasmin:5:Ps

    pbps(ii) = pb/Ps;

    popa(ii) = ...
        b + ...
        (1/Ra)*sqrt( ...
        (b^2 - 2*b + 1)/((1/Ra)^2) ...
        + 1 ...
        + pbps(ii)^-1*( ...
        -2*b + ...
        pbps(ii)^-1*(2*b - 1) ...
        ) ...
        );


    Vc_Fc(ii) = real( ...
        Ps*Po/pbps(ii) ...
        * sqrt( ...
        1 - ((pbps(ii) - b)/(1 - b))^2 ...
        ) ...
        * ( ...
        pbps(ii) ...
        + Po*(1 - Ra)/(Ps*Ra) ...
        - Fr_coef/Ra ...
        - (Po/(Ps*Ra))/popa(ii) ...
        ) ...
        );


    ii = ii + 1;

end


%% ------------------------------------------------
% Maximum flow condition
%% ------------------------------------------------

[~,I] = max(Vc_Fc);

pBpSbest = pbps(I);
p0pAbest = popa(I);


%% ------------------------------------------------
% Optimal friction force (retraction)
%% ------------------------------------------------

Fc_Opt = abs( ...
    Ab * ( ...
    Ps*( ...
    pBpSbest ...
    - (1/p0pAbest)*(Po/Ps)*(1/Ra) ...
    + (Po/Ps)*(1 - Ra)/Ra ...
    - Fr_coef/Ra ...
    ) ...
    ) ...
    );

            %% ------------------------------------------------
            % 5.8 Store observation
            %% ------------------------------------------------

            Data_set(cont,:) = [ ...
                j, ...
                v, ...
                Ps/1e5, ...
                Fc_Opt, ...
                Ud_ref];

            cont = cont + 1;

        end
    end
end

%% ============================================================
% 6. CLEAN DATASET AND REPORT SIZE
% ============================================================

% Remove any rows that were not populated.
Data_set = Data_set(1:cont-1,:);

nData = size(Data_set,1);

fprintf('\n');
fprintf('============================================\n');
fprintf('GENERATED DATASET\n');
fprintf('============================================\n');

fprintf('Total observations: %d\n',nData);
fprintf('Number of cylinders: %d\n',nCyl);

if nCyl > 0
    fprintf('Expected observations per cylinder: %d\n', ...
        nVel*nPressure);
end

%% ============================================================
% 7. PREPARE DATA FOR POOLED PARAMETER IDENTIFICATION
% ============================================================

X_pooled = Data_set(:,2:4);
Y_pooled = Data_set(:,5);

% Retain only observations with valid inputs and outputs.
validData = ...
    all(isfinite(X_pooled),2) ...
    & isfinite(Y_pooled) ...
    & X_pooled(:,3) > 0;

X_pooled = X_pooled(validData,:);
Y_pooled = Y_pooled(validData);

fprintf('Valid observations for fitting: %d\n', ...
    numel(Y_pooled));

%% ============================================================
% 8. POOLED PARAMETER IDENTIFICATION — EQUATION (8)
% ============================================================

B_pooled = lsqcurvefit( ...
    yfcn, ...
    b0, ...
    X_pooled, ...
    Y_pooled, ...
    [], ...
    [], ...
    options);

%% ============================================================
% 9. DISPLAY POOLED PARAMETERS
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('POOLED PARAMETERS — EQUATION (8)\n');
fprintf('============================================\n');

for k = 1:numel(B_pooled)
    fprintf('B%d = %.8f\n',k,B_pooled(k));
end

%% ============================================================
% 10. POOLED PREDICTIONS
% ============================================================

% In-sample predictions using the same observations used for fitting.
Ud_pooled_pred = yfcn(B_pooled,X_pooled);

Error_pooled = Ud_pooled_pred - Y_pooled;

%% ============================================================
% 11. POOLED MODEL PERFORMANCE
% ============================================================

% Mean absolute error
MAE_pooled = mean(abs(Error_pooled));

% Root mean squared error
RMSE_pooled = sqrt(mean(Error_pooled.^2));

% Mean absolute percentage error
nonzeroY = Y_pooled ~= 0;

if any(nonzeroY)
    MAPE_pooled = ...
        100*mean( ...
            abs(Error_pooled(nonzeroY)) ...
            ./abs(Y_pooled(nonzeroY)) ...
        );
else
    MAPE_pooled = NaN;
end

% Mean prediction error (bias)
Bias_pooled = mean(Error_pooled);

% Coefficient of determination
SSres_pooled = sum(Error_pooled.^2);

SStot_pooled = sum( ...
    (Y_pooled-mean(Y_pooled)).^2);

if SStot_pooled > 0
    R2_pooled = 1 - SSres_pooled/SStot_pooled;
else
    R2_pooled = NaN;
end

%% ============================================================
% 12. DISPLAY PERFORMANCE METRICS
% ============================================================

fprintf('\n');
fprintf('============================================\n');
fprintf('POOLED MODEL PERFORMANCE — IN-SAMPLE\n');
fprintf('============================================\n');

fprintf('Observations = %d\n',numel(Y_pooled));
fprintf('MAE          = %.6f\n',MAE_pooled);
fprintf('RMSE         = %.6f\n',RMSE_pooled);
fprintf('MAPE         = %.2f %%\n',MAPE_pooled);
fprintf('Bias         = %.6f\n',Bias_pooled);
fprintf('R2           = %.4f\n',R2_pooled);


