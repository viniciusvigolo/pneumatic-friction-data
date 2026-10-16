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
Bfr = [ ...
    0.542559631070737, ...
    0.0719617110574386, ...
   -0.0741569315759701, ...
    4.94342707811763, ...
   -0.0295866295873569, ...
   -0.729961349398097, ...
    1.21981500437611, ...
    1.44753814092612];

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
    X(:,1) .* ( ...
    B(1) .* (log(X(:,3)) + B(4).*X(:,2)) ...
    + B(2) .* (log(X(:,3)) + B(5).*X(:,2)).^2 ...
    + B(3) .* (log(X(:,3)) + B(6).*X(:,2)).^3 ...
    );

% Initial parameter estimates
b0 = [0.68; -0.25; 0.024; 0; 0; 0];

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

            %% ------------------------------------------------
            % 5.2 Friction force
            %% ------------------------------------------------

            Fat = ...
                exp(-Bfr(1)/( ...
                    v^Bfr(7) + ...
                    Dsum^Bfr(5) + ...
                    Bfr(6))) ...
                * exp(Bfr(2)*Dsum) ...
                * ( ...
                    Bfr(3)*Dsum ...
                    + Bfr(4)*log(Ps/1e5 + 1.013 + Bfr(8)) ...
                  );

            %% ------------------------------------------------
            % 5.3 Reference friction coefficient
            %% ------------------------------------------------

            Ud_ref = Fat/(Ps*Aa);

            %% ------------------------------------------------
            % 5.4 Minimum pressure ratio
            %% ------------------------------------------------

            papsmin = real( ...
                ( ...
                    b*Ra + ...
                    sqrt( ...
                        b^2*Ra^2 ...
                        + 5*b^2 ...
                        + 1 ...
                        - 4*b ...
                        - 2*b^3 ...
                        + Ra^2 ...
                        - 2*b*Ra^2 ...
                    ) ...
                )*Ra ...
                /(1 - 2*b + b^2 + Ra^2) ...
                * 1.00001 ...
            );

            %% ------------------------------------------------
            % 5.5 Evaluate operating points
            %% ------------------------------------------------

            % Pressure-ratio increment retained from the original
            % calculation.
            pa_vec = Ps*papsmin:5:Ps;

            nPoints = numel(pa_vec);

            paps = NaN(nPoints,1);
            p0pb = NaN(nPoints,1);
            Vc_Fc = NaN(nPoints,1);

            for ii = 1:nPoints

                paps(ii) = pa_vec(ii)/Ps;

                % Pressure ratio p0/pB
                discriminant = ...
                    b^2*paps(ii)^2 ...
                    - 2*b*paps(ii)^2 ...
                    + paps(ii)^2 ...
                    + paps(ii)^2*Ra^2 ...
                    - Ra^2 ...
                    + 2*b*Ra^2 ...
                    - 2*b*paps(ii)*Ra^2;

                p0pb(ii) = ...
                    ( ...
                        b*paps(ii) + ...
                        sqrt(discriminant) ...
                    )/paps(ii);

                % Flow-related operating-point expression
                Vc_Fc(ii) = real( ...
                    ( ...
                        Ps*Po*sqrt( ...
                            1 - ((paps(ii)-b)/(1-b))^2 ...
                        ) ...
                        .* ( ...
                            paps(ii) ...
                            - Po*(1-Ra)/Ps ...
                            - Ud_ref ...
                            - (Po*Ra/Ps) * ...
                            ( ...
                                paps(ii) ./ ...
                                ( ...
                                    b*paps(ii) + ...
                                    sqrt(discriminant) ...
                                ) ...
                            ) ...
                        ) ...
                    )/paps(ii) ...
                );

            end

            %% ------------------------------------------------
            % 5.6 Identify maximum-flow operating point
            %% ------------------------------------------------

            validPoints = isfinite(Vc_Fc);

            if ~any(validPoints)
                warning( ...
                    'No valid operating point for cylinder %d, v = %.3f m/s, Ps = %.1f bar.', ...
                    j, v, Ps/1e5);

                continue;
            end

            Vc_valid = Vc_Fc;
            Vc_valid(~validPoints) = -Inf;

            [~,I] = max(Vc_valid);

            pApSbest = paps(I);
            p0pBbest = p0pb(I);

            %% ------------------------------------------------
            % 5.7 Calculate optimal friction force
            %% ------------------------------------------------

            Fc_Opt = abs( ...
                Aa * ( ...
                    Ps*( ...
                        pApSbest ...
                        - (1/p0pBbest)*(Po/Ps)*Ra ...
                        - (Po/Ps)*(1-Ra) ...
                        - Ud_ref ...
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

