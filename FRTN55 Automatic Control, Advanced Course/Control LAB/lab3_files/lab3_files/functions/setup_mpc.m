function mpcobj = setup_mpc(lt, pumpConstraints, tankConstraints,q, r,...
    predictionHorizon, firstPenaltySample)

% Recover the discrete-time linear time-invariant system.
sys = lt.DtSys;

% Specify inputs, outputs and states
sys.InputName =  {'u_1', 'u_2'};
sys.OutputName = {'h_1', 'h_2', 'h_3', 'h_4'};
sys.StateName = {'u_1', 'u_2', 'h_3', 'h_4'};

% Specify groupings
sys.InputGroup.MV = [1 2];          % Manipulated variables, u1 and u2
sys.OutputGroup.MO = [1 2];         % Measured Outputs h_1, ..., h_2
sys.OutputGroup.UO = [3 4];         % Unknown outputs h_3, ..., h_4

% Get constraints in linearized variables
[MV, OV] = getConstraints(lt, pumpConstraints, tankConstraints);

% Specify time-dependent weights for output variables
Weights.OV = kron([0.0*ones(firstPenaltySample-1, 1); ...
    ones(predictionHorizon - firstPenaltySample + 1,1)],q);
% Weights on manipulated variables are just r
Weights.MV = r;
% Do not penalize rate
Weights.ManipulatedVariablesRate = [0. 0.];

controlHorizon = predictionHorizon;

mpcobj = mpc(sys,...
    lt.SampleTime,...
    predictionHorizon,...
    controlHorizon,...
    Weights,MV,OV);
end
