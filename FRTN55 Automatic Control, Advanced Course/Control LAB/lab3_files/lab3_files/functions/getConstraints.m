%GETCONSTRAINTS     Normalized constraints
%   [MV, OV] = GETCONSTRAINTS(lt, pumpConstraints, tankConstraints) returns
%   constraints MV, OV where
%   MV is a 1x2 struct with constraints on the Manipulated Variables,
%   OV is a 1x4 struct with constraints on the Output Variables and

function [MV, OV] = getConstraints(lt, pumpConstraints, tankConstraints)
arguments
    lt (1,1) LinearTankModel
    pumpConstraints (1,2) struct
    tankConstraints (1,4) struct
end
MV = MVConstraints(lt, pumpConstraints);
OV = OVConstraints(lt, tankConstraints);
end

%MVConstrainst  Normalized constraints on the manipulated variables
%   MV = MVConstraints(lt, pumpConstraints)
%   lt is a fully specified LinearTankModel object and
%   pumpConstraints is a 1x2 struct with the fields
%   'Min', 'Max', 'RateMin', 'RateMax'
%
%   MV is a 1x2 struct with the fields
%   'Min', 'Max', 'RateMin', 'RateMax'
function MV = MVConstraints(lt, pumpConstraints)
[u1Max, u2Max] = pumpConstraints.Max;
[u1Min, u2Min] = pumpConstraints.Min;
[u1RateMax, u2RateMax] = pumpConstraints.RateMax;
[u1RateMin, u2RateMin] = pumpConstraints.RateMin;

MV = struct('Max', {(u1Max - lt.LinPoint(5)) (u2Max - lt.LinPoint(6))},...
    'Min', {(u1Min - lt.LinPoint(5)) (u2Min - lt.LinPoint(6))}, ...
    'RateMax', {u1RateMax u2RateMax},...
    'RateMin', {u1RateMin u2RateMin});
end

%OVConstrainst  Normalized constraints on the manipulated variables
%   OV = MVConstraints(lt, pumpConstraints)
%   lt is a fully specified LinearTankModel object and
%   tankConstraints is a 1x4 struct with the fields
%   'Min', 'Max'
%
%   OV is a 1x4 struct with the fields
%   'Min', 'Max', 'RateMin', 'RateMax'
function OV = OVConstraints(lt, tankConstraints)
[h1Max, h2Max, h3Max, h4Max] = tankConstraints.Max;
[h1Min, h2Min, h3Min, h4Min] = tankConstraints.Min;

OV = struct('Max', num2cell(([h1Max h2Max h3Max h4Max] - lt.LinPoint(1:4))*...
    lt.MeasurementConstant), ...
    'Min', num2cell(([h1Min h2Min h3Min h4Min] - lt.LinPoint(1:4))*...
    lt.MeasurementConstant));
end