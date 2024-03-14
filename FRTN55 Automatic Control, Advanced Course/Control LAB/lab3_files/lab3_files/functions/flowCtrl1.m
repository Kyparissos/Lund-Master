function u = flowController(yyref)

persistent isInit yMin yMax y yref yold P I v lastTime

% parameters
K = 0.21;
Ti = 0.06;
Tr = Ti;
h = 0.01;

if isempty(isInit)
    isInit = 1;  % initiation flag
    yMin = 0.0;  % nominal min venturi pressure
    yMax = 10.0; % nominal max venturi pressure
    y = 0.0;
    yref = 0.5;  % necessary in order not to miscalibrate
    yold = 0.0;
    P = 0.0;
    I = 0.0;
    v  =0.0;
    lastTime = 0;
end

% extract raw inputs
yRaw = yyref(1);
yref = yyref(2);

%rescale input
yref = yref/10;

% venturi pressure to flow
y = sqrt(max((yRaw-yMin),0)/(yMax-yMin));

%calculate output
v = 0.0;
P = -K*y;
v = v+P;
v = v + I;
u = v;

if yref <= 0.0
    u = 0.0;
elseif yref >= 1.0
    u = 1.0;
end

% update state
I = I + h*K/Ti*(yref-y) + h*(1/Tr)*(u-v);
yold = y;

% rescale output
u = 10*u;

