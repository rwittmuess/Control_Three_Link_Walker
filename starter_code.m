% ME193B/293B Feedback Control of Legged Robots
% UC Berkeley

clear
%% Define symbolic variables for cofiguration variables and mechanical parameters

syms q1 q2 q3 x y real
syms dq1 dq2 dq3 dx dy real
syms u1 u2 real

% Position Variable vector
q = [x;y;q1;q2;q3];

% Velocity variable vector
dq = [dx;dy;dq1;dq2;dq3];

% State vector
s = [q;dq];

% Inputs
u = [u1;u2];

% parameters
          
lL = 1;
lT = 0.5;
mL = 5;
mT = 10;
JL = 0;
JT = 0;
mH = 15;
g = 9.81;

% # of Degrees of freedom
NDof = length(q);


%% Problem 1: Lagrangian Dynamics
%Find the CoM position of each link

% Torso
pComTorso = [x + lT*sin(q3);...
             y + lT*cos(q3)];

% Leg 1
pComLeg1 = [x - lL/2*cos(deg2rad(270) - (q1 + q3));...
            y - lL/2*sin(deg2rad(270) - (q1 + q3))];

% Leg 2
pComLeg2 = [x + lL/2*cos(q2 + q3 - deg2rad(90));...
            y - lL/2*sin(q2 + q3 - deg2rad(90))];


% Leg 1
pLeg1 = [x - lL*cos(deg2rad(270) - (q1 + q3));...
            y - lL*sin(deg2rad(270) - (q1 + q3))];

% Leg 2
pLeg2 = [x + lL*cos(q2 + q3 - deg2rad(90));...
            y - lL*sin(q2 + q3 - deg2rad(90))];

% Find the CoM velocity of each link

% Torso
dpComTorso = simplify(jacobian(pComTorso, q)*dq);

% Leg 1
dpComLeg1 = simplify(jacobian(pComLeg1, q)*dq);

% Leg 2
dpComLeg2 = simplify(jacobian(pComLeg2, q)*dq);


% Find absolute angular velocity associated with each link:
% Torso
dq3Absolute = dq3;
% Leg 1
dq1Absolute = dq3 + dq1;
% Leg 2
dq2Absolute = dq3 + dq2;


% Total Kinetic energy = Sum of kinetic energy of each link
% Torso
KETorso = 0.5*mT*dpComTorso(1)^2 + 0.5*mT*dpComTorso(2)^2 + 0.5*JT*dq3Absolute^2;

% Leg 1
KELeg1 = 0.5*mL*dpComLeg1(1)^2 + 0.5*mL*dpComLeg1(2)^2 + 0.5*JL*dq1Absolute^2;

% Leg 2
KELeg2 = 0.5*mL*dpComLeg2(1)^2 + 0.5*mL*dpComLeg2(2)^2 + 0.5*JL*dq2Absolute^2;

% Hip
KEHip = 0.5*mH*dx^2 + 0.5*mH*dy^2;

% Total KE
KE = simplify(KETorso + KELeg1 + KELeg2 + KEHip);


% Total potential energy = Sum of Potential energy of each link
% Torso
PETorso = mT*g*pComTorso(2);

%Leg 1
PELeg1 = mL*g*pComLeg1(2);

% Leg 2
PELeg2 = mL*g*pComLeg2(2);

% Hip
PEHip = mH*g*y;

% Total PE
PE = simplify(PETorso + PELeg1 + PELeg2 + PEHip);


% Lagrangian
L = KE - PE;


% Equations of Motion
% Actuated variables
qActuated = [q1;q2];

% D, C, G, and B matrices
[D, C, G, B] = LagrangianDynamics(KE, PE, q, dq, qActuated);


%%  Dynamics of Systems with Constraints
%Compute the Ground reaction Forces

% Compute the position of the stance foot (Leg 1) 
pst = [x - lL*cos(deg2rad(270) - (q1 + q3));...
       y - lL*sin(deg2rad(270) - (q1 + q3))];


% Compute the jacobian of the stance foot
JSt = jacobian(pst, q);


% Compute the time derivative of the Jacobian
dJSt = simplify(  jacobian(JSt*dq, q)  ) ;

% Compute the Stance Force
% [D   -JST';   [d2q ;     [-C*dq - G + B*u;
%  JSt   0]   *  FSt]   =   -dJSt * dq] ;
A = [D   -JSt';
     JSt  zeros(2)] ;
b = [-C*dq - G + B*u;
     -dJSt * dq] ;
d2q_FSt = A\b ;
FSt = simplify(  d2q_FSt(NDof+1:end)  ) ;


%% Impact Map

% Compute the swing leg position (leg 2)
pSw = [x + lL*cos(q2 + q3 - deg2rad(90));...
       y - lL*sin(q2 + q3 - deg2rad(90))];

JSw = jacobian(pSw, q);

% postImpact = [qPlus;F_impact];
% Here, q, dq represent the pre-impact positions and velocities
[postImpact] = ([D, -JSw';JSw, zeros(2)])\[D*dq;zeros(2,1)];

% Post Impact velocities
dqPlus = simplify(postImpact(1:NDof));

% Impact Force Magnitude
Fimpact = simplify(postImpact(NDof+1:NDof+2));


%% Other functions

% swing foot velocity
dpSw = JSw*dq;

%% Export functions
if ~exist('./gen')
    mkdir('./gen')
end
addpath('./gen')

matlabFunction(FSt,      'File', 'gen/Fst_gen', 'Vars', {s, u});
matlabFunction(dqPlus,   'File', 'gen/dqPlus_gen', 'Vars', {s});
matlabFunction(pSw,      'File', 'gen/pSw_gen', 'Vars', {s});
matlabFunction(dpSw,     'File', 'gen/dpSw_gen', 'Vars', {s});
matlabFunction(pst,      'File', 'gen/pSt_gen', 'Vars', {s});
matlabFunction(pComLeg1, 'File', 'gen/pComLeg1_gen', 'Vars', {s});
matlabFunction(pComLeg2, 'File', 'gen/pComLeg2_gen', 'Vars', {s});
matlabFunction(pComTorso,'File', 'gen/pComTorso_gen', 'Vars', {s});
matlabFunction(pLeg1,    'File', 'gen/pLeg1_gen', 'Vars', {s});
matlabFunction(pLeg2,    'File', 'gen/pLeg2_gen', 'Vars', {s});




%% [Part 1a] Compute the f and g vectors
% Fill up this


%% Change of Coordinates
% Transformation matrix:
T = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 1;
     0 0 0 1 1;
     0 0 0 0 1];
d = [0;
     0;
     -pi;
     -pi;
     0];

%% [Part 1b] Output dynamics
% Fill up this





%% [Part 1c] Lie Derivatives
% Fill up this






%% [Part 1d] Relabelling Matrix
% Fill up this


