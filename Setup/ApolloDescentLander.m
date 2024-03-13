%% Lunar Icecube Spacecraft System Parameters
% This function initializes Earth-Moon system parameters
function S = ApolloDescentLander(S)

% Mass Parameters:
S.WetMass = 15334; % kg
S.DPSPropellantMass = 8200; % kg
S.RCSPropellantMass = 287; % kg
S.DryMass = S.WetMass-S.DPSPropellantMass-S.RCSPropellantMass; % kg
S.InertiaTensor = 1.356*[25810 63 523;
                         63 27392 183;
                         523 183 26025];

% Angles:
S.GlidescopeAngle = 70; % degrees
S.MaximumPointingAngle = 90; % deg

% Thrust:
S.DPSMaximumThrust = 45040*.6; % N
S.DPSMinimumThrust = 45040*.1; % N
S.RCSMaximumThrust = 440; % N

% RCS Configuation
S.RCSThruster(1).Location = [1.6;1.6;1.65];
S.RCSThruster(1).Direction = [1;0;0];
S.RCSThruster(2).Location = [1.6;1.6;1.65];
S.RCSThruster(2).Direction = [0;1;0];
S.RCSThruster(3).Location = [1.6;1.6;1.65];
S.RCSThruster(3).Direction = [0;0;1];
S.RCSThruster(4).Location = [-1.6;1.6;1.65];
S.RCSThruster(4).Direction = [-1;0;0];
S.RCSThruster(5).Location = [-1.6;1.6;1.65];
S.RCSThruster(5).Direction = [0;1;0];
S.RCSThruster(6).Location = [-1.6;1.6;1.65];
S.RCSThruster(6).Direction = [0;0;1];
S.RCSThruster(7).Location = [1.6;-1.6;1.65];
S.RCSThruster(7).Direction = [1;0;0];
S.RCSThruster(8).Location = [1.6;-1.6;1.65];
S.RCSThruster(8).Direction = [0;-1;0];
S.RCSThruster(9).Location = [1.6;-1.6;1.65];
S.RCSThruster(9).Direction = [0;0;1];
S.RCSThruster(10).Location = [-1.6;-1.6;1.65];
S.RCSThruster(10).Direction = [-1;0;0];
S.RCSThruster(11).Location = [-1.6;-1.6;1.65];
S.RCSThruster(11).Direction = [0;-1;0];
S.RCSThruster(12).Location = [-1.6;-1.6;1.65];
S.RCSThruster(12).Direction = [0;0;1];

% Rates:
S.MaximumRotationalRate = 90*pi/180; % rad/s
S.MaximumVelocity = 200; % m/s

% Environment:
S.g = 1.625;%9.81;

% Propulsion System:
S.RCSSpecificImpulse = 290; % s
S.DPSSpecificImpulse = 311; % s
S.alphaRCS = 1/(S.RCSSpecificImpulse*9.81);
S.alphaDPS = 1/(S.DPSSpecificImpulse*9.81);

end