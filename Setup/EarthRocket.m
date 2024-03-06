%% Lunar Icecube Spacecraft System Parameters
% This function initializes Earth-Moon system parameters
function S = EarthRocket(S)

% Mass Parameters:
S.WetMass = 30000; % kg
S.DryMass = 22000; % kg
S.InertiaTensor = diag([4000000,4000000,100000]);
S.GimbaltoCOMDistance = 14; % m

% Angles:
S.MaximumGimbalAngle = 7; % degrees
S.MaximumPointingAngle = 70; % degrees
S.GlidescopeAngle = 70; % degrees

% Thrust:
S.MaximumThrust = 800000; % N
S.MinimumThrust = 320000; % N

% Rates:
S.MaximumRotationalRate = 90*pi/180; % rad/s
S.MaximumVelocity = 100; % m/s

% Environment:
S.g = 9.81;

% Propulsion System:
S.Isp = 282;
S.alpha = 1/(S.Isp*9.81);

end