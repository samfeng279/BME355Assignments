% h: HillTypeMuscle being simulated
% h has isometric F = 100, resting CE L = 1, resting SE length = 0.1
% hV: normalized velocity when a = 1, lM = 1, lT = 1.01
h = HillTypeMuscle(100, 1, 0.1);
h.plotCurves();
hV = h.getVelocity(1, 1, 1.01);
h.simulateIsometric();

% s: uncontrolled StabilityModel being simulated
s = StabilityModel();
s.simulate(0, 5);

% c: controlled StabilityModel being simulated
c = StabilityModel();
c.simulate(1, 10);