% h: HillTypeMuscle being simulated
% h has isometric F = 100, resting CE L = 1, resting SE length = 0.1
h = HillTypeMuscle(100, 1, 0.1);
h.plotCurves();
h.getVelocity(1, 1, 1.01);
h.simulateIsometric();

