% c: Circulation model with HR = 75 b/min, Emax = 2.0, Emin = 0.06
c = Circulation(75, 2.0, 0.06);

% simulate c for 10 seconds
c.simulate(10);