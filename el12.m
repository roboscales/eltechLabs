I = [0.48; 0.5; 0.6; 0.7; 0.8; 0.9];
Ui = [22.8; 23.3; 29.6; 34.1; 40.1; 42.9];
Un = [27.2; 26.7; 20.4; 16.0; 10.4; 6.5];
I2 = [0.3; 0.2; 0.4; 0.5; 0.7; 0.8];
I3 = [0.3; 0.25; 0.2; 0.15; 0.1; 0.05];

Ri = zeros(6, 1);

for i = 1:6
    Ri(i, 1) = ri(Ui(i,1), I(i,1));
end

display(Ri);

R2 = zeros(6, 1);

for i = 1:6
    R2(i, 1) = ri(Un(i,1), I2(i,1));
end

display(R2);

R3 = zeros(6, 1);

for i = 1:6
    R3(i, 1) = ri(Un(i,1), I3(i,1));
end

display(R3);

Rn = zeros(6, 1);

for i = 1:6
    Rn(i, 1) = rn2(R2(i,1), R3(i,1));
end

display(Rn);

P = zeros(6,1)

for l = 1:6
    P(l,1) = p(50, 50, Rn(l,1));
end

display (P);

Pn = zeros(6,1);

for m = 1:6
    Pn(m,1) = pn(Un(m,1), I(m,1));
end

display (Pn);

Cpd = zeros(6,1);

for n = 1:6
    Cpd(n,1) = cpd(Pn(n,1), P(n,1));
end

display (Cpd);