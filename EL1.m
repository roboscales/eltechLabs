%ok now

%begin data
I = [0.26; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8];
U = [12.7; 15.2; 18.3; 25.9; 30.1; 33.5; 39.2];
Un = [37.3; 34.8; 32.1; 24.2; 20.0; 16.8; 11.0];

display (I);
display (U);
display (Un);

Ri = zeros(7,1);

for i = 1:7
    Ri(i, 1) = ri(U(i,1), I(i,1));
end

display(Ri);

Rn = zeros(7,1);

for j = 1:7
    Rn(j, 1) = rn(Un(j,1), I(j,1));
end

display(Rn);

Re = zeros(7,1);

for k = 1:7
    Re(k, 1) = re(50, Rn(k,1));
end

display(Re);

P = zeros(7,1);

for l = 1:7
    P(l,1) = p(50, 50, Rn(l,1));
end

display (P);

Pn = zeros(7,1);

for m = 1:7
    Pn(m,1) = pn(Un(m,1), I(m,1));
end

display (Pn);

Cpd = zeros(7,1);

for n = 1:7
    Cpd(n,1) = cpd(Pn(n,1), P(n,1));
end

display (Cpd);

