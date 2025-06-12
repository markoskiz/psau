%% 3.19 (diskreten Diferencijalen Kompenzator)
clc
clear

syms s z Zc Pc Zd Pd T Kd

Gd(s) = (s+Zc)/(s+Pc);

Gdz(z) = Kd*(z-Zd)/(z-Pd);

% Zc -> nula na kontinualniot kompenzator
% Pc -> pol na kontinualniot kompenzator

% Zd -> nula na diskretniot kompenzator
% Pd -> pol na diskretniot kompenzator

% Gd(0) == Gdz(1)

solution = solve(Gd(0)==Gdz(1), Kd);
pretty(solution) % Kd -> statichko zasiluvanje
Kd = solution;

% z = e^(s*T)

Zd = exp(-Zc*T)
Pd = exp(-Pc*T)

Gd(z) = Kd * (z-Zd)/(z-Pd);
pretty(Gd)

%% 3.21 (diskreten Integralen Kompenzator)
clc
clear

syms s Zc Pc Zd Pd T z Kd

Gd(s) = Pc/Zc * (s+Zc)/(s+Pc)

Giz(z) = (1-Pd)/(1-Zd) * (z-Zd)/(z-Pd)

% Zc -> nula na kontinualniot kompenzator
% Pc -> pol na kontinualniot kompenzator

% Zd -> nula na diskretniot kompenzator
% Pd -> pol na diskretniot kompenzator

% z = exp(s*T)

Zd = exp(-Zc*T)
Pd = exp(-Pc*T)

Gi(z) = (1-Pd)/(1-Zd) * (z-Zd)/(z-Pd);
pretty(Gi)

%% 3.5
clc
clear

s = tf('s');

P = exp(-0.2*s)/(s+5);

K = 1;
R = K;

G0 = R*P;

%% a) K = ? t.sh. sistemot kje se naogja na gr. na stabilnost

nyquist(G0) % -> mozhe da se vnese zasiluvanje od 21.1 dB
hold on;

%% 20*log10(K) = 21.1 dB

K = 10^(21.1/20)
R = K;

G0 = R*P;

nyquist(G0);
legend('G0', 'G0_comp');

%% b) K = ? t.sh. Phi_rf = 50 deg.
clc

K = 1;
R = K;
G0 = R*P;

% Phi_rf = 180 + arg(G0) == 50 => arg(G0) = -130 deg.
margin(G0);
grid on;

%% -> 20*log10(K) = 18.4 dB

K = 10^(18.4/20)

R = K;
G0 = R*P;

margin(G0); % -> d = 2.67 dB

%% 3.15 (se bara kolku kje bide pr. op. na ZATVORENIOT sistem)
clc
clear

s = tf('s');

P = 200/((s+2)*(s+4)*(s+5))

Gd1 = (s+3)/(s+30);
Gd2 = (s+0.3)/(s+30);

margin (feedback(P, 1));
hold on;

% bez kompenzator -> w0 = 7.17 rad/s

margin(feedback(Gd1 * P, 1));
hold on;

% so Gd1 -> w0 = 3.38 rad/s

margin(feedback(Gd2 * P, 1));

% so Gd2 -> w0 = 12.8 rad/s

grid on;
legend('Gnc', 'Gd1', 'Gd2');

% za Gd1:

z_p1 = 3/30

% za Gd2:

z_p2 = 0.3/30

% so namaluvanje na soodnosot z/p na dif. komp., se zgolemuva propusniot
% opseg na sistemot

P_bw = bandwidth(feedback(P, 1))
P_Gd1_bw = bandwidth(feedback(Gd1 * P, 1))
P_Gd2_bw = bandwidth(feedback(Gd2 * P, 1))

%% 3.24.
clc
clear

s = tf('s');

K = 10;
p1 = 2;

G0 = K/(s^2 * (s+p1));

Gd = (s+0.3)/(s+30);

nyquist(Gd*G0);
hold on;

%% 3.28
clc
clear

syms R C s
Gc = simplify((1/(C*s))/(R + 1/(C*s))); % spec. sluchaj na integralen komp.
pretty(Gc)

s = tf('s');
P = 1/(s+1);

% T = RC = 1s
Gc1 = 1/(s+1);
G01 = P * Gc1;

% T = RC = 10s
Gc2 = 1/(10*s+1);
G02 = P * Gc2;

% T = RC = 100s
Gc3 = 1/(100*s+1);
G03 = P * Gc3;

margin(P);
hold on;
margin(G01);
hold on;
margin(G02);
hold on;
margin(G03);
grid on;
legend('P', 'T = 1s', 'T = 10s', 'T = 100s')

%% Otskochen odziv

step(feedback(P, 1))
hold on;
step(feedback(G01, 1))
hold on;
step(feedback(G02, 1))
hold on;
step(feedback(G03, 1))
grid on;
legend('P', 'T = 1s', 'T = 10s', 'T = 100s')

%% 3.33
clc
clear

s = tf('s');

% uslovi: Mp < 5%, Tp < 30s, Kp >= 6

P = 1000/((s+8)*(s+14)*(s+20))

% R1 e diferencijalen kompenzator -> otpagja
% R3 e obichen zasiulvach -> otpagja

% Proverka za R2 i R4:
R2 = (s+1)/(s + 0.074)
R4 = (20*s+1)/(100*s+1) % so R4 ima mnogu golemo vreme na porast -> otpagja

step(P*R2)
hold on;
step(P*R4)
legend('Gc(s) = R2(s)', 'Gc(s) = R4(s)')
grid on;

% Proverka dali R2 go ispolnuva uslovot Kp >= 6:

Kp = dcgain(P*R2)

% => R2 gi ispolnuva site uslovi