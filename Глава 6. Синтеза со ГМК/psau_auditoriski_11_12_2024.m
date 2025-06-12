%% Potsetuvanje za GMK
clc
clear

s = tf('s');

G0 = (s+1)/((s+2)*(s+3))

rlocus(G0) % naredba za iscrtuvanje na GMK

% GMK se skicira za otvoreniot sistem pri zasiluvanje K = 1 !

%% 6.29
clc
clear

s = tf('s');

G0 = (s+2)/(s^2 + 4*s + 5)

rlocus(G0)

K = 2;
pole(feedback(K*G0, 1))

%% 6.28
clc
clear

s = tf('s');

G0 = 1/((s+1)*(s+15)^2)

Kp_new = 20;

Kp = dcgain(G0)

K = Kp_new/Kp

%% za Kp >= 20 potrebno e K >= 4500 (uslov I)
clc

% K_gr -> vrednost na K za koja sistemot se naogja na gr. na stabilnost
% K_konstruktivno -> potrebnata vrednost na K za da se dobie baranata rez. na zasiluvanje

% rezervata na zasiluvanje d e definirana kako: d = K_gr/K_konstr
% od tuka sleduva deka K_konstruktivno = K_gr/d

d_new = 1.5;

rlocus(G0)

%% sistemot se naogja na gr. na stabilnost pri K_gr ~= 7680

K_konstr = 7680/d_new

% za d >= 1.5 potrebno e K <= 5120 (uslov 2)
% od tuka sleduva 4500 <= K <= 5120

%% proverka

K = 4750
Kp = dcgain(K*G0)

[Gm, Pm, Wcg, Wcp] = margin(K*G0);

d = Gm

%% 6.25
clc
clear

s = tf('s');

G0 = 1/(s*(s^2 + 2*s + 5))

rlocus(G0) % sleduva deka pri a = 10 sistemot e na gr. na stabilnost

%% 6.23
clc
clear

s = tf('s');

G0 = (s^2 + 10*s + 74)/(s^2 * (s+4))

rlocus(G0)
sgrid(0.512, 50) % za wn vrednosta e proizvolna

% za K = 49.5 se dobiva vrednost na ceta = 0.512
% dominatnite polovi na sistemot kje bidat -4.66 +- j7.82

%% 6.24
clc
clear

syms K s G0

G = K/((s+25)^2 + K);

G0 = solve(G == G0/(1+G0), G0);

% G0 se dobiva deka e G0 = K/(s+25)^2

s = tf('s');

G0 = 1/(s+25)^2

rlocus(G0)
sgrid(0.707, 50)
ylim([-30 30])

% za ceta = 0.707 potrebno e K da bide K = 625 so shto se dobivaat d.p.
% s_1/2 = -25 +- j25

K = 625;
pole(feedback(G0*K, 1))

%% 6.30 - neidealen I-kompenzator
clc
clear

s = tf('s');

G0 = 20/(s*(s+10)^2)

Kv = dcgain(s*G0)
e_stac = 1/Kv

e_stac_new = e_stac/100
Kv_new = 1/e_stac_new

gain_needed = Kv_new/Kv % od tuka go naogjame potrebnoto zasiluvanje

% gain_needed e baraniot soodnost pomegju nulata i polot z_p na komp.

% Gi(s) = (s+z)/(s+p); z > p

z_p = gain_needed;

z = 0.1; % proizvolna vrednost vo blizina na imaginarnata oska
p = z * 1/z_p; % se presmetuva spored potrebniot soodnos i proizvolnata nula

Gi = (s+z)/(s+p)

rlocus(G0)
hold on;
rlocus(G0*Gi) % mozhe da se uvidi deka GMK e skoro nepromeneto

%% proverka
clc

Kv = dcgain(s*G0*Gi)

e_stac = 1/Kv

%% 6.32 - neidealen I-kompenzator
clc
clear

s = tf('s');
G0 = 1/(s*(s+5))

% pri sinteza so GMK sekogash prvo se nagoduva PREODEN odziv, a potoa
% STACIONAREN

rlocus(G0)

ceta = 0.592;
wn = 50;

sgrid(ceta, wn)
ylim([-7, 7])

K_dobieno = 17.8

% za da C = 0.592 potrebno e K da bide 17.8, a pritoa se dobivaat
% dominantni polovi s12 = -2.5 +- j3.4

%%
clc

Kv = dcgain(s*K_dobieno*G0)
e_stac = 1/Kv

% za presmetanata vrednost na K = 17.8 se dobiva ev(inf) = 28%, shto e
% mnogu pogolema od baranata vrednost od ev(inf) = 2%

e_stac_posakuvano = 0.02

e_stac_dobieno = 5/K_dobieno

gain_needed = e_stac_dobieno/e_stac_posakuvano

%% gain_needed e vsushnost potrebniot odnos z/p na kompenzatorot

z_p = gain_needed

z = 0.1; % proizvolna vrednost vo blizina na Im-oskata
p = z * 1/(z_p);

Gi = (s+z)/(s+p);

rlocus(G0)
hold on;
rlocus(G0*Gi)
sgrid(0.592, 50)
ylim([-7, 7])

%% Proverka
clc

Kv_new = dcgain(s*K_dobieno*G0*Gi)
e_stac_new = 1/Kv_new

%% Sporedba na odzivi
K = 17.8
G0 = 1/(s*(s+5));
Gi = (s+z)/(s+p)

t = 0:0.01:100;

lsim(feedback(G0, 1), t, t);
hold on;
lsim(feedback(K*G0, 1), t, t);
hold on;
lsim(feedback(K*G0*Gi, 1), t, t);
legend('G0', 'K*G0', 'K*G0*Gi')