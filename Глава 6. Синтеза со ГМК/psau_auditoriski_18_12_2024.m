%% 6.34 - idealen I-kompenzator
clc
clear

s = tf('s');

G0 = 1/((s+2)*(s+4)*(s+8))

rlocus(G0)
sgrid(0.18, 50)

K_dobieno = 316;

% za da C = 0.18 potrebno e K = 316, a pritoa se dobivaat dominantni polovi
% s12 = -1.01 +- j5.54

%% idealen I-kompenzator -> G(i) = (s+zi)/s; pi=0; z_p = z/0 = inf.
% za da e(inf) = 0, potrebno e da se vovede astatizam od prv red

z = 0.01; % proizvolna vrednost vo blizina na Im-oskata

Gi = (s+z)/s

rlocus(G0)
hold on;
rlocus(G0*Gi)
sgrid(0.18, 50)

pole(feedback(K_dobieno*G0*Gi, 1))

Kp = dcgain(K_dobieno*G0*Gi)

e_stac = 1/(1+Kp)
%% 6.37 - idealen D-kompenzator
clc
clear

s = tf('s');

G0 = 1/(s*(s+3)*(s+5))

rlocus(G0);
sgrid(0.504, 50);

% za ceta = 0.504, dominantniot par na polovi s12 = -0.94 +- j1.61
% pri K = 21.3

K = 21.3

% zd = 0.623

Gd = (s + 0.623)

rlocus(G0 * Gd);
sgrid(0.504, 50);

%% od rlocus za ceta = 0.504 -> Kd = 44.2;

Gd = Gd * 44.2;

pole(feedback(G0*Gd, 1))

step(feedback(G0*Gd, 1))

Gk = zpk(feedback(K*G0, 1))

G0_d = zpk(feedback(G0*Gd, 1))

% bidejkji tretiot pol vo 0.4964 go narushuva povedeniot se aproksimira so
% sistem od vtor red

G_approx = 44.2/(s^2 + 7.504*s + 55.48) % se ponishtuvaat polot vo 0.4964 i nulata vo 0.623

Gk_0 = dcgain(Gk); % Gk(0)
G_approx_0 = dcgain(G_approx); % G_approx(0)

approx_gain = Gk_0/G_approx_0; % Gk(0) = G_approx(0)

G_approx = G_approx * approx_gain

step(Gk);
hold on;
step(G0_d)
hold on;
step(G_approx);
legend('Gk', 'G0 * Gd', 'G_approximated')
hold off;

%% 6.40 - neidealen ID-kompenzator
clc
clear

s = tf('s');

Mp = 20;
C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2) % ceta preku Mp%

G0 = 1/(s*(s+5))

rlocus(G0);
sgrid(C, 50);
ylim([-7 7])

% za ceta = 0.456 dominantnite polovi s1/s = -2.5 +- j4.88
% pri K = 30

K = 30;

%% Neidealen D-kompenzator
clc

Gd = (s+6)/(s+11.35)

rlocus(G0*Gd)
sgrid(C, 50)

%% za ceta = 0.456 Kd = 127, no se zema 113 za da se zachuva Mp% = 20%

Gd = 127 * (s+6)/(s+11.35)

step(feedback(G0*Gd, 1))

%% Neidealen I-kompenzator
clc

Gi = (s+0.1)/(s+0.01)

Kv = dcgain(s*G0*Gd)
e_stac = 1/Kv

Kv_comp = dcgain(s*G0*Gd*Gi)
e_stac_comp = 1/Kv_comp

rlocus(G0*Gd);
hold on
rlocus(G0*Gd*Gi)
sgrid(0.456, 50)
hold off;

%% g) aproksimacija so sistem od vtor red
clc

Gk = zpk(feedback(G0 * K, 1))

G0_id = zpk(feedback(G0 * Gd * Gi, 1))

G_approx = 127/(s^2 + 9.919*s + 119.4) % se ponishtuvaat polot vo 6.34 i nulata vo 6; se ponishtuvaat polot vo 0.1007 i nulata vo 0.1

G0_id_0 = dcgain(G0_id); % G0_id(0)
G_approx_0 = dcgain(G_approx); % G_approx(0)

approx_gain = G0_id_0/G_approx_0; % G0_id(0) = G_approx(0)

G_approx = G_approx * approx_gain

step(Gk);
hold on;
step(G0_id);
hold on;
step(G_approx);
legend('G', 'G0 * Gid', 'G_approximated')
hold off;
