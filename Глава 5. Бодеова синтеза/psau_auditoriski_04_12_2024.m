%% 5.22
clc
clear

syms s K

G0 = K/(s*(s+10))

% a) e(inf) < 5% => e(inf) < 0.05
% b) Phirf > 45 deg.

% a)

Kv = limit(s*G0, s, 0)

e_stac = 1/Kv

sol = solve(e_stac <= 0.05, 'ReturnConditions', true);

sol.conditions

% K mora da e >= 200 za da bide ispolneto baranjeto pod a)

s = tf('s');

K = 200;
G0 = K/(s*(s+10))
%%
margin(G0);
grid on;

%% se proektira dif. kompenzator Gd(s) = p/z * (s+z)/(s+p); z < p

Phirf_wanted = 45;
Phirf_current = 38.7;
corr_factor = 5; % korekcionen faktor so cel da se nadomesti slabeenjeto od komp., mozhe da bide od 5 do 12 stepeni

% maksimalnata pozitivna faza shto ja vnesuva diff. komp.
Phimax = Phirf_wanted - Phirf_current + corr_factor

% soodnos nula/pol -> z/p
z_p = (1 - sind(Phimax))/(1 + sind(Phimax))

% magnitudata na amplitudniot dijagram pri frekvencija wmax
mag_wmax = 20*log10(1/sqrt(z_p))

% wmax se opredeluva od bodeovite dijagrami (preku mag_wmax)
wmax = 11;

% wmax = sqrt(p*z) => sqrt(p*z) * sqrt(z/p) = z;
% sqrt(p*z) * 1/sqrt(z/p) = sqrt(p*z) * sqrt(p/z) = p

z = wmax*sqrt(z_p);
p = z/z_p;

Gd = p/z * (s+z)/(s+p);

bode(G0, G0*Gd, Gd)
legend('G0_{nc}', 'G0*Gd', 'Gd')
grid on;

%% proverka za a)

Kv = dcgain(G0*Gd*s)

e_stac = 1/Kv

%% 5.23
clc
clear

s = tf('s');
G0 = 72/((s+1)*(s+3)^2)

% a) ist propusen opseg w0/ista presechna frekvencija w1 => Kd = sqrt(p/z)
% b) Phirf ~= 45

margin(G0);
grid on;

% treba da se zgolemi Phi_rf do 45 a pritoa w1 da ostane 3.38 rad/s
%% se proektira dif. kompenzator od II tip
wmax = 3.38;

Phirf_wanted = 45;
Phirf_current = 9.66;
corr_factor = 0; % od 5 do 12

Phimax = Phirf_wanted - Phirf_current + corr_factor;

z_p = (1 - sind(Phimax))/(1 + sind(Phimax))

z = wmax*sqrt(z_p);
p = z/z_p;

Gd = sqrt(1/z_p) * (s+z)/(s+p)

bode(G0, G0*Gd, Gd)
legend('G0_{nc}', 'G0*Gd', 'Gd')
grid on;

%% 5.27
clc
clear

s = tf('s');

K = 1;
G0 = K * 100/(s*(s+50)*(s+100))

Mp = 10;

C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2) % ceta preku Mp%

Phi_rf = 100*C % aproksimacija na Phirf preku ceta

% Phirf = 180 + arg(G0) = 59.115 => arg(G0) = -120.885

margin(G0)
grid on;

%% za da Phi_rf = 59.115 => K = 10^(60/20)

K = 10^(59.4/20);
G0 = K * 100/(s*(s+50)*(s+100))

Kv = dcgain(s*G0)

e_stac = double(1/Kv)

e_stac_komp = e_stac/10

Kv_new = 1/e_stac_komp % ova ni ukazhuva deka vo sistemot treba da dodademe zasiluvanje od 10 (Kv_new e 10 pati pogolemo od Kv)

%%
K = 10^(59.4/20) * 10;
G0 = K * 100/(s*(s+50)*(s+100))

margin(G0);
grid on;

% bidejkji sistemot preminuva vo nestabilen, potrebno e da se proektira
% integralen komp. Gi(s) = p/z * (s+z)/(s+p); z > p

%%
corr_factor = 2;

wmax = 15.7 % so korekcionen faktor od 2 stepeni, preku arg(G0) (od prethodno)
mag_wmax = 21;

z_p = 10^(mag_wmax/20);

z = wmax/10; % po konvencija, mozhe da bide 10, 15 itn.
p = z/z_p % z * 1/z/p = z * p/z = p

Gi = 1/z_p * (s+z)/(s+p)

bode(G0, G0*Gi, Gi)
legend('G0_{nc}', 'G0*Gi', 'Gi')
grid on;

%% proverka za e(inf)
clc

K = 10^(59.4/20) * 10;
G0 = K * 100/(s*(s+50)*(s+100))

Gi = p/z * (s+z)/(s+p)

Kv = dcgain(s*G0*Gi)
e_stac = double(1/Kv)

%% 5.29
clc
clear

syms s K

G0 = K/(s*(s+5)^2)

Kv = limit(s*G0, s, 0)

K_sol = double(solve(Kv > 5))

s = tf('s');

K = K_sol
G0 = K/(s*(s+5)^2)

margin(G0);
grid on;

%%
wmax = 1;

z = wmax/10; % po konvencija
z_p = 10^(13.8/20);
p = z/z_p;

Gi = p/z * (s+z)/(s+p);

bode(G0, G0*Gi, Gi)
legend('G0_{nc}', 'G0*Gi', 'Gi')
grid on;

%% proverka za Kv

Kv = dcgain(s*G0*Gi)

%% 5.30
clc
clear

syms s K
G0 = K/(s*(s+8)*(s+30))

Mp = 10;
Tr = 0.6;

C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2) % ceta preku Mp%

Phi_rf = 100*C % Phirf preku ceta

% od Phi_rf zakluchuvame deka arg(G0) ~= -121 deg.

w0 = (2.2/(Tr * C)) * sqrt((1 - 2*C^2) + sqrt(4*C^4 - 4*C^2 + 2)) % propusen opseg preku vreme na porast Tr i ceta

w1 = w0 % se usvojuva w1 <= w0 <= 2*w1, najchesto se zema w1 = 0.8*w0 !

Kv = limit(s*G0, s, 0);

K_sol = double(solve(Kv == 10))

s = tf('s');

K = K_sol;
G0 = K/(s*(s+8)*(s+30));

margin(G0)
grid on;
%% Gdi(s) = Gd(s) * Gi(s) = (s + zd)/(s + pd) * (s + zi)/(s + pi)
clc

Phirf_wanted = Phi_rf;
Phirf_current = 180 - 146; % 180 - arg[G0(w1)] -> arg(G0(w1)) = -146 -> se bara preku novoto w1 dobieno od w0
corr_factor = 12; % od 5 do 12

Phimax = Phirf_wanted - Phirf_current + corr_factor;

B = (1 - sind(Phimax))/(1 + sind(Phimax))
Y = 1/B

% prvo se proektira Gi(s)
zi = w1/10;
pi = zi/Y;

Gi = (1/Y) * (s + zi)/(s + pi);

% se proektira Gd(s)
zd = w1 * sqrt(1/Y);
pd = zd * Y;

Gd = Y * (s+zd)/(s+pd);

Gid = Gi * Gd;

bode(G0, G0*Gid, Gid)
legend('G0_{nc}', 'G0*Gid', 'Gid')
grid on;

%% proverka

Kv = dcgain(s*G0*Gid)
step(feedback(G0, 1), feedback(G0*Gid, 1))

% paramtetar | proektno baranje | realizirano
%     Kv              10               10
%     Mp              10%              8.38%
%     Tp              0.6s             0.527s