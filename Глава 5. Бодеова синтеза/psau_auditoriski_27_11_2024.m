%% 5.12
clc
clear

s = tf('s');

K = 10^(3.64/20);
R = K;
P = exp(-2*s)/(s+1);

G0 = R*P;

margin(G0);
grid on;

% 20*log10(K) = 3.64 dB => K = 10^(3.64/20)

%% 5.15
clc
clear

s = tf('s');

K = 10000;
G0 = K/((s+5)*(s+20)*(s+50));

margin(G0);
grid on;

% 20*log10(K) = 99.7 dB => K = 10^(99.7/20)
% a) intervalot pri koj sistemot kje ni bide stabilen e K E (0, 96605)

% b) K = 10000

% d = 19.7 dB; Phi_rf = 92.9 stepeni

%% 5.16
clc
clear

s = tf('s');

K = 10^(43/20);
G0 = K/(s*(s+10))

margin(G0);
grid on;

% Phi_rf = 180 + arg(G0) = 45 stepeni => arg(G0) = -135 stepeni.
% 20*log10(K) = 43 dB => K = 10^(43/20)

% maksimalnata dozvolena vrednost na K = 141.3

%% 5.17
clc
clear

s = tf('s');

K = 60.5;
G0 = K/(s*(s+5)^2)

margin(G0);
grid on;

% Phi_rf = 180 + arg(G0) = 45 stepeni => arg(G0) = -135 stepeni.
% 20*log10(K) = 35.5 => K = 10^(35.5/20)

% uslovot e ispolnet za maskimalna vrednost na K = 60.5

%% 5.18
clc
clear

s = tf('s');

K = 20 / 10^(5.74/20);
G0 = K/(s*(s+2)*(s+5))

margin(G0)
grid on;
hold on;

% Phi_rf = 180 + arg(G0) = 55 stepeni => arg(G0) = -125 stepeni.
% 20*log10(K) = 5.74 => K = 10^(5.74/20)

%% 5.18 nachin 2
clc
clear

s = tf('s');

K = 10;
G0 = K/(s*(s+2)*(s+5))

margin(G0)
grid on;
hold on;

% Phi_rf = 180 + arg(G0) = 55 stepeni => arg(G0) = -125 stepeni.
% 20*log10(K) = 20 => K = 10^(20/20)

%% 5.19
clc
clear

s = tf('s');

K = 585;
G0 = K * 100/(s*(s+36)*(s+100))

Mp = 9.5;

C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2)

Phi_rf = atan2d(2*C, sqrt(-2*C^2 + sqrt(1 + 4*C^4)))

margin(G0);
grid on;

% Phi_rf = 59.1621 => arg(G0) = 59.1621 - 180 = -120.8379
% 20*log10(K) = 55.5 dB => K = 10^(55.5/20)

step(feedback(G0, 1))

%% 5.20
clc
clear

s = tf('s');

K = 10^(17.4/20);
G0 = K * 25/(s+25) * 10/(s*(s+10))

Mp = 24;

C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2)

Phi_rf = atan2d(2*C, sqrt(-2*C^2 + sqrt(1 + 4*C^4)))

margin(G0);
grid on;

% Phi_rf = 44.3743 => arg(G0) = 44.3743 - 180 = -135.6257
% 20*log10(K) = 17.4 dB => K = 10^(17.4/20)

step(feedback(G0, 1))

%% Konstanti na greshka so bode

%% Sistem od tip 0
clc
clear

s = tf('s');

r = 0;

p = 3;
z = 5;

G0 = (s+z)/(s^r * (s+p)*(s+1))

bode(G0);
grid on;

Kp_db = 4.44

% 20*log10(Kp) = Kp_db => Kp = 10^(Kp_db/20)

Kp_bode = 10^(Kp_db/20)

Kp = dcgain(G0)

%% Sistem od tip 1
clc
clear

s = tf('s');

r = 1;

p = 2;
z = 6;

G0 = (s+z)/(s^r * (s+p)*(s+1))

freq = 10^-4:0.01:10^4;
mag = 20.*log10(z/p) - (20*r).*log10(freq);

Kv = mean(freq(round(mag, 1) == 0))

plot(freq, mag, 'r--')
hold on;

plot(Kv, mag(round(mag, 1) == 0), 'mo', 'MarkerSize', 5)

bode(G0);
grid on;

Kv = dcgain(s * G0)

%% Sistem od tip 2
clc
clear

s = tf('s');

r = 2;

p = 2;
z = 1;

G0 = (s+z)/(s^r * (s+p))

freq = 10^-3:0.001:10^3;
mag = 20.*log10(z/p) - (20*r).*log10(freq);

Ka = mean(freq(round(mag, 1) == 0))^2

plot(freq, mag, 'r--')
hold on;

plot(sqrt(Ka), mag(round(mag, 1) == 0), 'mo', 'MarkerSize', 5)

bode(G0)
grid on;

Ka = dcgain(s^2 * G0)
