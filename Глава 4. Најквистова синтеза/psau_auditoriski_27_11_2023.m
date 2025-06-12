%% 4.1
clc
clear

s = tf('s');

G0 = 2/(s+1);

nyq = nyquistplot(G0)
setoptions(nyq, 'ShowFullContour', 'off')
hold on;

x = 0:0.01:2*pi;
plot(sin(x), cos(x), 'k-.', 'LineWidth', 2);
hold on;

xt = 0.5;
yt = -0.865;

plot(xt, yt, 'ro', 'LineWidth', 2)
hold on;

line([0, xt], [0, yt], 'color', 'c', 'LineWidth', 2);
hold on;

line([xt, xt], [yt 0], 'color', 'm', 'LineWidth', 2);
line([0 xt], [0 0], 'color', 'm', 'LineWidth', 2);

argG0 = atan2d(yt, xt)

Phir = 180 + argG0
%% 4.9
clc
clear

s = tf('s');

K = 1;

G0 = (9*K)/((s+1)*(s^2 + 3*s + 9));

nyq = nyquistplot(G0);
setoptions(nyq, 'ShowFullContour', 'off');
hold on;

% 20*log10(K) = 12.7Db => K = 10^(12.7/20)

K = 10^(12.7/20)

%% 4.11
clc
clear

s = tf('s');

K = 1;

G0 = (K*exp(-0.1*s))/(s+5);

nyq = nyquistplot(G0);
setoptions(nyq, 'ShowFullContour', 'off');
hold on;

%% 4.13
clc
clear

syms s

P = 1/(s*(1 + s/8)*(1+s/20));

R1 = (1+100*s)/(1+120*s);
R2 = ((1+s)*(1+20*s))/((1+s/0.01)*(1+s/50));
R3 = 100;
R4 = (100*(1+s)*(1+s/5))/((1+s/0.1)*(1+s/50));

% 1) Kv >= 100
% Kv = lim(s->0) [s * G0(s)]

G01 = P * R1;
G02 = P * R2;
G03 = P * R3;
G04 = P * R4;

Kv1 = limit(s*G01, s, 0);
Kv2 = limit(s*G02, s, 0);
Kv3 = limit(s*G03, s, 0);
Kv4 = limit(s*G04, s, 0);

Kv_vals = [Kv1 Kv2 Kv3 Kv4] % => samo R3 i R4 go zadovoluvaat baranjeto 1)

% 2) d >= 10 dB

s = tf('s');

P = 1/(s*(1 + s/8)*(1+s/20));
R3 = 100;
R4 = (100*(1+s)*(1+s/5))/((1+s/0.1)*(1+s/50));

G03 = P * R3;
G04 = P * R4;

nyq = nyquistplot(P);
hold on;
nyq = nyquistplot(G03);
hold on;
nyq = nyquistplot(G04);
setoptions(nyq, 'ShowFullContour', 'off');
legend('P bez komp.', 'komp. so R3', 'komp. so R4');
hold off;

%% 4.18
clc
clear

s = tf('s');

K = 100;
G0 = K/(s*(s+1)*(s+5));

Gd = (s+1)/(s+10)
Gi = (s+10)/(s+1)

nyq = nyquistplot(G0);
setoptions(nyq, 'ShowFullContour', 'off');
hold on;
nyquistplot(G0 * 10)
hold on;
nyquistplot(G0 * Gd)
hold on;
nyquistplot(G0 * Gi)
legend('G0 bez komp.', 'komp. so zasiluvanje', 'komp. so Gd', 'komp. so Gi');
hold off;
%% 4.20
clc
clear

s = tf('s');

P = 4/(s*(s^2 + 3.2*s + 64))

nyq = nyquistplot(P);
setoptions(nyq, 'ShowFullContour', 'off');
hold on;

z = 10;
p = 1;

K = p/z;
Gc = 1 * (s+1.4)/(s+0.1)
nyquistplot(Gc * P);
legend('Nekompenziran', 'Kompenziran');
hold off;

%% Bodeovi dijagrami
clc

margin(P);
hold on;
margin(Gc*P);
grid on;
hold off;

%% Otskochen odziv
clc

step(feedback(P, 1));
hold on;
step(feedback(Gc*P, 1));
hold off;

%% 4.24
clc
clear

z = tf('z');

K = [1, 0.01, 0.05, 0.1, 0.05, 10, 50];

for i=1:1:length(K)
    G0 = (K(i)*(z+1)^3)/((z-1)*(z+0.25)*(z+0.6));
    nyq = nyquistplot(G0);
    hold on;
end

setoptions(nyq, 'ShowFullContour', 'off');
legend('K = 1', 'K = 0.01', 'K = 0.05', 'K = 0.1', 'K = 0.5', 'K = 10', 'K = 50')