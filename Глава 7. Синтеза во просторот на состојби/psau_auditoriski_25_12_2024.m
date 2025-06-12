%% Primer 1
clc
clear

Mp = 20.8;
C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2)
wn = 2.23

s12 = -C*wn + sqrt(-1)*wn*sqrt(1-C^2)

A = [-5 1 0; 0 -2 1; 0 0 -1]
B = [0 0 1]'
C = [-1 1 0]
D = [0]

sys = ss(A,B,C,D)

[num, den] = ss2tf(A, B, C, D)

G_nc = tf(num, den)

syms s k1 k2 k3

as = collect((s+1-1j*2)*(s+1+1j*2)*(s+5)) % se dodava tret pol 5 do 10 pati polevo od ostanatite za da ima pomalo vlijanie

ak = collect(det([s+5 -1 0; 0 s+2 -1; k1 k2 s+1+k3]))

sol = solve(coeffs(ak, s) == coeffs(as, s), [k1 k2 k3]);

%% Sporedba na odzivi
clc

k = double([sol.k1 sol.k2 sol.k3])
A_c = A - B*k

[num_c, den_c] = ss2tf(A_c, B, C, D)

G_c = tf(num_c, den_c)

pole(G_c)

step(G_nc)
hold on;
step(G_c)
legend('Nekompenziran', 'Kompenziran')

%% Primer 2
clc
clear

Mp = 20.8;
C = (-log(Mp/100))/sqrt(pi^2 + (log(Mp/100))^2)
wn = 2.23

s12 = -C*wn + sqrt(-1)*wn*sqrt(1-C^2)

A = [-5 1 0; 0 -2 1; 0 0 -1]
B = [0 0 1]'
C = [-1 1 0]
D = [0]

sys = ss(A,B,C,D)

[num, den] = ss2tf(A, B, C, D)

G_nc = tf(num, den)

syms s kw1 kw2 kw3

Aw = [0 1 0; 0 0 1; -10 -17 -8];
Bw = [0 0 1]';

kw = [kw1 kw2 kw3];

as = collect((s+1-1j*2)*(s+1+1j*2)*(s+5))

ak = collect(det(s*eye(3) - (Aw - Bw*kw)))

sol = solve(coeffs(ak, s) == coeffs(as, s), [kw1 kw2 kw3]);

kw = double([sol.kw1 sol.kw2 sol.kw3])

%% Proverka za upravlivost
clc

Qcv = ctrb(A, B)

fprintf('Proverka za upravlivost:\n')
if rank(Qcv) == rank(A)
   fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => originalniot sistem e kompletno upravliv\n', rank(Qcv), length(A), length(A))
else
    fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => originalniot sistem NE e kompletno upravliv\n', rank(Qcv), length(A), length(A))
end

Qcw = ctrb(Aw, Bw)

fprintf('Proverka za upravlivost:\n')
if rank(Qcw) == length(Aw)
   fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => transformiraniot sistem e kompletno upravliv\n', rank(Qcv), length(A), length(A))
else
    fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => transformiraniot sistem NE e kompletno upravliv\n', rank(Qcv), length(A), length(A))
end

%% Opredeluvanje na P

P = Qcv * inv(Qcw)

kv = kw * inv(P) % se poklopuva so dobienoto k od primer 1!

%% Primer 3
clc
clear

syms s

as = collect((s+1)^2 * (s+2))
aso = collect((s+1*10)^2 * (s+2*10))

A = [0 1 0; 0 0 1; -4 -6 -4];
B = [0 0 1]';
C = [1 0 1];
D = 0;

% Opredeluvanje na zasiluvanja na sistemot {matrica K}
syms s k1 k2 k3

K = [k1 k2 k3]

ak = collect(det(s*eye(3) - (A - B*K)))

sol = solve(coeffs(ak, s) == coeffs(as, s), [k1 k2 k3]);

K = double([sol.k1 sol.k2 sol.k3])

% Opredeluvanje na zasiluvanja na observerot {matrica L}
Aw = A.';
Bw = C.';
Cw = B.';
Dw = D;

syms Lw1 Lw2 Lw3

Lw = [Lw1 Lw2 Lw3].';

aLw = collect(det(s*eye(3) - (Aw - Lw*Cw)))

sol = solve(coeffs(aLw, s) == coeffs(aso, s), [Lw1 Lw2 Lw3]);

Lw = double([sol.Lw1 sol.Lw2 sol.Lw3])'

%% Proverka za nabljudlivost
clc

Qov = obsv(A, C)

fprintf('Proverka za nabljudlivost:\n')
if rank(Qov) == rank(A)
   fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => originalniot sistem e kompletno nabljudliv\n', rank(Qov), length(A), length(A))
else
    fprintf('Rangot na Qcv e %d, dimenziite na A se %dx%d => originalniot sistem NE e kompletno nabljudliv\n', rank(Qov), length(A), length(A))
end

Qow = obsv(Aw, Cw)

fprintf('Proverka za nabljudlivost:\n')
if rank(Qow) == length(Aw)
   fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => transformiraniot sistem e kompletno nabljudliv\n', rank(Qow), length(A), length(A))
else
    fprintf('Rangot na Qcw e %d, dimenziite na A se %dx%d => transformiraniot sistem NE e kompletno nabljudliv\n', rank(Qow), length(A), length(A))
end

%% Opredeluvanje na P

P = inv(Qov) * Qow

Lv = P * Lw

K = [-2 -1 0];

polovi_observer = eig(A - Lv*C)
polovi_sistem = eig(A - B*K)