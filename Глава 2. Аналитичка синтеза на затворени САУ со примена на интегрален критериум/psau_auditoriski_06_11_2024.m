%% 2.12
clc
clear

syms K1 K2 s

P = 1/(4*s^3 + 9*s^2 + 2*s + 1);

disp('P(s)=');
pretty(P)

R = K1 + K2/s;

disp('R(s)=');
pretty(R)

N = 1/s;

% Ey = 1/(1+G0) * Y <=> [1 - Gy] * Y
% En = -Gn * N <=> -P/(1 + R*P) * N

En = simplify(-P/(1 + R*P) * N);

disp('E(s)=');
pretty(collect(En))

a4 = 4;
a3 = 9;
a2 = 2;
a1 = K1 + 1;
a0 = K2;

b3 = 0;
b2 = 0;
b1 = 0;
b0 = -1;

J = simplify((b3^2*(a0*a1*a2 - a0^2*a3) + a0*a1*a4*(b2^2 - 2*b1*b3) + a0*a3*a4*(b1^2 - 2*b0*b2))/(2*a0*a4*(a1*a2*a3 - a0*a3^2 - a1^2*a4)) + ((a2*a3*a4 - a1*a4^2)*b0^2)/(2*a0*a4*(a1*a2*a3 - a0*a3^2 - a1^2*a4)));

disp('J=');
pretty(J)

dJK1 = simplify(diff(J, K1));
dJK2 = simplify(diff(J, K2));

disp('dJK1=');
pretty(dJK1)

disp('dJK2=');
pretty(dJK2)

solution = solve([dJK1 == 0, dJK2 == 0], [K1 K2]);

sol_1 = [solution.K1(1), solution.K2(1)]
sol_2 = [solution.K1(2), solution.K2(2)]

% se dobivaat tochkite (7/2, 0) i (2, 1/9)

%% Parametarska oblast na stabilnost
clc

[b, a] = numden(En);

disp('a(s)=')
pretty(collect(a))

% Routh-Hurwitz

% a4 a2 a0
% a3 a1 a0
% B1 B2  0
% C1  0  0
% D1  0  0

B1 = simplify(-det([a4 a2; a3 a1])/a3);
B2 = simplify(-det([a4 a0; a3 0])/a3);

C1 = simplify(-det([a3 a1; B1 B2])/B1);

D1 = simplify(-det([B1 B2; C1 0])/C1);

routh = [a4 a2 a0; a3 a1 a0; B1 B2 0; C1 0 0; D1 0 0];

disp('Routh Table:');
pretty(routh);

r_conditions = simplify(routh > 0);

disp('Routh Conditions:');
pretty(r_conditions);

cond_K2(K1) = solve(r_conditions(4, 1), K2)

% -1 < K1 < 7/2
% K2 < (4*K1)/27 - (4*K1^2)/81 + 7/81

%% Plotiranje na zavisnosta pomegju K1 i K2
clc

cond_K1 = -1:0.01:7/2;

cond_K2r = cond_K2(cond_K1);

plot(cond_K1, cond_K2r, 'color', '#0072BD', 'LineWidth', 3);
hold on;
xline(-1, '-..', 'color', '#77AC30', 'LineWidth', 3);
hold on;
xline(7/2, '-..', 'color', '#77AC30', 'LineWidth', 3);
hold on;
plot(sol_1(1), sol_1(2), '.', 'color', '#7E2F8E', 'MarkerSize', 30);
plot(sol_2(1), sol_2(2), '.', 'color', '#7E2F8E', 'MarkerSize', 30);
hold on;
area(cond_K1(cond_K2r > 0), cond_K2r(cond_K2r > 0), 'FaceColor', '#EDB120', 'FaceAlpha', 0.25)
axis([-1.5 4 -0.05 0.2]);
grid on;
