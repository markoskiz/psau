Set-ExecutionPolicy Bypass -Scope Process -Force; `
[System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072; `
iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))






choco install lanmessenger -y

choco uninstall lanmessenger -y


posle WinKey + R -> pisi %appdata% izbrisi sve od Lan ms app
num = 5; 
den = [15 32 15];
G = tf(num, den);

%% edinecen vlez 
[y, t] = step(G, 25);
figure;
plot(t, y, 'b', 'LineWidth', 1.5); hold on;
title('odziv');
xlabel('Vreme (s)');
ylabel('Izlez');
grid on;

% gain Kp
Kp = y(end);

%% Tu
dy = diff(y) ./ diff(t);
[~, idx] = max(dy);

t_inflection = t(idx);
y_inflection = y(idx);
slope = dy(idx);

Tu = t_inflection - y_inflection / slope;

%% Tg
Tg = (Kp - y_inflection) / slope;

%% 2. FOPDT aproksimacija
G_fopdt = tf(Kp, [Tg 1], 'InputDelay', Tu);


[y_fopdt, t_fopdt] = step(G_fopdt, 25);
plot(t_fopdt, y_fopdt, 'r--', 'LineWidth', 1.5); % overlay

legend('originalen sistem', 'FOPDT aproksimacija');
