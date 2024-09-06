% Parâmetros do sistema
m = 10; % Massa
zeta_sub = 0.4; % Amortecimento subcrítico
omega_d = 1; % Frequência natural amortecida (1 rad/s)

% Cálculo de k e b
omega_n = omega_d / sqrt(1 - zeta_sub^2); % Frequência natural não amortecida
k = m * omega_n^2; % Constante da mola
b_sub = 2 * zeta_sub * sqrt(m * k); % Subamortecido

% Cálculo de b crítico
b_crit = 2 * sqrt(m * k); % Criticamente amortecido

% Entrada degrau
u = 1;

% Função do sistema massa-mola-amortecedor
function dx = mass_spring_damper(t, x, u, b, m, k)
    dx = zeros(2, 1);
    dx(1) = x(2); % Velocidade
    dx(2) = (u - (b/m) * x(2) - (k/m) * x(1)); % Aceleração com entrada degrau
end

% Condições iniciais
x0 = [0; 0]; % Posição inicial = 0 m, Velocidade inicial = 0 m/s

% Simulação
t = 0:0.01:10;
[t, x] = ode45(@(t, x) mass_spring_damper(t, x, u, b_sub, m, k), t, x0);

% Plot
plot(t, x(:, 1), 'r', t, x(:, 2), 'b');
xlabel('Tempo (s)');
ylabel('Posição (m) e Velocidade (m/s)');
legend('Posição', 'Velocidade');
