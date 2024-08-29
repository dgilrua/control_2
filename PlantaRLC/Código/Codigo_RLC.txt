%% Graficación del sistema con solución a las E.D.O's con "ode45"
clc; format long

% Define timespan
tspan = [0 0.8];

% Define initial conditions
x0 = [0; 0];

% Solve the system of equations
[t, x] = ode45(@system, tspan, x0);

% Plot the solutions
plot(t, x(:,2), 'r', t, x(:,1),'b')
legend('x2: Voltaje en el Capacitor','Datos reales'); xlabel('Tiempo'); ylabel('Voltaje')

%% Modelo del sistema

% Parámetros del circuito
R1=1000;     % [Ohm]
R2=220;      % [Ohm]
R3=1000;     % [Ohm]
L=1000e-6;   % [mHenrrios]
C=470e-6;    % [Faradios]
V=4;         % [Voltios] 

% Matrices
A=[-R2/L -1/L; 1/C (1/C)*((-1/R1)-(1/R3))]; % [iL ; Vc]
B=[1/L; 1/(R1*C)];
CM=[0 1];
D=0;

% Create state space model
sys = ss(A,B,CM,0);
Ts = 7e-3;

% Convert to discrete, where dt is the discrete time-step (in seconds)
d_sys = c2d(sys, Ts,'zoh');
[Ad, Bd, Cd, Dd] = ssdata(d_sys);

%% Transfer Function
clc;

[numTF, denTF] = ss2tf(A,B,CM,D);
Transfer_Fuction = tf(numTF,denTF);

%% ESTABILIDAD
clc; syms s

% ESTABILIDAD EXTERNA
disp('Estabilidad externa (T.F)');
Polos1 = pole(Transfer_Fuction);
raices = solve(Polos1, s);
disp('Los polos del sistema son:');
disp(Polos1);
fprintf('\n')

% ESTABILIDAD INTERNA - Método del determinante
disp('Estabilidad interna - Método del determinante');
Polos2 = det(s*eye(size(A)) - A);
disp(['El polinomio es:', ' ', char(Polos2)]);
raices = solve(Polos2, s);
disp(['Los polos del sistema son:', ' ', char(raices)]);

%% CONTROLABILIDAD
clc;

% Análisis de Controlabilidad en la entrada
fprintf('CONTROLABILIDAD EN LA ENTRADA')
Controlabilidad2 = ctrb(A,B);
Tamano2 = size(Controlabilidad2);   
Dimesion2 = rank(Controlabilidad2);
disp('Por el método del comando "ctrb()":');
disp(['El tamaño de la matriz es ','[', num2str(Tamano2),']', ' y el rango es ', num2str(Dimesion2)]);
% Para que el sistema sea completamente controlable la matriz de controlabilidad debe ser de rango completo
%% Datos Arduino - Lazo abierto

%Importación de datos reales
data_x = datos(:, 1);        %Selección la primera columna
DatosReales = data_x(1:80);     % Crea una variable para almacenar los datos

%% Gráficas
close all;

t2 = linspace(0, 0.4, length(DatosReales));
x0=[0 0];
[y,t1,x1] = lsim(sys, 4*ones(size(t2)),t2,x0);

% Graficar la respuesta al impulso del sistema discreto
figure(1)
plot(t1,y,'r')
hold on
step(4*d_sys,'b');
%plot(t, x(:,2), 'r')
legend('Continuo', 'Discreto')
title('Respuesta al Impulso del Sistema Discreto (ZOH)');
xlabel('Muestras');
ylabel('Voltaje [V]');
grid on

figure(2)
plot(t2, DatosReales, 'g')
legend('Datos reales')
title('Datos del sistema real');
ylabel('Voltaje [V]');
grid on

figure(3)
step(4*d_sys,'b');
hold on
plot(t2, DatosReales, 'g')
legend('Sist. Discreto', 'Datos reales')
title('Sist. Discreto vs Datos reales');
ylabel('Voltaje [V]');
grid on

%% ASIGNACIÓN DE POLOS
close all; clc;

%Condiciones de diseño
mp = 8;       % sobrenivel
ts = 0.5;      % tiempo de establecimiento
z = sqrt(log(mp)^2)/sqrt(pi^2 + log(mp)^2); % Amortiguamiento
wn = 4/(ts*z);                              % Frecuencia
denominador = [1 2*wn*z wn^2];
p = roots(denominador);

% Discretización de los polos deseados
pd= exp(p.*Ts);

% ASIGNACIÓN DE POLOS
k = acker(Ad,Bd,pd); % Vector de ganacias para la retroalimentación de estado

% Valor de referencia
reff=3;

% Nuevas matrices del sistema en lazo cerrado (controlado)
% f: Feedforward (Retroalimentación)
Af= Ad-(Bd*k);
H = inv(Cd*inv(eye(2)-Af)*Bd);   % Precompensación
Bf = Bd*H*reff;
sysDiscreto = ss(Af,Bf,CM,D,Ts);

%% GRÁFICAS - CONTROL POR RETRO DE ESTADO & PRECOMPENSACIÓN
close all; clc;

figure(1)
step(sysDiscreto)
title('Sistema discretizado en lazo cerrado (Controlado)')
ylabel('Voltaje en el Capacitor [V]'); xlabel('Tiempo')
legend('Sistema con controlador')


%Importación de datos reales
data_x1 = datosRetro_Precom(:, 1);        %Selección la primera columna
DatosReales2 = data_x1(1:80);     % Crea una variable para almacenar los datos

t3 = linspace(0, 0.7, length(DatosReales2));
x0=[0 0];
[y,t4,x1] = lsim(sys, 3*ones(size(t3)),t3,x0);

figure(2)
step(sysDiscreto)
hold on
plot(t4, DatosReales2, 'g')
grid on
title('Sistema discretizado en lazo cerrado (Controlado)')
ylabel('Voltaje en el Capacitor [V]'); xlabel('Tiempo')
legend('Sist. Discreto', 'Datos reales')
xlim([0 0.7]);

%% Esfuerzo de control
% Esfuerzo de control: Voltaje y Corriente dados por el DAC
clear up x
close all; clc;

x(:,1)=[0 0]';
for i=1:120
    up(i)=-k*x(:,i)+H*3;
    x(:,i+1)=Ad*x(:,i)+Bd*up(i);
end

figure (1)
stairs(up)
title('Esfuerzo de control - Sistema discretizado')
ylabel('Voltaje entregado por el DAC [V]'); xlabel('Tiempo')
legend('Sistema con controlador')
xlimit=[0 120];

figure (2)
stairs(x(1,:))
title('Esfuerzo de control - Sistema discretizado')
ylabel('Corriente en el circuito [A]'); xlabel('Tiempo')
legend('Sistema con controlador')

figure (3)
stairs(x(2,:))
title('Esfuerzo de control - Sistema discretizado')
ylabel('Voltaje en el Capacitor [V]'); xlabel('Tiempo')
legend('Sistema con controlador')

%% LAZO ABIERTO VS LAZO CERRADO

% Vector de tiempo
t = 0:Ts:0.9;
x0=[0 0];

u1=ones(size(t));
%u1(60:80)=1.8;
u1(60:end)=1.3;

u2=1.5*ones(size(t));
%u2(60:80)=1.8;
u2(60:end)=1.8;


[y1, t1] = lsim(sysDiscreto,u1,t,x0);
[y2, t2] = lsim(d_sys, u2,t,x0);

% Gráfico
figure;
hold on;
plot(t1,y1,'r');
plot(t2,y2,'b');
title('Voltaje en el Capacitor - Sistema discretizado');
ylabel('Voltaje [V]'); xlabel('Tiempo')
legend('Sistema con controlador','Sistema en lazo abierto');
%% ACCIÓN INTEGRAL
clc; close all;

pNuevo=[p(1),p(2),10*p(1)];                         
poloDiscreto_AccionInt = exp(real(pNuevo).*Ts);


Af_Intg = [Ad [0;0];-Cd 1];
Bf_Intg = [Bd; 0];

% Vector de ganacias
K2 = acker(Af_Intg,Bf_Intg,poloDiscreto_AccionInt);
K_Retro = K2(1:length(K2)-1); % Ganacia de Retroalimentación
K_Intg = K2(length(K2));      % Ganacia Integral

% Matrices aumentadas discretas en lazo cerrado
Af_Exp = [Af_Intg - Bf_Intg*K2];
Bf_Exp = [zeros(rank(A),1);1];
Cf_Exp = [Cd 0];

% Space State
sysDiscreto_Int = ss(Af_Exp,Bf_Exp,Cf_Exp,D,Ts);
step(4*sysDiscreto_Int)

%%
clear x2 up2 error; clc; close all;

% Se debe expandir el vector porque el control integral es un variabe de estado más que se añade al sistema
x2(:,1) = [0 0 0]'; % x2 = [IL Vc Error]
r = 3;
pert=0;

for i = 1:400
    if i >= 200
        pert = 0.5;    % Se crea una perturbacción en un momento arbitrario
    end
    up2(i) = -K2*x2(:,i) + pert;
    x2(:,i+1) = Af_Intg*x2(:,i) + Bf_Intg*up2(i) + [0 0 1]'*r;
end

figure (1)
stairs(up2)
title('Esfuerzo de control - Sistema discretizado')
ylabel('Voltaje entregado por el DAC [V]'); xlabel('Tiempo')
legend('Sistema con controlador')

figure (2)
stairs(x2(2,:))
title('Sistema discretizado controlado')
ylabel('Voltaje en el Capacitor [V]'); xlabel('Tiempo')
legend('Sistema con controlador')

%% Gráficas - Control Integral

data_Int = datos_Intg(:, 1);

% Crea una variable para almacenar los datos
DatosReales = data_Int(:);
t2 = linspace(0, 1, 1840);

figure(1)
plot(t2, DatosReales, 'r')
legend('Voltaje en el capacitor')
title('Voltaje en el capacitor frente a perturbarciones');
ylabel('Voltaje [V]');

figure(2)
plot(t2(450:550),DatosReales(560:660))
legend('Voltaje en el capacitor')
title('Zoom - Voltaje en el capacitor frente a perturbarciones');
ylabel('Voltaje [V]');

%% Función: Solucionar E.D.O's con ode45

% Define the system function
function dxdt = system(t, x)
    R1=1000;    % [Ohm]
    R2=220;     % [Ohm]
    R3=1000;    % [Ohm]
    L=1000e-6;  % [mHenrrios]
    C=470e-6;   % [Faradios]
    Vi=4;       % [Voltios] 
    dxdt = [x(1)*(-R2/L) - x(2)/L + Vi/L; x(1)/C - x(2)*((R1+R3)/(R1*R3*C)) + Vi/(R1*C)];
end