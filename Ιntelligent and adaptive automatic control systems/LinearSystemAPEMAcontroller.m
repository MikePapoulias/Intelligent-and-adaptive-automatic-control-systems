format long
clear workspace
close all
clc 


% time duration
t = 0:0.01:10;

% Parameters of real Model:
C = 1; 
G = 10;
M = 1;
kp = 1/M;

% Parameters of reference Model: 
am1 = 0;  am2 = 1;  am3 = -9;  am4 = -6;
bm1 = 0;  bm2 = 9;
km = 9;

% input r(t)
r = @(t) 1;

% Array F
F = -1;
% Array g
g = 1;
% Pole p0
p0 = 1;

% Array Gamma
g1 = 10000; g2 = 10000; g3 = 10000; g4 = 1;
Gamma = [g1 0 0 0 ; 0 g2 0 0 ; 0 0 g3 0 ; 0 0 0 g4];

% Initial Conditions:
w1_0 = 0;
w2_0 = 0;
f_0 = [0 ; 0 ; 0 ; 0];
th_0 = [0 ; 0 ; 0 ; 0];
x_0 = [0; 0];  % Real System
z_0 = [2; 1];  % Reference System 
y0 = [th_0' w1_0 w2_0 f_0' x_0' z_0'];


% Solve of differential equations:
options = odeset(Refine=100);
[t,diaf] = ode45(@(t,diaf) odefun(t,diaf,r,F,Gamma,C,G,M,g,p0,kp,km,am1,am2,am3,am4,bm1,bm2), t, y0, options);


% Controller's parameters:
w1 = diaf(:,1);
w2 = diaf(:,2);
f1 = diaf(:,3);
f2 = diaf(:,4);
f3 = diaf(:,5);
f4 = diaf(:,6);
th1 = diaf(:,7);
th2 = diaf(:,8);
th3 = diaf(:,9);
th4 = diaf(:,10);
y = diaf(:,11); % y = x1
x2 = diaf(:,12);
ym = diaf(:,13); % ym = z1
z2 = diaf(:,14);


% Graphs Results:
figure
hold on
plot(t, w1, 'blue');
plot(t, w2, 'red');
xlabel('t (sec)')
ylabel('ω_1,  ω_2');
title('ω table parameters')
legend('ω_1', 'ω_2');

figure
hold on
plot(t, f1, 'blue');
plot(t, f2, 'red');
plot(t, f3, 'green');
plot(t, f4, 'm');
xlabel('t (sec)')
ylabel('φ_1,  φ_2,  φ_3,  φ_4');
title('φ table parameters')
legend('φ_1', 'φ_2', 'φ_3', 'φ_4');

figure
hold on
plot(t, th1, 'blue');
plot(t, th2, 'red');
plot(t, th3, 'green');
plot(t, th4, 'm');
xlabel('t (sec)')
ylabel('θ_1,  θ_2,  θ_3,  θ_4');
title('θ table parameters')
legend('θ_1', 'θ_2', 'θ_3', 'θ_4');

figure
hold on
plot(t, y, 'blue');
plot(t, x2, 'red');
xlabel('t (sec)')
ylabel('x_1(t),   x_2(t)')
title('Real System State Variables x_1, x_2')
legend('x_1', 'x_2');

figure
hold on
plot(t, ym, 'blue');
plot(t, z2, 'red');
xlabel('t (sec)')
ylabel('z_1(t),   z_2(t)')
title('Reference System State Variables z_1, z_2')
legend('z_1', 'z_2');

figure
hold on
plot(t,y,'blue');
plot (t, ym,'red');
xlabel('t (sec)')
ylabel('y(t),   y_m(t)')
title('Outputs y, y_m')
legend('y', 'ym')

figure
plot(t, y-ym,'green');
xlabel('t (sec)')
ylabel('ε(t)')
title('Error ε = y-y_m')
legend('ε')






% Controller
function Controller = odefun(t,diaf,r,F,Gamma,C,G,M,g,p0,kp,km,am1,am2,am3,am4,bm1,bm2)
    % ======================================================
    % diaf_array:
    % (1) --> w1    (2) --> w2
    % (3) --> f1    (4) --> f2   (5) --> f3    (6) --> f4
    % (7) --> th1   (8) --> th2  (9) --> th3   (10) --> th4
    % (11) --> x1   (12) --> x2  (13) --> z1   (14) --> z2
    % ======================================================
    A = [0 1; 0 -(C/M)]; % real (non-linear) system for simulation 
    B = [0; 1/M];
    Am = [am1 am2; am3 am4];
    Bm = [bm1; bm2];

    w = [diaf(1)' ; diaf(2)' ; diaf(11) ; r(t)];
    f = [diaf(3) ; diaf(4) ; diaf(5) ; diaf(6)];
    th = [diaf(7) ; diaf(8) ; diaf(9) ; diaf(10)];
    x = [diaf(11); diaf(12)];
    z = [diaf(13); diaf(14)];
    epsilon = diaf(11) - diaf(13);
        
    dth = -Gamma*epsilon*f*sign(kp/km);
    u = th'*w + dth'*f;
    
    dx = A*x + B*(u - G*sin(diaf(11))); % real (non-linear) system for simulation 
    dz = Am*z + Bm*r(t);
    dw1 = F*diaf(1) + g*u;
    dw2 = F*diaf(2) + g*diaf(11);
    df = -p0*f + w;

    Controller = [dw1; dw2; df(1); df(2); df(3); df(4); dth(1); dth(2); dth(3); dth(4); dx(1); dx(2); dz(1); dz(2)];
end





