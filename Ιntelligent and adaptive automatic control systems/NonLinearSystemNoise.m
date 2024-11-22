format long
clear workspace
close all
clc 


% time duration
t = 0:0.01:120;

% Parameters of real Model:
C = 1; 
G = 10;
M = 1;

% Parameters of reference Model: 
aref1 = 0;  aref2 = 1;  aref3 = -1;  aref4 = -1.4;
bref1 = 0;  bref2 = 1;

% input r(t)
r = @(t) 1;

% Noise d(t)
d = @(t) noise(t,20);

% Array Q
Q = eye(2);
% Solve Lyapunov Equarion (Aref'*P + P*Aref = -Q):
Aref = [aref1 aref2 ; aref3 aref4];
P = lyap(Aref',Q);

% Profits g1, g2, g3:
g1 = 1000;
g2 = 10;
g3 = 1000;

% Initial Conditions:
x_0 = [0 ; 0];   % Real System
xref_0 = [2 ; 1]; % Reference System 
k_0 = [0 ; 0];
L_0 = 0;
N_0 = 0;
y0 = [x_0' xref_0' k_0' L_0 N_0];



% Solve of differential equations:
options = odeset(Refine=100);
[t,diaf] = ode45(@(t,diaf) odefun(t,diaf,r,C,G,M,aref1,aref2,aref3,aref4,bref1,bref2,g1,g2,g3,P,d), t, y0, options);




% Controller's parameters:
x1 = diaf(:,1);
x2 = diaf(:,2);
xref1 = diaf(:,3);
xref2 = diaf(:,4);
k1 = diaf(:,5);
k2 = diaf(:,6);
L = diaf(:,7);
N = diaf(:,8);


% Graphs Results:
figure
hold on
plot(t,k1,'blue');
plot(t,k2,'red');
xlabel('t (sec)')
ylabel('k(t) = [k1(t) k2(t)]')
title('k table parameters')
legend('k1','k2')

figure
plot(t,L,'blue');
xlabel('t (sec)')
ylabel('L(t)')
title('L parameter')
legend('L')

figure
plot(t,N,'blue');
xlabel('t (sec)')
ylabel('N(t)')
title('N parameter')
legend('N')

figure
hold on
plot(t,x1,'blue');
plot(t,xref1,'red');
plot(t, x1-xref1,'green');
xlabel('t (sec)')
ylabel('e_1(t)')
title('Error e_1 = x1-x_r_e_f1')
legend('x_1','x_r_e_f_1','e_1')

figure
hold on
plot(t,x2,'blue');
plot(t,xref2,'red');
plot(t, x2-xref2,'green');
xlabel('t (sec)')
ylabel('e_2(t)')
title('Error e_2 = x2-x_r_e_f2')
legend('x_2','x_r_e_f_2','e_2')





% Controller
function Controller = odefun(t,diaf,r,C,G,M,aref1,aref2,aref3,aref4,bref1,bref2,g1,g2,g3,P,d)
    % ======================================================
    % diaf_array:
    % (1) --> x1      (2) --> x2
    % (3) --> xref1   (4) --> xref2  
    % (5) --> k1      (6) --> k2       
    % (7) --> L 
    % (8) --> N 
    % ======================================================
    A = [0 1; 0 -(C/M)];
    B = [0; 1/M];
    Aref = [aref1 aref2; aref3 aref4];
    Bref = [bref1; bref2];

    x = [diaf(1); diaf(2)];
    xref = [diaf(3); diaf(4)];
    e = x - xref;
    f = -G*sin(diaf(1));
    Fi = sin(diaf(1));
    k = [diaf(5); diaf(6)];
    L = diaf(7);
    N = diaf(8);
        
    u = -k'*x -L*r(t) - N*Fi;
    
    dx = A*x + B*(u + d(t) + f);
    dxref = Aref*xref + Bref*r(t);
    dk = -g1*Bref'*P*e*x'*signNumber(L);
    dL = -g2*Bref'*P*e*r(t)'*signNumber(L);
    dN = -g3*Bref'*P*e*Fi'*signNumber(L);

    Controller = [dx(1); dx(2); dxref(1); dxref(2); dk(1); dk(2); dL; dN];
end



function result = signNumber(arg)
    if (arg<0.0000001 && arg>-0.0000001)
        if(arg>0)
            s = 1;
        else
            s = -1;
        end
    else
        s = sign(arg);
    end

    result = s;
end



function output_signal = noise(t, A)
    % Initialize of output:
    output_signal = zeros(size(t));
    % Noise
    for i = 1:length(t)
        if ((t(i) >= 50 && t(i) < 55) || (t(i) >= 60 && t(i) < 65) || (t(i) >= 80 && t(i) < 85))
            output_signal(i) = A;
        end
    end
end


