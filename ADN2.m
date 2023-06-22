%% © Deepa, PhD Scholar - IIT Roorkee
%% ___AD Noise in QSS (security check) by Grover's Algorithm __________
clc;
clear all;

%% _________________Pauli Matrices___________________
I=[1 0;0 1];
X=[0 1;1 0];
Y=[0 -1i;1i 0];
Z=[1 0;0 -1];
H=(1/sqrt(2))*[1,1;1,-1];
CX=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0];
CZ=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,-1];

% ____________________Bell Basis____________________
a= [1;0];       % |0>
b= [0;1];       % |1>
aa=kron(a,a);   % |00>
ab=kron(a,b);   % |01>
ba=kron(b,a);   % |10>
bb=kron(b,b);   % |11>



aaa = kron(aa,a);      % |000>
aab = kron(aa,b);      % |001>
aba = kron(ab,a);      % |010>
abb = kron(ab,b);      % |011>
baa = kron(ba,a);      % |100>
bab = kron(ba,b);      % |101>
bba = kron(bb,a);      % |110>
bbb = kron(bb,b);      % |111>

% n=0;
% syms n;
FF1 =[];
FF2 =[];
FF3 =[];
FF4 =[];
FF5 =[];
FF6 =[];
FF7 =[];
FF8 =[];

GG = [];
mm = [];

for n =0:0.05:1

A = -1+3*sqrt(1-n) +3*(1-n) + (1-n)^(3/2);

%B = 6*A-4-12*sqrt(1-n)-12*(1-n)-4*(1-n)^(3/2);

C1 = (A+4)*aaa ;
C2 = (A-4*sqrt(1-n))*(-aab + aba + baa);
C3 = (-A+4*(1-n))*(abb+bab); 
C4 = (A-4*(1-n))*(bba);
C5 = (-A+4*(1-n)^(3/2))*(bbb);

C = C1+C2+C3+C4+C5;
CC = C*C';
D = (-3*aaa-aab+aba-abb+baa-bab+bba-bbb);
DD = (n^3)*(D*D');

rho = (1/128)*(CC+DD);

F1 = rho(1,1); 
F2 = rho(2,2);
F3 = rho(3,3);
F4 = rho(4,4);
F5 = rho(5,5);
F6 = rho(6,6);
F7 = rho(7,7);
F8 = rho(8,8);

FF1=[FF1 F1];
FF2=[FF2 F2];
FF3=[FF3 F3];
FF4=[FF4 F4];
FF5=[FF5 F5];
FF6=[FF6 F6];
FF7=[FF7 F7];
FF8=[FF8 F8];


end

p=20;
NoiseP=[0:1/p:1];

%% Plotting
figure
xlim([0 1])
ylim([0 1])
plot(NoiseP,FF1,'b--^')

hold on
plot(NoiseP,FF2,'k-.*')

hold on
plot(NoiseP,FF3,'m--.')

hold on
plot(NoiseP,FF4,'g--d')

hold on
plot(NoiseP,FF5,'c--*')

hold on
plot(NoiseP,FF6,'b-.')

hold on
plot(NoiseP,FF7,'r--p')

hold on
plot(NoiseP,FF8,'-o')

hold on
legend('|000>','|001>','|010>','|100>','|011>','|101>','|110>','|111>')
ylabel('Fidelity (Probability)') 
xlabel('Noise parameter (0\leq\mu\leq1)') 