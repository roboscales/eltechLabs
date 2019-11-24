%общие данные
Em1 = 141;
Em2 = 179;
psi = deg2rad(-40); % тут поменяешь на свое значение угла(буковка пси)
e1 = Em1/sqrt(2);
e2 = (Em2/sqrt(2))*exp(1i*psi);
disp(e1);
disp(e2);
z1 = 1i*10; %здесь свои значения понапишешь
z2 = 30; %здесь
z3 = 1i*40; %здесь
z4 = 50;%здесь
z5 = -1i*60;%здесь

z = [1i*10; 30; 1i*40; 50; -1i*60];
ZI = [1 1 -1 0 0; 0 0 1 -1 -1; z1 -z2 0 0 0; 0 z2 z3 z4 0; 0 0 0 -z4 z5];%здесь поменяешь на свою матрицу
sum = e1 + e2;
E = [0; 0; -e1; sum; 0];

I = linsolve(ZI, E);

display(sum,'Cumma dvux tokov');

display('По законам Ома-Кирхгофа');
display(I);

disp('Pokazatelnaya forma I');
exForm(I);

U = z .* I;

display('Compute a voltage');
display(U);
exForm(U);

%Метод контурных токов
R11 = z1+z2;%here
R22 = z2+z3+z4;%here
R33 = z4+z5;%here
R12 = z2;%here
R21 = z2;%here
R23 = z4;%here
R32 = z4;%here
R13 = 0; %here
R31 = 0;%here
EK1 = e1;%here
EK2 = -e1-e2;%here
EK3 = 0;%here
ZIK = [R11 -R12 -R13; -R21 R22 -R23; -R31 -R32 R33];
EK = [EK1; EK2; EK3];

IK = linsolve(ZIK, EK);
display(IK);
display(EK2);

Ik = zeros(length(I));

Ik(1) = -IK(1);
Ik(2) = -IK(2)+IK(1);
Ik(3) = -IK(2);
Ik(4) = IK(3)-IK(2);
Ik(5) = -IK(3);

display('Метод контурных токов');
display(Ik);

G11=((1/z1)+(1/z2)+(1/z3));
G12 = 1/z3;
G21 = 1/z3;
G22 = ((1/z3)+(1/z4)+(1/z5));
I11 = (e1/z2)-(e2/z3);
I22 = (e2/z3);

%матрица коэффициентов для потенциалов
P = [G11 -G12; -G21 G22];
display(P);

PE = [exForm(G11) exForm(-G12); exForm(-G21) exForm(G22)];
display(PE);
PI = [I11; I22];

display(PI);
exForm(PI);

del = det(P);
display(del);
exForm(del);

del11 = ((-1)^(1+1)) * G22;
del22 = ((-1)^(2+2)) * G11;
del12 = ((-1)^(1+2)) * (-G12);
del21 = ((-1)^(1+2)) * (-G21);

display(del11);
display(del22);
display(del12);
display(del21);

hi1 = ((del11/del)*I11) + ((del12/del)*I22);
display(hi1);
hi2 = ((del21/del)*I11) + ((del22/del)*I22);
display(hi2);

hi3 = 0;

U1 = hi3 - hi1;
U2 = e1- (hi1 - hi3);
U3 = e2-(hi2-hi1);
U4 = hi2-hi3;
U5 = hi2 - hi3;

IP = zeros(length(I));

IP(1) = U1/z1;
IP(2) = U2/z2;
IP(3) = U3/z3;
IP(4) = U4/z4;
IP(5) = U5/z5;

display('Метод узловых потенциалов');
display(IP);

U =[I(1)*z1;I(2)*z2;(I(3)*z3);(I(4)*z4);I(5)*z5];

%%векторная диаграмма
compass(I);
hold on;
compass(U, '--');

%%checking balance of power
sum1 = (I(1)^2)*z1+(I(2)^2)*z2+(I(3)^2)*z3+(I(4)^2)*z4+(I(5)^2)*z5;
SUMI = I.^2;
SUMI = SUMI.*z;
display(SUMI);
display(sum1);
sum2 = e1*I(2)+e2*I(3);
display(sum2);
Pis = real(sum2);
Qis = imag(sum2);
Pprost = real(sum1);
Qprost = imag(sum1);

dQ = (abs(Qis - Qprost)/Qis)*100;
dP = (abs(Pis - Pprost)/Pis)*100;

display('Оценки погрешностей:');
display(dQ);
display(dP);

%what does the AMP and VOLT say
%i'm so bored! i want to invent a new function!
A = comtopres(I(2))*sqrt(2);
display(A);
V = comtopres(I(4)*z4)*sqrt(2);
display(V);

c = cos(atan(imag(I(2)*z2)/real(I(2)*z2)-atan(imag(I(3))/real(I(3)))));
hc = rad2deg(atan(imag(I(2)*z2)/real(I(2)*z2)-atan(imag(I(3))/real(I(3)))));
display(hc);

W = comtopres(I(2)*z2)*comtopres(I(3))*2*c;
display(W);



disp('AMP I');
IA = abs(I).*sqrt(2);
display(IA);

disp('AMP U');
UA = abs(U).*sqrt(2);
display(UA);

%%method of equal generator
% z45 = (z4*z5)/(z4+z5);
% 
% zq1 = z2+z3+z45;
% I21 = e1/zq1;
% 
% % zq2 = z3+(z2*(z4*z5)/(z4+z5))/(z2+(z4*z5)/(z4+z5));
% %zq2 = (z3+z2)*z45/(z3+z2+z45);
% zq2 = z2+z3+z45;
% I22 = e2/zq2;
% 
% IN2 = I21+I22;
% 
% zeq = z2*(z3+z45)/(z2+z3+z45);
% %zeq = z2+z3+(z4*z5)/(z4+z5);
% 
% IE1 = (IN2*z2)/(zeq+z1);
% 
% disp(IE1);

%first ist
ze1 = z2+z3+(z4*z5/(z4+z5));
I21 = e1/ze1;

exForm(z4*z5);
exForm(z4+z5);

display(ze1);
exForm(ze1);

display(I21);
exForm(I21);

display(I22);
exForm(I22);

ze2 = z2+z3+(z4*z5/(z4+z5));
I22 = e2/ze2;

IC2 = I21+I22;
zq = (z2*(z3+(z4*z5/(z4+z5)))/(z2+z3+(z4*z5/(z4+z5))));

display(IC2);
exForm(IC2);

EDS = e1-z2*IC2;

display(EDS);
exForm(EDS);

exForm(z2*(z3+(z4*z5/(z4+z5))));
exForm((z2+z3+(z4*z5/(z4+z5))));

display(zq);
exForm(zq);

IE1 = (e1-z2*IC2)/(zq+z1);

exForm(zq+z1);

disp(IE1)