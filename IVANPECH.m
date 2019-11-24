%общие данные
Em1 = 141;
Em2 = 179;
psi = deg2rad(80); % тут поменяешь на свое значение угла(буковка пси)
e1 = Em1/sqrt(2);
e2 = (Em2/sqrt(2))*exp(1i*psi);
disp(e1);
disp(e2);
z1 = -1i*10; %здесь свои значения понапишешь
z2 = 1i*30; %здесь
z3 = 1i*40; %здесь
z4 = 50;%здесь
z5 = 60;%здесь


ZI = [1 1 -1 0 0; 0 0 1 -1 -1; z1 -z2 0 0 0; 0 z2 z3 z4 0; 0 0 0 -z4 z5];%здесь поменяешь на свою матрицу
sum = e1 - e2;
E = [0; 0; sum; e2; 0];

I = linsolve(ZI, E);

disp(sum);
disp('По законам Ома-Кирхгофа');
display(I);

UOK = zeros(5,1);

UOK(1) = I(1)*z1;
UOK(2) = I(2)*z2;
UOK(3) = I(3)*z3;
UOK(4) = I(4)*z4;
UOK(5) = I(5)*z5;

display(UOK);

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
EK1 = e2-e1;%here
EK2 = -e2;%here
EK3 = 0;%here

ZIK = [R11 -R12 -R13; -R21 R22 -R23; -R31 -R32 R33];
display(ZIK);
EK = [EK1; EK2; EK3];
display(EK);

IK = linsolve(ZIK, EK);
display(IK);

Ik = zeros(length(I));

Ik(1) = -IK(1);
Ik(2) = IK(1)-IK(2);
Ik(3) = -IK(2);
Ik(4) = -IK(2)+IK(3);
Ik(5) = -IK(3);

display('Метод контурных токов');
display(Ik);

%%метод узловых потенциалов
G11=((1/z1)+(1/z2)+(1/z3));
display(G11);
exForm(G11);
G12 = 1/z3;
display(G12);
exForm(G12);
G21 = 1/z3;
G22 = ((1/z3)+(1/z4)+(1/z5));
display(G22);
exForm(G22);
I11 = (e1/z1)+(e2/z2);
display(I11);
exForm(I11);
I22 = 0;
display(I22);
%exForm(I22);
%матрица коэффициентов для потенциалов
P = [G11 -G12; -G21 G22];
PI = [I11; I22];

del = det(P);
display(del);
exForm(del);

del11 = ((-1)^(1+1)) * G22;
del22 = ((-1)^(2+2)) * G11;
del12 = ((-1)^(1+2)) * (-G12);
del21 = ((-1)^(1+2)) * (-G21);

display(del11);
exForm(del11);
display(del22);
exForm(del22);
display(del12);
exForm(del12);
display(del21);
exForm(del21);

hi1 = ((del11/del)*I11) + ((del12/del)*I22);
display(hi1);
exForm(hi1);
hi2 = ((del21/del)*I11) + ((del22/del)*I22);
display(hi2);
exForm(hi2);

hi3 = 0;

UM = zeros(5,1);

UM(1)= e1 - (hi1 - hi3);
UM(2) = e2 - (hi1 - hi3);
UM(3) = (hi1-hi2);
UM(4) = (hi2-hi3);
UM(5) = hi2 - hi3;

display(UM);

IP = zeros(length(I));

IP(1) = UM(1)/z1;
IP(2) = UM(2)/z2;
IP(3) = UM(3)/z3;
IP(4) = UM(4)/z4;
IP(5) = UM(5)/z5;

display('Метод узловых потенциалов');
display(IP);

U =[I(1)*z1;I(2)*z2;(I(3)*z3);(I(4)*z4);I(5)*z5];

%%векторная диаграмма
compass(I);
hold on;
compass(U, '--');

%%checking balance of power
sum1 = (I(1)^2)*z1+(I(2)^2)*z2+(I(3)^2)*z3+(I(4)^2)*z4+(I(5)^2)*z5;
display(sum1);
exForm(sum1);
sum2 = e1*I(1)+e2*I(2);
display(sum2);
exForm(sum2);
Pis = real(sum2);
Qis = imag(sum2);
Pprost = real(sum1);
Qprost = imag(sum1);

display(Pis);
exForm(Pis);
display(Qis);
exForm(Qis);
display('Slagaymie')
display((I(1)^2)*z1);
exForm((I(1)^2)*z1);
display((I(2)^2)*z2);
exForm((I(2)^2)*z2);
display((I(3)^2)*z3);
exForm((I(3)^2)*z3);
display((I(4)^2)*z4);
exForm((I(4)^2)*z4);
display((I(5)^2)*z5);
exForm((I(5)^2)*z5);

display(e1*I(3));
exForm(e1*I(3));
display(e2*I(4));
exForm(e2*I(4));

display(Pprost);
exForm(Pprost);
display(Qprost);
exForm(Qprost);

dQ = 100*(abs(Qis - Qprost)/Qis);
dP = 100*(abs(Pis - Pprost)/Pis);

display(abs(Qis - Qprost));
exForm(abs(Qis - Qprost));
display(abs(Pis - Pprost));
exForm(abs(Pis - Pprost));

display('Оценки погрешностей:');
display(dQ);
exForm(dQ);
display(dP);
exForm(dP);

%%what does the AMP and VOLT say
%%i'm so bored! i want to invent a new function!
A = comtopres(I(2))*sqrt(2);
display(A);
V = comtopres(I(4)*z4)*sqrt(2);
display(V);

c = cos(atan(imag(I(2)*z2)/real(I(2)*z2)-atan(imag(I(3))/real(I(3)))));
pc = rad2deg(atan(imag(I(2)*z2)/real(I(2)*z2)-atan(imag(I(3))/real(I(3)))));
display(pc);

W = comtopres(I(2)*z2)*comtopres(I(3))*c*2;
display(comtopres(I(2)*z2));
display(comtopres(I(3)));
display(W);

%%method of equal generator
% z145 = ((1/z1)+(1/z4)+(1/z5))^(-1);
% ze1 = z3+z145;
% I11 = e1/ze1;
% ze2 = z4+(z3*(z1+z5)/(z1+z3+z5));
% I42 = e2/ze2;
% %U12 = I42 * (z3*(z1+z5))/(z1+z3+z5);
% U12 = e2 - I42*z4;
% I12 = U12/(z1+z3);
% IC = I11+I12;
% Eeg = IC*z1;
% ;
% Zeg = (z3*(z1+z45))/(z1+z3+z45);
% IE2 = Eeg/(Zeg+z2);
% display(IE2);
z45 = z4*z5/(z4+z5);
z15 = z1*z5/(z1+z5);

ze1 = z1+(z2*(z3+z4))/(z2+z3+z4);
display(ze1);
I11 = e1/ze1;
display(I11);
I41 = (e1-I11*z1)/(z3+z4);

ze2 = z2 + (z1*(z3+z4)/(z1+z3+z4));
display(ze2);
I22 = e2/ze2;
%display(I42);
U = I42*(z5*(z3+z1)/(z1+z3+z5));
display(U);
I42 = (e2-I22*z2)/(z4+z3);
display(I42);


IE = I42 + I41;
display(IE);


ze = (z4*(z3+((z1*z2)/(z1+z2))))/(z4+z3+((z1*z2)/(z1+z2)));
IE2 = z4*IE/(ze+z5);

display(z1*IE);

display(ze);
display(z1*(z3+((z4*z5)/(z4+z5))));

display(IE2);
display(ze+z2);
Im = zeros(5,1);
Um = zeros(5,1);

for c = 1:5
    Im(c) = comtopres(I(c))*sqrt(2);
end

display(Im);

for c = 1:5
    Um(c) = comtopres(UOK(c))*sqrt(2);
end

display(Um);

Ia = zeros(5,1);
Ua = zeros(5,1);

for c = 1:5
    Ia(c) = comtopres(I(c));
end

display(Ia);

for c = 1:5
    Ua(c) = comtopres(UOK(c));
end

display(Ua);

Ih = zeros(5,1);
Uh = zeros(5,1);

for c = 1:5
    Ih(c) = rad2deg(angle(I(c)));
end

display(Ih);

for c = 1:5
    Uh(c) = rad2deg(angle(UOK(c)));
end

exForm(I);


display(Uh);

exForm(UOK);

ze1 = z1+(z2*(z3+z4))/(z2+z3+z4);
exForm(z2*(z3+z4));
exForm(z2+z3+z4);
display(ze1);
exForm(ze1);

I11 = e1/ze1;
exForm(e1);

display(I11);
exForm(I11);

I41 = (e1-I11*z1)/(z3+z4);
exForm(e1-I11*z1);
exForm(z3+z4);
display(I41);
exForm(I41);

display('Second ist of EDS');
ze2 = z2 + (z1*(z3+z4)/(z1+z3+z4));
exForm(z1*(z3+z4));
exForm(z1+z3+z4);
display(ze2);
exForm(ze2);
disp('Tok22');
I22 = e2/ze2;
display(I22);
exForm(I22);

%display(I42);
%U = I42*(z5*(z3+z1)/(z1+z3+z5));
%display(U);
disp('Iskomiy Tok');
I42 = (e2-I22*z2)/(z4+z3);
exForm(e2-I22*z2);
exForm(z4+z3);
display(I42);
exForm(I42);


IE = I42 + I41;
display(IE);
exForm(IE);


ze = (z4*(z3+((z1*z2)/(z1+z2))))/(z4+z3+((z1*z2)/(z1+z2)));
IE2 = z4*IE/(ze+z5);

display(z4*IE, 'EDS');
exForm(z4*IE);

disp('PORA');
exForm(z4*(z3+((z1*z2)/(z1+z2))));
exForm(z4+z3+((z1*z2)/(z1+z2)));

display(ze);
exForm(ze);
display(z1*(z3+((z4*z5)/(z4+z5))));

display(IE2);
display(ze+z5);
exForm(ze+z5);