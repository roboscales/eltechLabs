%общие данные
Em1 = 141;
Em2 = 179;
psi = deg2rad(-50); % тут поменяешь на свое значение угла(буковка пси)
e1 = Em1/sqrt(2);
e2 = (Em2/sqrt(2))*exp(1i*psi);
disp(e1);
disp(e2);
z1 = 1i*10; %здесь свои значения понапишешь
z2 = 30; %здесь
z3 = -1i*40; %здесь
z4 = 50;%здесь
z5 = -1i*60;%здесь
z = [1i*10 30 -1i*40 50 1i*60];

ZI = [1 -1 -1 0 0; 0 0 1 1 -1; -z1 -z2 0 0 0; 0 z2 -z3 z4 0; 0 0 0 -z4 -z5];%здесь поменяешь на свою матрицу
sum = -e1 + e2;
E = [0; 0; 0; sum; -e2];

I = linsolve(ZI, E);

disp(sum);
disp('По законам Ома-Кирхгофа');
display(I);

UOK = zeros(5,1);

for c=1:5
    UOK(c) = I(c)*z(c);
end    

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
EK1 = 0;%here
EK2 = -e1+e2;%here
EK3 = -e2;%here

ZIK = [R11 -R12 -R13; -R21 R22 -R23; -R31 -R32 R33];
display(ZIK);
EK = [EK1; EK2; EK3];
display(EK);

IK = linsolve(ZIK, EK);
display(IK);

Ik = zeros(length(I));

Ik(1) = -IK(1);
Ik(2) = -IK(1)+IK(2);
Ik(3) = -IK(2);
Ik(4) = IK(2)-IK(3);
Ik(5) = -IK(3);

display('Метод контурных токов');
display(Ik);

%%метод узловых потенциалов
G11=((1/z1)+(1/z2)+(1/z3));
display(G11);
G12 = 1/z3;
display(G12);
G21 = 1/z3;
G22 = ((1/z3)+(1/z4)+(1/z5));
display(G22);
I11 = -(e1/z3);
display(I11);
I22 = (e1/z3)+(e2/z4);
display(I22);
%матрица коэффициентов для потенциалов
P = [G11 -G12; -G21 G22];
PI = [I11; I22];

del = det(P);
display(del);

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

UM = zeros(5,1);

UM(1)= hi3 - hi1;
UM(2) = hi1 - hi3;
UM(3) = e1 - (hi2-hi1);
UM(4) = e2 - (hi2-hi3);
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
sum2 = e1*I(3)+e2*I(4);
display(sum2);
Pis = real(sum2);
Qis = imag(sum2);
Pprost = real(sum1);
Qprost = imag(sum1);

display(Pis);
display(Qis);
display((I(1)^2)*z1);
display((I(2)^2)*z2);
display((I(3)^2)*z3);
display((I(4)^2)*z4);
display((I(5)^2)*z5);

display(e1*I(3));
display(e2*I(4));

display(Pprost);
display(Qprost);

dQ = 100*(abs(Qis - Qprost)/Qis);
dP = 100*(abs(Pis - Pprost)/Pis);

display(abs(Qis - Qprost));
display(abs(Pis - Pprost));

display('Оценки погрешностей:');
display(dQ);
display(dP);

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

%meqg
zq2 = z4+z5;
I42 = e2/zq2;
% eeg = 

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
% z45 = z4*z5/(z4+z5);
% z15 = z1*z5/(z1+z5);
% 
% ze1 = z3+z1+z45;
% display(ze1);
% I11 = e1/ze1;
% display(I11);
% 
% ze2 = z4 + (z5*(z3+z1)/(z1+z3+z5));
% display(ze2);
% I42 = e2/ze2;
% display(I42);
% U = I42*(z5*(z3+z1)/(z1+z3+z5));
% display(U);
% I12 = U/(z1+z3);
% display(I12);
% 
% 
% IE = I12 - I11;
% display(IE);
% 
% 
% ze = (z1*(z3+((z4*z5)/(z4+z5))))/(z1+z3+((z4*z5)/(z4+z5)));
% IE2 = z1*IE/(ze+z2);
% 
% display(IE2);
% 
% display(z1*IE);
% 
% display(ze);
% display(z1*(z3+((z4*z5)/(z4+z5))));
% 
% display(IE2);
% display(ze+z2);
% Im = zeros(5,1);
% Um = zeros(5,1);



% for c = 1:5
%     Im(c) = comtopres(I(c))*sqrt(2);
% end
% 
% display(Im);
% 
% for c = 1:5
%     Um(c) = comtopres(UOK(c))*sqrt(2);
% end
% 
% display(Um);
% 
