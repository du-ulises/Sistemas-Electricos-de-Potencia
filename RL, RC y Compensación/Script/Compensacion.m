%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              COMPENSACIÓN                               %
%                       Autor:Diego Villegas Govea                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

fprintf('Cálculo de un Sistema Eléctrico de AC\n\n');
disp('+---------------+');
disp('|               |');
disp('|               R1');
disp('|               |');
disp('V1              |');    
disp('|               L1');
disp('|               |');
disp('|               |');
disp('+---------------+');
fprintf('\n');

%PARÁMETROS DE ENTRADA
f=input('Ingrese el valor de la frecuencia (Hz): ');

V1=input('Ingrese la magnitud de V1 (Volts): ');
angV1=input('Ingrese el ángulo de V1: ');

R1=input('Ingrese el valor de R1 (Ohm): ');
Xvalor=input('Indique si es reactancia o valor inductivo de L1(0-Reactancia, 1-Valor inductivo): ');
if(Xvalor==0)
    Xreac=input('Ingrese el valor de la reactancia de L1: ');
elseif (Xvalor==1)
    Xvi=input('Ingrese el valor inductivo de L1(mH): ');    
end

% f=60;
% V1=220;
% angV1=0;
% R1=18.7;
% Xvalor=1;
% Xvi=50;

fprintf('\n');
disp('+----------------+');
disp('|                |');
disp('|                |');
disp('|                |');
disp('V1               Z1');    
disp('|                |');
disp('|                |');
disp('|                |');
disp('+----------------+');
fprintf('\n');

%IMPEDANCIA
if(Xvalor==0)
    Z1=complex(R1,Xreac);
elseif (Xvalor==1)
    Xvi=Xvi*10^-3;
    XL=2*pi*f*Xvi; %Valor inductivo
    Z1=complex(R1,XL);
end

[realV1, imagV1]=pol2cart((angV1*pi/180),V1);
disp(['V1 = ' num2str(complex(realV1,imagV1))]);
disp(['V1 = ' num2str(V1) ' |_ ' num2str(angV1) '° V']);
fprintf('\n');

disp(['Z1 = ' num2str(Z1)]);
[angZ1,magZ1] = cart2pol(real(Z1),imag(Z1));
angZ1=angZ1*180/pi;
disp(['Z1 = ' num2str(magZ1) ' |_ ' num2str(angZ1) '° ohms']);
fprintf('\n');

disp(['La impedancia equivalente es: ' num2str(Z1)]);
disp(['Zeq = ' num2str(magZ1) ' |_ ' num2str(angZ1) '° ohms']);
fprintf('\n');

%CORRIENTE
[realV1, imagV1]=pol2cart((angV1*pi/180),V1);

I1=complex(realV1, imagV1)/Z1;
disp(['I1 = ' num2str(I1)]);
[angI1,magI1] = cart2pol(real(I1),imag(I1));
angI1=angI1*180/pi;
disp(['I1 = ' num2str(magI1) ' |_ ' num2str(angI1) '° A']);
fprintf('\n');

%POTENCIA APARENTE
S=V1*magI1;
disp(['La potencia aparente es: ' num2str(S) ' VA']);

%POTENCIA ACTIVA
P=S*cos(angI1*pi/180);
disp(['La potencia activa es: ' num2str(P) ' W']);

%POTENCIA REACTIVA
Q=S*sin(angI1*pi/180);
if Q<0
    disp(['La potencia reactiva es: ' num2str(-Q) ' VAR (atraso)']);
else
    disp(['La potencia reactiva es: ' num2str(Q) ' VAR']);
end
fprintf('\n');

%FACTOR DE POTENCIA
FP=P/S;
disp(['El factor de potencia es: ' num2str(FP)]);

%SEÑALES SENOIDALES
w=2*pi*f;
t=0:0.0001:(1/f)*5;

Vm=V1*cos(w*t+angV1*pi/180);

Im=magI1*cos(w*t+angI1*pi/180);

%GRÁFICOS
%Gráfica de comparación de voltaje vs corriente
subplot(3,3,[1,6]), [hAx,hLine1,hLine2] = plotyy(t,Vm,t,Im,'plot'); 
title('SISTEMA ELÉCTRICO NO COMPENSADO'); 
xlabel('Tiempo'); 
ylabel(hAx(1),'Voltaje (V)'); 
ylabel(hAx(2),'Corriente (A)'),grid
V=['V = ' num2str(V1) ' V'];
I=['I = ' num2str(round(magI1,2)) ' A'];
set(hLine1, 'Color', [21/255, 67/255, 96/255]);
set(hLine1, 'LineWidth',2);
set(hLine2, 'Color', [129/255, 0, 0]);
set(hLine2, 'LineWidth',2);
set(hAx, {'ycolor'}, {[21/255, 67/255, 96/255]; [129/255, 0, 0]})
legend(V, I);
                    
%Gráfica de potencias
x=[0 100]; 
yS=[S S];
subplot(3,3,7), plot(x,yS,'LineWidth',2, 'Color', [255/255, 87/255, 51/255]); 
title('Potencia Aparente'); xlabel('Tiempo'); ylabel('VA'),grid
Etiq=['S = ' num2str(round(S,2)) ' VA'];
legend(Etiq);

yP=[P P];
subplot(3,3,8), plot(x,yP,'LineWidth',2, 'Color', [250/255, 210/255, 50/255]); 
title('Potencia Activa'); xlabel('Tiempo'); ylabel('W'),grid
Etiq=['P = ' num2str(round(P,2)) ' W'];
legend(Etiq);

yQ=[Q Q];
subplot(3,3,9), plot(x,yQ,'LineWidth',2, 'Color', [60/255, 150/255, 0]); 
title('Potencia Reactiva'); xlabel('Tiempo'); ylabel('VAR'),grid
Etiq=['Q = ' num2str(round(Q,2)) ' VAR'];
legend(Etiq);

%COMPENSAR EL SISTEMA
fprintf('\n');
opc=input('\n¿Deseas compensar el sistema? [si,no]: ','s');
if(opc=='si')


disp('   +------------+-------------+');
disp('   |            |             |');
disp(['   |          ' num2str(R1) 'ohm       ' num2str(R1) 'ohm']);
disp('   |            |             |');
disp([num2str(V1) '<' num2str(angV1) '°V' '         |             |']);    
if(Xvalor==0)
    C1=-Xreac;
    disp(['   |         '  num2str(round(C1,2)) 'i        ' num2str(Xreac) 'i']);
elseif (Xvalor==1)
    C1=1/(2*pi*f*XL*10^-6);
    disp(['   |         ' num2str(round(C1,2)) 'uF        ' num2str(Xvi*10^3) 'mH']);
end
disp('   |            |             |');
disp('   |            |             |');
disp('   +------------+-------------+');
fprintf('\n');
disp('+-------+-------+');
disp('|       |       |');
disp('|       |       |');
disp('V1      Z1      Z2');    
disp('|       |       |');
disp('|       |       |');
disp('+-------+-------+');
fprintf('\n');

[realV1, imagV1]=pol2cart((angV1*pi/180),V1);
disp(['V1 = ' num2str(complex(realV1,imagV1))]);
disp(['V1 = ' num2str(V1) ' |_ ' num2str(angV1) '° V']);
fprintf('\n');

Z2=Z1;
Z1=complex(R1,-XL);

disp(['Z1 = ' num2str(Z1)]);
[angZ1,magZ1] = cart2pol(real(Z1),imag(Z1));
angZ1=angZ1*180/pi;
disp(['Z1 = ' num2str(magZ1) ' |_ ' num2str(angZ1) '° ohms']);
fprintf('\n');

disp(['Z2 = ' num2str(Z2)]);
[angZ2,magZ2] = cart2pol(real(Z2),imag(Z2));
angZ2=angZ2*180/pi;
disp(['Z2 = ' num2str(magZ2) ' |_ ' num2str(angZ2) '° ohms']);
fprintf('\n');

disp('+-------------+');
disp('|             |');
disp('|             |');
disp('V1            Z3');    
disp('|             |');
disp('|             |');
disp('+-------------+');
fprintf('\n');

%Z3=Z1||Z2
Z3=(Z1*Z2)/(Z1+Z2);

disp(['Z3 = ' num2str(Z3)]);
[angZ3,magZ3] = cart2pol(real(Z3),imag(Z3));
angZ3=angZ3*180/pi;
disp(['Z3 = ' num2str(magZ3) ' |_ ' num2str(angZ3) '° ohms']);
fprintf('\n');

disp(['La impedancia equivalente es: ' num2str(Z3)]);
disp(['Zeq = ' num2str(magZ3) ' |_ ' num2str(angZ3) '° ohms']);
fprintf('\n');

%Corriente
[realV1, imagV1]=pol2cart((angV1*pi/180),V1);

I1=complex(realV1, imagV1)/Z3;
disp(['I1 = ' num2str(I1)]);
[angI1,magI1] = cart2pol(real(I1),imag(I1));
angI1=angI1*180/pi;
disp(['I1 = ' num2str(magI1) ' |_ ' num2str(angI1) '° A']);
fprintf('\n');

%POTENCIA APARENTE
S=V1*magI1;
disp(['La potencia aparente es: ' num2str(S) ' VA']);

%POTENCIA ACTIVA
P=S*cos(angI1*pi/180);
disp(['La potencia activa es: ' num2str(P) ' W']);

%POTENCIA REACTIVA
Q=S*sin(angI1*pi/180);
if Q<0
    disp(['La potencia reactiva es: ' num2str(-Q) ' VAR (atraso)']);
else
    disp(['La potencia reactiva es: ' num2str(Q) ' VAR']);
end
fprintf('\n');

%FACTOR DE POTENCIA
FP=P/S;
disp(['El factor de potencia es: ' num2str(FP)]);

%SEÑALES SENOIDALES
w=2*pi*f;
t=0:0.0001:(1/f)*5;

Vm=V1*cos(w*t+angV1*pi/180);

Im=magI1*cos(w*t+angI1*pi/180);

%GRÁFICOS
figure(2)
%Gráfica de comparación de voltaje vs corriente
subplot(3,3,[1,6]), [hAx,hLine1,hLine2] = plotyy(t,Vm,t,Im,'plot'); 
title('SISTEMA ELÉCTRICO COMPENSADO'); 
xlabel('Tiempo'); 
ylabel(hAx(1),'Voltaje (V)'); 
ylabel(hAx(2),'Corriente (A)'),grid
V=['V = ' num2str(V1) ' V'];
I=['I = ' num2str(round(magI1,2)) ' A'];
set(hLine1, 'Color', [21/255, 67/255, 96/255]);
set(hLine1, 'LineWidth',2);
set(hLine2, 'Color', [129/255, 0, 0]);
set(hLine2, 'LineWidth',2);
set(hAx, {'ycolor'}, {[21/255, 67/255, 96/255]; [129/255, 0, 0]})
legend(V, I);
                    
%Gráfica de potencias
x=[0 100]; 
yS=[S S];
subplot(3,3,7), plot(x,yS,'LineWidth',2, 'Color', [255/255, 87/255, 51/255]); 
title('Potencia Aparente'); xlabel('Tiempo'); ylabel('VA'),grid
Etiq=['S = ' num2str(round(S,2)) ' VA'];
legend(Etiq);

yP=[P P];
subplot(3,3,8), plot(x,yP,'LineWidth',2, 'Color', [250/255, 210/255, 50/255]); 
title('Potencia Activa'); xlabel('Tiempo'); ylabel('W'),grid
Etiq=['P = ' num2str(round(P,2)) ' W'];
legend(Etiq);

yQ=[Q Q];
subplot(3,3,9), plot(x,yQ,'LineWidth',2, 'Color', [60/255, 150/255, 0]); 
title('Potencia Reactiva'); xlabel('Tiempo'); ylabel('VAR'),grid
Etiq=['Q = ' num2str(round(Q,2)) ' VAR'];
legend(Etiq);
end