%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Señales Periódicas                            % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

fprintf('Señal periódica de f(x)=x^2, [-pi, pi]\n');
N1=input('Número de armónicas: ');

fprintf('\nSeñal periódica de f(x)=Vm*Sen(wt), [0, 3*pi]\n');
N=input('Número de armónicas: ');
vm=input('Ingrese la amplitud del rectificador trifásico: ');

Vo=(2*(pi*pi))/3;
an=((4*((-1)^N1))/(N1*N1)); % Vector de componentes coseno
wt=-pi:0.0001:pi; % Vector wt
vo=Vo+an*cos(N1*wt); % Vector vo
subplot(2,1,1);
plot(wt,vo,'r','LineWidth',1.5),grid
title('Señal periódica de f(x)=x^2, [-pi, pi]'); 
xlabel('Tiempo (rad/seg)'); 
ylabel('Amplitud (V)'); 
legend('f(x)=x^2');
set(gca,'XTick',-2*pi:pi/4:2*pi) ;
set(gca,'XTickLabel',{'-2*pi','-7*pi/4','-3*pi/2','-5*pi/4','-pi','-3*pi/4','-pi/2','-pi/4','0','pi/4','pi/2','3*pi/4','pi','5*pi/4','3*pi/2','7*pi/4','2*pi'});

Vo=3*sqrt(3)*vm/(2*pi);
n=1:N;
an=-3*vm/(2*pi)*((cos((1+3*n)*5*pi/6)-cos((1+3*n)*pi/6))./(1+3*n)+(cos((1-3*n)*5*pi/6)-cos((1-3*n)*pi/6))./(1-3*n));
bn=3*vm/(2*pi)*((sin((1-3*n)*5*pi/6)-sin((1-3*n)*pi/6))./(1-3*n)-(sin((1+3*n)*5*pi/6)-sin((1+3*n)*pi/6))./(1+3*n));
wt=0:0.001:3*pi;
vo=Vo+an*cos(3*n'*wt)+bn*sin(3*n'*wt); % Vector vo
subplot(2,1,2);
plot(wt,vo,'g','LineWidth',1.5),grid
title('Señal periódica de f(x)=Vm*sen(wt), [0, 3*pi]'); 
xlabel('Tiempo (rad/seg)'); 
ylabel('Amplitud (V)'); 
etiq=['f(x)=' num2str(vm) '*sen(wt)'];
legend(etiq);
set(gca,'XTick',0:pi/6:3*pi) ;
set(gca,'XTickLabel',{'0','pi/6','pi/3','pi/2','2*pi/3','5*pi/6','pi','7*pi/6','4*pi/3','3*pi/2','5*pi/3','11*pi/6','2*pi','13*pi/6','7*pi/3','5*pi/2','8*pi/3','17*pi/6','3*pi'});
