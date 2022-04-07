% Hypersonics Aerothermodynamics
% Developed by Alisson Vinicius Brito Lopes - 18/11/2014
% Computing shock wave angle under chemical equilibrium conditions

M1c = 25; % Select a Mach number from 4 up to 25
theta = 20; % Inlet ramp angle, typical of waverider application

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1; i=1; ii=1; jj=1; 
tol=1e-6; % tol
theta=theta*(pi()/180);
gammac = 1.4; 
b=0.001; % error tol for the iterative while loop
T0 = 151.78; % Reference Temperature paper Tannehil
R=287.17; 
T1c = 226.650; % Temperature free stream at 30 Km
cp1c = 1006; 
a1c = sqrt(gammac*R*T1c);
rho1c = 1.8012e-2; % Density free stream 
p1c = 1171.8; % Pressure free stream  N/m²
h1c = cp1c*T1c; % Enthalpy free stream  J/kg
u1c = M1c*a1c;

%%%%%%%% Caloric perfect gas calculation %%%%%%%%%%%%%%%%%%%%
% Iterative process to compute Beta angle

beta = 18*(pi/180); % Initial guessing
a=tan(theta);

while b <= a
b(i)=2*cot(beta)*(((M1c^2*(sin(beta))^2)-1)/(M1c^2*(gammac+cos(2*beta))+2)); % Equation to be equaled using the 'new' Beta
beta = beta + 0.0001; % Angle increment
i=i+1;
end

betac=beta;
Abeta=betac*(180/pi()); 
Mn1=M1c*sin(betac); %% Mach Normal to shockwave - Eq. 4.7 Anderson (compressible Flow)
rho2c=rho1c*(((gammac+1)*(Mn1^2))/(((gammac-1)*(Mn1^2)) +2)); %% Eq. 4.8
p2c=p1c*(1+(2*gammac/(gammac+1))*((Mn1^2)-1));  %% Eq. 4.9
Mn2=(((Mn1^2)+(2/(gammac-1)))/(((2*gammac/(gammac-1))*(Mn1^2)) -1))^0.5;  %% Eq. 4.10
T2c=T1c*(p2c/p1c)*(rho1c/rho2c); %% Eq. 4.11
M2c=(Mn2/sin(betac-theta));  %% Eq. 4.12
a2c=sqrt(gammac*R*T2c);

% Normal velocity calculation 
Vn1=u1c*sin(betac); 
% un2=(rho1c*un1)/rho2c
Vn2=(Vn1*(tan(betac-theta))/tan(betac));
rz=0.1; % First estimative of Pho 2
rho2c=rho1c/rz; 

% Considering Vn,1/ Vn,2
p2= p1c + rho1c*(Vn1^2)*(1-(rho1c/rho2c)); % Eq. 2.40   Aanderson guessing
h2= h1c +((Vn1^2)/2)*(1-((rho1c/rho2c)^2)); % Eq. 2.42  Aanderson guessing
aa(1)=rho1c/rho2c;
it(1)=aa(1);
betann=0.5; % Arbitrated value for new Beta to run first 'while' loop
betan(1)=betac;

%%%%  PRINCIPAL LOOP (GLOBAL) - ITERATIVE PROCEDURE TO UPDATE THE BETA ANGLE %%%%%%%%

while ((betann)>=tol)
[d,YT,XT,ZT,X,A,B,C,rho2c,p2,h2,gamma,it,aa,x]= tannehilcalculo(tol,theta,p1c,h1c,rho1c,Vn1,h2,aa,it,rho2c,p2,ii,jj,T0);
betan(j+1)=x;
BETA=(180/pi())*(betan(j+1));
betacn=betan(j+1);
Vn1=u1c*sin(betacn); 
betann=abs(betan(j+1)-betan(j));
j=j+1;
end

% De posse de todos os parametros já convergidos, inclusive conhecendo-se o
% ângulo beta real da onda de choque é momento de realizar o calculo da
% tempratura
% [YT,XT,ZT,d]=tabelatanehill(rho2c,p2,ii,jj,T0)

T2=10^(( d(1) + d(2)*YT + d(3)*ZT + d(4)*YT*ZT + d(5)*(ZT^2) + ((d(6) + d(7)*YT + d(8)*ZT + d(9)*YT*ZT + d(10)*(ZT^2))/(1 + exp(d(11)*(ZT + d(12)))))));
T2=T2*T0;
Resultados=[BETA,p2,rho2c,T2] % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Foi encontrado o valor de pho1/pho2 baseando-se nas velocidades normais , de
%   posse dessa valor é momento de atualizar o ângulo beta, pois a convergência no passo 
%   anterior somente foi obtida graças utilizando o ângulo desatualizada.
%   O procedimento de atualização do ângulo é conforme o loop abaixo, e aplicando a Eq. 14.21
%   do livro Hypersonic High Temperature do Anderson. 
