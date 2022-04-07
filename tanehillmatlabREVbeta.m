% Hypersonics Aerothermodynamics
% Developed by Alisson Vinicius Brito Lopes - 18/11/2014
% Computing shock wave angle under chemical equilibrium conditions

M1c = 25; % Variar o Valor do número de Mach para Mach 4, 7 e 25
theta = 20; % DEVE SER INFORMADO O ÂNGULO DA CUNHA - neste caso 20 graus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1; i=1; ii=1; jj=1; % Contadores a serem utilizados durante os procedimento iterativos
tol=1e-6; % Serve como tolerância e calculo de erros nos laços "while"
theta=theta*(pi()/180); % Transformando o ângulo da cunha para radianos
gammac = 1.4; % Razão entre os calores especificos igual 1.4 para o caso caloricamente perfeito
b=0.001; % Variavel que serve como cálculo de erro em procedimento iterativo envolvendo o "while"
T0 = 151.78; % Temperatura de referência do paper do Tannehil
R=287.17; % Constante Equivalente do ar 
T1c = 226.650; % Temperatura do escoamento Livre para a Altitude de 30 Km
cp1c = 1006; % Calor especifico para o ar, valor dado no enunciado do problema J/kgK
a1c = sqrt(gammac*R*T1c); % Velocidade do Som para as condições de vôo informada 
rho1c = 1.8012e-2; % Densidade da corrente livre
p1c = 1171.8; % Pressão da corrente livre N/m²
h1c = cp1c*T1c; % Entalpia para corrente libre J/kg
u1c = M1c*a1c;% Velocidade da corrente livre  

%%%%%%%% CALCULO PARA O GÁS CALORICAMENTE PERFEITO %%%%%%%%%%%%%%%%%%%%
% Processo iterativo para calculo do angulo Beta caloricamente perfeito

beta = 18*(pi/180); % Chute inicial para o ângulo beta, Chutar um VALOR menor que o ângulo da cunha !!!
a=tan(theta); % Tangente do ângulo teta

while b <= a
b(i)=2*cot(beta)*(((M1c^2*(sin(beta))^2)-1)/(M1c^2*(gammac+cos(2*beta))+2)); % Equação a ser igualada por meio do novo beta
beta = beta + 0.0001; % Incremento no ângulo para obter novo beta
i=i+1;
end

betac=beta;% ãngulo beta em radianos
Abeta=betac*(180/pi()); % Angulo beta em graus
Mn1=M1c*sin(betac); %% Número de Mach Normal a onda - Eq. 4.7 Anderson (compressible Flow)
rho2c=rho1c*(((gammac+1)*(Mn1^2))/(((gammac-1)*(Mn1^2)) +2)); %% Eq. 4.8
p2c=p1c*(1+(2*gammac/(gammac+1))*((Mn1^2)-1));  %% Eq. 4.9
Mn2=(((Mn1^2)+(2/(gammac-1)))/(((2*gammac/(gammac-1))*(Mn1^2)) -1))^0.5;  %% Eq. 4.10
T2c=T1c*(p2c/p1c)*(rho1c/rho2c); %% Eq. 4.11
M2c=(Mn2/sin(betac-theta));  %% Eq. 4.12
a2c=sqrt(gammac*R*T2c);

% Calculo da Velocidade Normal ao Shock
Vn1=u1c*sin(betac); % Velocidade
% un2=(rho1c*un1)/rho2c
Vn2=(Vn1*(tan(betac-theta))/tan(betac));
rz=0.1; % Primeira estimativa do Pho 2
rho2c=rho1c/rz; % Primeira estimativa do rho2

% Considerando a Velocidade Normal Vn,1/ Vn,2
p2= p1c + rho1c*(Vn1^2)*(1-(rho1c/rho2c)); % Eq. 2.40   ANDERSON chute
h2= h1c +((Vn1^2)/2)*(1-((rho1c/rho2c)^2)); % Eq. 2.42  ANDERSON chute
aa(1)=rho1c/rho2c;
it(1)=aa(1);
betann=0.5; % Valor arbitrado para o beta novo, afim de rodar a primeira iteração no laço "while"
betan(1)=betac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOOP PRINCIPAL (GLOBAL) - PROCEDIMENTO ITERATIVO PARA CALCULO %%%%%%%%
%%%%%%%%%%%%%%%%%%% DE ATUALIZAÇÃO DO ÂNGULO BETA %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ((betann)>=tol)
[d,YT,XT,ZT,X,A,B,C,rho2c,p2,h2,gamma,it,aa,x]= tannehilcalculo(tol,theta,p1c,h1c,rho1c,Vn1,h2,aa,it,rho2c,p2,ii,jj,T0);
betan(j+1)=x;
BETA=(180/pi())*(betan(j+1));
betacn=betan(j+1);
Vn1=u1c*sin(betacn); % Velocidade
betann=abs(betan(j+1)-betan(j));
j=j+1;
end

% De posse de todos os parametros já convergidos, inclusive conhecendo-se o
% ângulo beta real da onda de choque é momento de realizar o calculo da
% tempratura
% [YT,XT,ZT,d]=tabelatanehill(rho2c,p2,ii,jj,T0)

T2=10^(( d(1) + d(2)*YT + d(3)*ZT + d(4)*YT*ZT + d(5)*(ZT^2) + ((d(6) + d(7)*YT + d(8)*ZT + d(9)*YT*ZT + d(10)*(ZT^2))/(1 + exp(d(11)*(ZT + d(12)))))));
T2=T2*T0;
Resultados=[BETA,p2,rho2c,T2] % VETOR QUE IMPRIME OS RESULTADOS FINAIS DESTE TRABALHO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Foi encontrado o valor de pho1/pho2 baseando-se nas velocidades normais , de
%   posse dessa valor é momento de atualizar o ângulo beta, pois a convergência no passo 
%   anterior somente foi obtida graças utilizando o ângulo desatualizada.
%   O procedimento de atualização do ângulo é conforme o loop abaixo, e aplicando a Eq. 14.21
%   do livro Hypersonic High Temperature do Anderson. 
