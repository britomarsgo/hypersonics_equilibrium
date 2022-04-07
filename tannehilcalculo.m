%%%% Função TANNEHILL e MUGGE 

function [d,YT,XT,ZT,X,A,B,C,rho2c,p2,h2,gamma,it,aa,x]= 
tannehilcalculo(tol,theta,p1c,h1c,rho1c,Vn1,h2,aa,it,rho2c,p2,ii,jj,T0)

while (aa (ii) >=tol)

[c,X,Y,Z,YT,XT,ZT,d]=tabelatanehill(rho2c,p2,ii,jj,T0);
gamma = c(1) + c(2)*Y(ii) + c(3)*Z(ii) + c(4)*Y(ii)*Z(ii) + ((c(5) + 
c(6)*Y(ii) + c(7)*Z(ii) + c(8)*Y(ii)*Z(ii))/(1 + exp(c(9)*(X(ii) + 
c(10)*Y(ii) + c(11)))));
rho2(ii) = (p2/h2)*(gamma/(gamma-1)) ;% 
p2= p1c + rho1c*((Vn1)^2)*(1-(rho1c/rho2(ii))) ;% 
h2= h1c +(((Vn1)^2)/2)*(1-((rho1c/rho2(ii))^2)) ;
rho2c=rho2(ii);
it(ii+1)=rho1c/rho2(ii);
aa(ii+1)=abs(it(ii+1)-it(ii));
ii=ii+1;
end

% Resolvendo a equação do segundo Grau A*tan(beta)² + B*tan(beta) + C 
% = 0, onde C=tan(theta)

% Algoritmo (Método) de baskara foi utlizado.

X= rho1c/rho2c; % Razão entre densidades reais, calculada pelo metodo 
de tannehill
A=X*(tan(theta));
B=X-1;
C=tan(theta) ;
% Xbeta=atan((-B -sqrt((B^2) -4*A*C))/(2*A)) + pi()*n
delta=(B^2) - 4*A*C;
x1 = (-B + sqrt(delta))/(2*A);
x2= (-B -sqrt(delta))/(2*A);
tanbetanew=[x1,x2] ;
tanbetanew=(180/pi())*atan(tanbetanew);

    if x1>x2 % Esse Loop escolhe sempre a menor raiz, que corresponde ao 
    % valor próximo ao valor de beta caloricamente perfeito, o ãngulo muito 
    % grande não tem sentido físico
    x=atan(x2);
    end
end