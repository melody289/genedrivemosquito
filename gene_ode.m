function dy = gene_ode(t,y, val)
% y = [ J_ww(1) , J_dw(2), J_dd(3), A_ww(4) , A_dw(5), A_dd(6)]
beta = val(1);
r = val(2);
mu_0 = val(3);
mu_1 = val(4);
mu_2 = val(5);
alpha = val(6);
delta = val(7);
x = val(8);
f =  val(9);
fd = val(10);

%initial matrix
dy = zeros(6,1);
%
dy(1) =  beta*(1-r)*r*( y(4)^2 + x*(f +1)*y(5)*y(4) + (x^2)*f*y(5)^2)/(r*(y(4)+ y(5)+y(6)) + delta) - (mu_1 + alpha )*y(1) - mu_0*(y(1)+y(2)+y(3))*y(1);

dy(2) = (beta*(1-r)*r/(r*(y(4)+ y(5)+y(6)) + delta) )*((1-x)*(1+f)*y(5)*y(4) +(1+fd)*y(6)*y(4) + 2*x*(1-x)*f*y(5)^2 + x*(f+fd)*y(5)*y(6))  - (mu_1 + alpha )*y(2) - mu_0*(y(1)+y(2)+y(3))*y(2);
dy(3) = (beta*(1-r)*r/(r*(y(4)+ y(5)+y(6)) + delta) )*( f*(1-x)^2*y(5)^2  + (1-x)*(f+fd)*y(5)*y(6)+ fd*y(6)^2 )  - (mu_1 + alpha )*y(3) - mu_0*(y(1)+y(2)+y(3))*y(3);

dy(4) = alpha*y(1) - mu_2*y(4);
dy(5) = alpha*y(2) - mu_2*y(5);
dy(6) = alpha*y(3) - mu_2*y(6);
end
