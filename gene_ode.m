function dy = gene_ode(t,y)
% y = [ J_ww(1) , J_dw(2), J_dd(3), A_ww(4) , A_dw(5), A_dd(6)]
global beta r x f delta mu_0 mu_1 mu_2 alpha

%initial matrix
dy = zeros(6,1);
%Goodwin model for
dy(1) =  beta*(1-r)*r*( y(4)^2 + x*(f +1)*y(5)*y(4) + (x^2)*f*y(5)^2)/(r*(y(4)+ y(5)+y(6)) + delta)- ( mu_1 + mu_0*y(1))*y(1) - alpha*y(1);

dy(2) = (beta*(1-r)*r/(r*(y(4)+ y(5)+y(6)) + delta) )*((1-x)*(1+f)*y(5)*y(4) +y(6)*y(4) + 2*x*(1-x)*f*y(5)^2 + x*f*y(5)*y(6))  - (mu_1 + mu_0*y(2))*y(2) - alpha*y(2);
dy(3) = ((1-x)*f*beta*(1-r)*r*y(5)/(r*(y(4)+ y(5)+y(6)) + delta))*( (1-x)*y(5)  + y(6) ) - (mu_1 + mu_0*y(3))*y(3) - alpha*y(3);

dy(4) = alpha*y(1) - mu_2*y(4);
dy(5) = alpha*y(2) - mu_2*y(5);
dy(6) = alpha*y(3) - mu_2*y(6);
end