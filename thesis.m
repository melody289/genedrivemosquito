%% Figure 5.2a

% Anopheles parameters
beta = 33.3333; 
r = .5;
%x = Parm3(3,i);
%f =  Parm3(4,i);
mu_1 = .2;
mu_0 = 5e-2;
alpha = 0.1020;
mu_2 = 0.0435;
delta2 = 200:762;
    

% Just putting the equilibriums in pieces
aa = (- mu_0*mu_2^2 ).*delta2 - alpha^2*mu_2 *r - alpha* mu_1*mu_2 *r + alpha^2*r*(1-r)*beta ;
   


bb =  sqrt((  (mu_0*mu_2^2).*delta2 + alpha^2* mu_2 *r +alpha*mu_1*mu_2* r - alpha^2*r*(1-r)*beta).^2 - (4*r*alpha).*delta2.*(mu_0*mu_2^3*(mu_1 + alpha)) );

dd = (2*mu_0*mu_2^2*r);

% This are the two equilibria for a given delta value
eq1a = (aa + bb)./dd;
  
eq2a = (aa - bb)./dd;

%
plot(delta2, eq1a,'b', delta2, eq2a,'--b')

xlabel('Allee Constant: \delta' ,'FontSize', 25)

ylabel('Equilibrium Value' ,'FontSize', 25)
hold on
plot([200 1400], [0 0], 'b','LineWidth', 1.4)

%
axis([200,1000, -150, 1500])
hold on
plot(762, 157.7488, '*')
%
legend('Stable','Unstable', 'Bifurcation')

%% Figure 5.2b and Latin hypercube sampling
% Now we choose a sample size
samp = 225;
% This is the latin hypercube command in matlab. When I run it in MatLab
% R2016a, I use the default kernel when I first open matlab.

% Here is the seed

s = RandStream('mt19937ar','Seed',0, 'NormalTransform', 'Ziggurat');
RandStream.setGlobalStream(s)

Parm = lhsdesign(7,samp);


% % Order is beta, r, mu0, mu1, mu2, alpha, delta

% This is my lower and upper bounds
lb = [ 50/4, 0.5,  log(1e-6),  0.05, 0.0333, 0.0333, 50 ]';
ub = [ 122/2, 0.55, log(.1),  0.8,  0.125, 0.125, 250 ]';

% Adjusting random sample to fit bounds
Parm200 = bsxfun(@plus, lb, bsxfun(@times, Parm, (ub-lb)));

% Getting rid of log on mu_0
Parm200(3,:) = exp(Parm200(3,:));




% First I need to refind all equilibriums
eq1 = zeros(225, 1);
eq2 = zeros(225, 1);
for i = 1:samp
    
beta = Parm200(1,i);
r = Parm200(2,i);

mu_0 = Parm200(3,i);
mu_1 = Parm200(4,i);
mu_2 = Parm200(5,i);

alpha =Parm200(6,i);
delta =   Parm200(7,i);
aa = - mu_0*mu_2^2 *delta - alpha^2*mu_2 *r - alpha* mu_1*mu_2 *r + alpha^2*r*(1-r)*beta ;
   


bb =  sqrt((  mu_0*mu_2^2 *delta + alpha^2* mu_2 *r +alpha*mu_1*mu_2* r - alpha^2*r*(1-r)*beta)^2 - 4*r*alpha*delta*mu_0*mu_2^3*(mu_1 + alpha) );

dd = (2*mu_0*mu_2^2*r);


eq1(i) = (aa + bb)/dd;
  
eq2(i) = (aa - bb)/dd;
end



% Finding positive equilibria
kpos = 0;
kimag = 0;
kneg = 0;

for i = 1:samp
  if( isreal(eq1(i) ) )
   if(  eq1(i)> 0 )
       kpos = [kpos, i ];
   elseif (eq1(i) < 0)
       kneg = [kneg, i ] ;
   end
  else
      kimag = [kimag,i];
  end
end
kneg(1) = [];
kpos(1) = [];
kimag(1) = [];

%
Parm = Parm200(:,kpos);
eq1 = eq1(kpos);
eq2 = eq2(kpos);

% These are the ones I will not use, because they have negative and
% imaginary equilbria

Parmnegi = Parm200(:,[kneg, kimag]);


%
bif = zeros(1,192);
del = zeros(1,192);
% I am finding the bifurcation point for each set of parameters, but solving for delta value.
for i =1:192 
% beta   
beta = Parm(1,i);
% r
r = Parm(2,i);
% mu0
mu_0 = Parm(3,i);
%mu1
mu_1 = Parm(4,i);
% mu2
mu_2 = Parm(5,i);
% alpha
alpha = Parm(6,i);

a1 = (  mu_0*mu_2^2)^2;
b1 = 2*mu_0*mu_2^2*(  alpha^2* mu_2 *r +alpha*mu_1*mu_2* r - alpha^2*r*(1-r)*beta)- 4*r*alpha*mu_0*mu_2^3*(mu_1 + alpha);
c1 = (  alpha^2* mu_2 *r +alpha*mu_1*mu_2* r - alpha^2*r*(1-r)*beta)^2;
del(i) = (-b1 + sqrt(b1^2 -4*a1*c1))/ (2*a1);

bb1 = - mu_0*mu_2^2 *del(i) - alpha^2*mu_2 *r - alpha* mu_1*mu_2 *r + alpha^2*r*(1-r)*beta;
ac1 = sqrt((  mu_0*mu_2^2 *del(i) + alpha^2* mu_2 *r +alpha*mu_1*mu_2* r - alpha^2*r*(1-r)*beta)^2 - 4*r*alpha*del(i)*mu_0*mu_2^3*(mu_1 + alpha) );
d = 2*mu_0*mu_2^2*r;

if(abs(ac1)< 1e-4)
bif(i) = -bb1/d;
else
    i
end

end


% Plotting bifurcation with delta value
plot(log10(del),log10(bif), 'o') 
hold on 
%title(' ', 'FontSize', 20)
ylabel('Bifurcation Value' ,'FontSize', 20)

xlabel('Allee Constant: \delta' ,'FontSize', 20)

% axis([ 0 1274, 0 0.5])

% Finding linear regression line of log values
xr = [ ones(192,1), log10(del)'];
line = xr\log10(bif)';

xr2 = [ ones(192,1), del'];
line2 = xr2\bif';

%
hold on
t = [2:8];
plot(t, line(1) + line(2).*t)

%% Figure 5.2c

hold on
plot(log10(Parm(3,:)), log10(bif), 'o') 
hold on 
%title(' ', 'FontSize', 20)
ylabel('Bifurcation Value' ,'FontSize', 20)

xlabel('Density dependent death: \mu_0' ,'FontSize', 20)

% axis([ 0 1274, 0 0.5])

% Finding linear regression line of log values
xr3 = [ ones(192,1), log10(Parm(3,:))'];
line3 = xr3\log10(bif)';

xr4 = [ ones(192,1), Parm(3,:)'];
line4 = xr4\bif';
hold on
t = [-6,-1];
y = line3(1) + line3(2).*t;
plot(t, y)

%% Figure 5.3
% This is the step size I will use for the histogram
stepn = (ub-lb)./20;
for i = 1:7
edges(i,:) = lb(i):stepn(i):ub(i);
end

% mu0 is on a log scale so I do this separately
 step = (log10(ub(5))+6 )/20 ;
stepup= -6: step:log10(ub(5));
stepup = 10.^(stepup);
edges5 = stepup;
%%
 ss2 = zeros(7,2);
  words = { ' Egg laying Rate', ' Male Mosquitoes' ,   'Density Dependent Death Rate of Larvae' ,  'Constant Death Rate of Larvae' ,  'Adult Death Rate', 'Transition Rate from Larvae to Adult' , 'Allee Constant'   };
% % % Order is beta, r, mu0, mu1, mu2, alpha, delta
  xwords =   { '\beta', 'r' ,   '\mu_0' , '\mu_1' ,  '\mu_2' , '\alpha' , '\delta'  };

for i = 1:2



 

 hh = figure


%subplot(1,2,1)  
h1 = histogram(Parm(i,:),edges(i,:));
hold on
  h2= histogram(Parmnegi(i,:),edges(i,:));




% This was to make it easier to compare to normalize the data
h1.Normalization = 'probability';
h2.Normalization = 'probability';

figure
stairs(h1.Values, 'LineWidth', 2)
hold on
stairs(h2.Values, 'LineWidth', 2)

title( words(i), 'Fontsize', 12);
legend('Positive Equilibria', 'Negative or Imaginary Equilibria');
xlabel( xwords(i), 'Fontsize', 20)
ylabel( 'Proportion', 'Fontsize', 15)

%set(gca, 'XTick',  [5:4:21], 'XTickLabel',  [  '10^{-5}','10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}' ])
%set(gca, 'XTick',  [1.5:4:21], 'XTickLabel',  [ .03:.02:.11])
set(gca, 'FontSize', 15);
close hh
end

hh= figure
i =3;

%subplot(1,2,1)  
h1 = histogram(Parm(i,:),edges5);
hold on
  h2= histogram(Parmnegi(i,:),edges5);



% This was to make it easier to compare to normalize the data
h1.Normalization = 'probability';
h2.Normalization = 'probability';


set(gca,'xscale','log') 
 figure
stairs(h1.Values, 'LineWidth', 2)
hold on
stairs(h2.Values, 'LineWidth', 2)


title( words(i), 'Fontsize', 12);
legend('Positive Equilibria', 'Negative or Imaginary Equilibria');
xlabel( xwords(i), 'Fontsize', 20)
ylabel( 'Proportion', 'Fontsize', 15)


set(gca, 'FontSize', 15);
close hh

for i = 4:7



 

 hh = figure


%subplot(1,2,1)  
h1 = histogram(Parm(i,:),edges(i,:));
hold on
  h2= histogram(Parmnegi(i,:),edges(i,:));




% This was to make it easier to compare to normalize the data
h1.Normalization = 'probability';
h2.Normalization = 'probability';

figure
stairs(h1.Values, 'LineWidth', 2)
hold on
stairs(h2.Values, 'LineWidth', 2)

title( words(i), 'Fontsize', 12);
legend('Positive Equilibria', 'Negative or Imaginary Equilibria');
xlabel( xwords(i), 'Fontsize', 20)
ylabel( 'Proportion', 'Fontsize', 15)

%set(gca, 'XTick',  [5:4:21], 'XTickLabel',  [  '10^{-5}','10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}' ])
%set(gca, 'XTick',  [1.5:4:21], 'XTickLabel',  [ .03:.02:.11])
set(gca, 'FontSize', 15);
close hh
end




%% Figure 5.4a
% Note: this takes some over 15 minutes to run

%  Matrix for samples without Allee effect of Latin hypercube 190 samples
% Set up
% For IV

in = (eq2./eq1);
e = zeros(192,1);
for i = 1:192
    
  beta = Parm(1,i);
r = Parm(2,i);

mu_0 = Parm(3,i);
mu_1 = Parm(4,i);
mu_2 = Parm(5,i);

alpha =Parm(6,i); 
    

% This is the equilibrium when there is no Allee effect i.e. delta = 0
   e(i) =  ((1-r)*alpha^2*beta - mu_1*mu_2*alpha - mu_2*alpha^2) /(mu_2^2*mu_0);
    
    
end

% Making sure that the Allee threshold/Carrying capacity < 15 percent
[rel, place2] = sort(in);

[noa, place] = sort(e);
toobig(1) = find(place==place2(end));
toobig(2) = find(place==place2(end-1));
%
% Removing the two too large Relative Allee effects
place(toobig) = [];


% Matrix for samples
xx  = linspace(.01,.2,15);
zz = zeros(15,192);

%
for j = 1:15;
tic
for i= 1:192;
%187
iV = [ 0,0, 0,  eq1(i)*in(i)^.05 , 0, 10];
%parm = [ 46.7863    0.5124    0.01    0.4311  158.0851    0.5590    0.3337    0.0987    0.0399]

 % beta   
parm(1) = Parm(1,i);
% r
parm(2) = Parm(2,i);
% mu0
parm(3) = Parm(3,i);
%mu1
parm(4) = Parm(4,i);
% mu2
parm(5) = Parm(5,i);
% alpha
parm(6) = Parm(6,i);
% delta
parm(7) =  Parm(7,i);
% x
parm(8) = xx(j);
% f
parm(9) = 0.5;
% fd
parm(10) = 0;

tspan = [0,5000];
%solving the ODE 
options = odeset('NonNegative',1:6);


[t,y] = ode45(@gene_ode,tspan, iV, options,parm);
%semilogy(t,sum(y(:,4:6),2), 'LineWidth', 2);
%hold on
%legend('J_{WW}','J_{DW}','J_{DD}','A_{WW}','A_{DW}','A_{DD}','Location','South')
 zz(j,i) = sum(y(end,4:6));
end
toc
end

%

% I log scale it so that the color scale is easier to read

zzl = log10(zz);

%

figure
imagesc(zzl(:,place))
c = colorbar;
caxis([0 , 5.5])
%
% converting it from log scale
newtic = 0:5;
newtic = 10.^newtic;

set(c,'YTick',0:5,'YTickLabel', newtic)

%
 xti = ones(1,15) - round(xx, 2);
k = [1:47:188, 190];
set(gca, 'XTick',  k, 'XTickLabel',  round(e(place(k)))) % 10 ticks 
set(gca, 'YTick', 1:3:15, 'YTickLabel', xti(1:3:15) ) % 20 ticks
title('No Allee Effect', 'FontSize', 20)
%xlabel(' Allee Threshold/ Carrying Capacity' ,'FontSize', 18)

xlabel('Equilibria' ,'FontSize', 20, 'FontWeight', 'bold')
colormap copper
ylabel('Gene Drive' ,'FontSize', 20, 'FontWeight', 'bold')


%% Figure 5.4b
% Note: this takes some over 15 minutes to run

% Matrix for samples WITH Allee effect of Latin hypercube 190 samples

zz2 = zeros(15,192);



%
for j = 1:15;
tic
for i= 1:192;
%187
iV = [ 0,0, 0,  eq1(i)*in(i)^.05 , 0, 10];
%parm = [ 46.7863    0.5124    0.01    0.4311  158.0851    0.5590    0.3337    0.0987    0.0399]

 % beta   
parm(1) = Parm(1,i);
% r
parm(2) = Parm(2,i);
% mu0
parm(3) = Parm(3,i);
%mu1
parm(4) = Parm(4,i);
% mu2
parm(5) = Parm(5,i);
% alpha
parm(6) = Parm(6,i);
% delta No Allee effect
parm(7) =  0;
% x
parm(8) = xx(j);
% f
parm(9) = 0.5;
% fd
parm(10) = 0;

tspan = [0,5000];
%solving the ODE 
options = odeset('NonNegative',1:6);


[t,y] = ode45(@gene_ode,tspan, iV, options,parm);
%semilogy(t,sum(y(:,4:6),2), 'LineWidth', 2);
%hold on
%legend('J_{WW}','J_{DW}','J_{DD}','A_{WW}','A_{DW}','A_{DD}','Location','South')
 zz2(j,i) = sum(y(end,4:6));
end
toc
end

%

% I log scale it so that the color scale is easier to read

zzl2 = log10(zz2);

%

figure
imagesc(zzl2(:,place))
c = colorbar;
caxis([0 , 5.5])
%
% converting it from log scale
newtic = 0:5;
newtic = 10.^newtic;

set(c,'YTick',0:5,'YTickLabel', newtic)

%
 xti = ones(1,15) - round(xx, 2);
k = [1:47:188, 190];
set(gca, 'XTick',  k, 'XTickLabel',  round(e(place(k)))) % 10 ticks 
set(gca, 'YTick', 1:3:15, 'YTickLabel', xti(1:3:15) ) % 20 ticks
title('Allee Effect', 'FontSize', 20)
%xlabel(' Allee Threshold/ Carrying Capacity' ,'FontSize', 18)

xlabel('Equilibria' ,'FontSize', 20, 'FontWeight', 'bold')
colormap copper
ylabel('Gene Drive' ,'FontSize', 20, 'FontWeight', 'bold')


%% Figure 5.5a and b
% Note this takes over 15 minutes to run

%Intitializing
zxf2 = zeros(11,11);
zxfA2 = zeros(11,11);
tic
for i = 1:11;
    
for j = 1:11;
 
    
%%% order of parameter
%%%% beta r mu_0 mu_1 mu_2 alpha delta x f fd 
%parm = [33.3333, 0.5, 5e-2, 0.2, 0.0435,0.1020, 0, 0.5, 0.5, 0];
parm = [ 33.3333    0.5    0.05    0.2  0.0435      0.102  0  0.5  .2 0];

iV = [ 0,0, 0, 50 , 0, 10];

parm(8) = xg(i);
parm(9) = f(j); 





tspan = [0,5000];
%solving the ODE 
options = odeset('NonNegative',1:6);

[t,y] = ode45(@gene_ode,tspan, iV, options,parm);

%
%semilogy(t,y, 'LineWidth', 2);

 zxf2(i,j) = sum(y(end,4:6));
end

end
toc

for i = 1:11;
    tic
for j = 1:11;
 
%%% order of parameter
%%%% beta r mu_0 mu_1 mu_2 alpha delta x f fd 
%parm = [33.3333, 0.5, 5e-2, 0.2, 0.0435,0.1020, 0, 0.5, 0.5, 0];
parm = [ 33.3333    0.5    0.05    0.2  0.0435      0.102  150  0.5  .2 0];

iV = [ 0,0, 0, 50 , 0, 10];

parm(8) = xg(i);
parm(9) = f(j); 





tspan = [0,500];
%solving the ODE 
options = odeset('NonNegative',1:6);

[t,y] = ode45(@gene_ode,tspan, iV, options,parm);

%
%semilogy(t,y, 'LineWidth', 2);

 zxfA2(i,j) = sum(y(end,4:6));
end
toc
end

%%

figure
imagesc(zxf2)
c = colorbar;
ylabel(c,'Population at Equilbrium','FontSize',20)

%
set(gca, 'XTick',  1:2:11, 'XTickLabel',  f(1:2:11)) % 10 ticks 
set(gca, 'YTick', [1:2:11], 'YTickLabel', 1-xg([1:2:11]) )

xlabel('Fertility of Heterozygote' ,'FontSize', 15)

ylabel('Gene Drive' ,'FontSize', 15)
title('Varying fertility and Gene Drive' ,'FontSize', 15)
colormap bone
hold on
plot([ 1.5,1.5,2.5, 2.5, 3.5, 3.5, 4.5, 4.5, 5.5, 5.5, 6.5,6.5,7.5,7.5, 9.5, 9.5, 10.5, 10.5, 11.5],[ .5,2.5, 2.5, 4.5, 4.5, 5.5, 5.5, 6.5, 6.5, 7.5,7.5,8.5,8.5, 9.5, 9.5, 10.5, 10.5, 10.5, 10.5], 'w','LineWidth', 5)
plot([1.5,11.5],[1.5, 1.5],'w','LineWidth', 5)

hold on
plot([ 1.5,1.5,2.5, 2.5, 3.5, 3.5, 4.5, 4.5, 5.5, 5.5, 6.5,6.5,7.5,7.5, 9.5,  11.5],[ .5,2.5, 2.5, 4.5, 4.5, 4.5, 4.5, 5.5, 5.5, 5.5, 5.5,5.5, 5.5, 4.5, 4.5, 4.5], ':w','LineWidth', 5)



% Now we find the difference with Allee
difA = (zxf2 - zxfA2)./zxf2;

figure
imagesc(difA)
c = colorbar;
ylabel(c,'Population Difference','FontSize',20)

%
set(gca, 'XTick',  1:2:11, 'XTickLabel',  f(1:2:11)) % 10 ticks 
set(gca, 'YTick', [1:2:11], 'YTickLabel', 1-xg([1:2:11]) )

xlabel('Fertility of Heterozygote' ,'FontSize', 15)

ylabel('Gene Drive' ,'FontSize', 15)
title('Change in Population due to Allee Effect' ,'FontSize', 15)

colormap copper
hold on
plot([ 1.5,1.5,2.5, 2.5, 3.5, 3.5, 4.5, 4.5, 5.5, 5.5, 6.5,6.5,7.5,7.5, 9.5, 9.5, 10.5, 10.5, 11.5],[ .5,2.5, 2.5, 4.5, 4.5, 5.5, 5.5, 6.5, 6.5, 7.5,7.5,8.5,8.5, 9.5, 9.5, 10.5, 10.5, 10.5, 10.5], ':w','LineWidth', 5)
plot([1.5,11.5],[1.5, 1.5],':k','LineWidth', 5)

hold on
plot([ 1.5,1.5,2.5, 2.5, 3.5, 3.5, 4.5, 4.5, 5.5, 5.5, 6.5,6.5,7.5,7.5, 9.5,  11.5],[ .5,2.5, 2.5, 4.5, 4.5, 4.5, 4.5, 5.5, 5.5, 5.5, 5.5,5.5, 5.5, 4.5, 4.5, 4.5], 'k','LineWidth', 5)


%% Figure 5.6

% All lines were found in mathematica, but plotted here
figure
subplot(1,3,1)
plot([-.05, 1.05], [.25, .25], 'k', 'LineWidth', 3)
hold on 
axis([ -.05, 1.05, -.05, 1.05])
plot( [.25, 1.05],[.25, 1.05], 'k', 'LineWidth', 3)
title(' x = 0.5')
ylabel('f_d: Fertility of Homozygous Gene Drive')

set(gca, 'FontSize', 18)
%
subplot(1,3 ,2)
plot([-.05, 1.05], [.25, .25], 'k', 'LineWidth', 3)
hold on 
axis([ -.05, 1.05, -.05, 1.05])
plot( [.25, .25],[-.05, 1.05], 'k', 'LineWidth', 3)
title(' x = 0.8')
xlabel('f: Fertility of Heterozygotes')
set(gca, 'FontSize', 18)

subplot(1,3 ,3)
plot([-.05, 1.05], [.25, .25], 'k', 'LineWidth', 3)
hold on 
axis([ -.05, 1.05, -0.05, 1.05])
plot( [0.0526, 0.0526],[-.05, 1.05], 'k', 'LineWidth', 3)

title(' x = 0.95')

set(gca, 'FontSize', 18)

%% Figure 5.7

samp = 11; 

ff = linspace(0,1,samp); % values 0:.1:1 for f 
fdd = linspace(1,0,samp); % values 0:.1:1 for fd 
zziv = zeros(samp,samp);
zziv2 = zeros(samp,samp);

for i = 1:samp;

for j = 1:samp

%%% Initial conditions
mos = 2000;
test = 10;
iV = [ 0,0, 0, mos-test , 0, test];

%%% Order of parameter
%%%% is beta, r, mu0, mu1, mu2, alpha, delta, x, f, fd 
parm = [33.3333, 0.5, 0.05, 0.2, 0.0435,0.1020, 150, 0.5, 0.5, 0];


%%% Over-writing some of the parameters that are varied
parm(8) = .8;  % This is gene drive, if you want 90, make it =0.9
parm(9) = ff(j);% The heterozygous fertiility f
parm(10) = fdd(i); % The homozygous fertility fd


 tspan = [0 10000]; % use this for better accuracy, allows matlab to chose time step
%tspan = 0:10000;% use to run faster (note may have to choose values > 2000 for some situations)

%%% solving the ODE 
% options = odeset('NonNegative',1:6); % use this for best accuracy but
% much slower
options = [];%odeset('NonNegative',1:6); % use this for faster but less accurate near zero


[t,y] = ode45(@gene_ode,tspan, iV, options,parm);


zziv(i,j) = sum(y(end,4:6));

end

end


for i = 1:samp;

for j = 1:samp

%%% Initial conditions
mos = 10;
test = 2000;
iV = [ 0,0, 0, mos-test , 0, test];

%%% Order of parameter
%%%% is beta, r, mu0, mu1, mu2, alpha, delta, x, f, fd 
parm = [33.3333, 0.5, 0.05, 0.2, 0.0435,0.1020, 150, 0.5, 0.5, 0];


%%% Over-writing some of the parameters that are varied
parm(8) = .8;  % This is gene drive, if you want 90, make it =0.9
parm(9) = ff(j);% The heterozygous fertiility f
parm(10) = fdd(i); % The homozygous fertility fd


 tspan = [0 10000]; % use this for better accuracy, allows matlab to chose time step
%tspan = 0:10000;% use to run faster (note may have to choose values > 2000 for some situations)

%%% solving the ODE 
% options = odeset('NonNegative',1:6); % use this for best accuracy but
% much slower
options = [];%odeset('NonNegative',1:6); % use this for faster but less accurate near zero


[t,y] = ode45(@gene_ode,tspan, iV, options,parm);


zziv2(i,j) = sum(y(end,4:6));

end

end

imagesc(zziv)
c = colorbar;
caxis([1 , 1650])
 %colormap(map);
set(gca,'fontsize',14,'xtick',1:11,'xticklabel',0:.1:1,'ytick',1:11,'yticklabel',fliplr(0:.1:1))

ylabel('f_d: Fertility of Homozygous Gene Drive')
xlabel('f: Fertility of Heterozygotes')

title('Initial Amount of Wild Mosquitoes greater than Gene Drive')
ylabel(c,'Population at Equilibrium','FontSize',20)
hold on
plot([-.05, 1.05], [.25, .25], 'w', 'LineWidth', 3)
hold on 
axis([ -.05, 1.05, -.05, 1.05])
plot( [.25, .25],[-.05, 1.05], 'k', 'LineWidth', 3)

figure
imagesc(zziv2)
c = colorbar;
caxis([1 , 1650])
 %colormap(map);
set(gca,'fontsize',14,'xtick',1:11,'xticklabel',0:.1:1,'ytick',1:11,'yticklabel',fliplr(0:.1:1))

ylabel('f_d: Fertility of Homozygous Gene Drive')
xlabel('f: Fertility of Heterozygotes')

title('Initial Amount of Gene Drive Mosquitoes greater than wild')

ylabel(c,'Population at Equilibrium','FontSize',20)
hold on
plot([-.05, 1.05], [.25, .25], 'w', 'LineWidth', 3)
hold on 
axis([ -.05, 1.05, -.05, 1.05])
plot( [.25, .25],[-.05, 1.05], 'w', 'LineWidth', 3)
