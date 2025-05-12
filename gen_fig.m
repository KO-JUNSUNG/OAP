% Fig 1 3 anchors 
% Fig 2 7 acnhors
close all
N=3;
xy=exp(1j*2*pi*(0:N-1)/N);
x=real(xy);
y=imag(xy);
x=x/norm(x);
y=y/norm(y);

figure(1)
plot([0 x],[0 y], 'b*','LineWidth',2)
hold on
plot([-1,1], [0,0],':k')
plot([0,0],[-1,1],':k')
r=sqrt(x(1)^2+y(1)^2);
 plot(r*exp(1j*linspace(0,2*pi)),':k')
hold off
axis equal
axis ([-1 1 -1 1])
xlabel('x')
ylabel('y')
title('Optimal anchor placement')
%%
figure(2)

N=3;
xy=exp(1j*2*pi*(0:N-1)/N);
x1=0.5*real(xy);
y1=0.5*imag(xy);

N=4;
xy=exp(1j*2*pi*(0:N-1)/N);
x2=real(xy);
y2=imag(xy);

x=[x1 x2];
y=[y1 y2];
x=x/norm(x);
y=y/norm(y);
r1=sqrt(x(1)^2+y(1)^2);
r2=sqrt(x(end)^2+y(end)^2);
plot([0 x], [0,y],'b*','LineWidth',2)
hold on
plot([-1,1], [0,0],':k')
plot([0,0],[-1,1],':k')
 plot(r1*exp(1j*linspace(0,2*pi)),':k')
 plot(r2*exp(1j*linspace(0,2*pi)),':k')
hold off
xlabel('x')
ylabel('y')

axis equal
axis ([-1 1 -1 1])
title('Optimal anchor placement')
%%
figure(3)
theta=0.5*asin(2-sqrt(3));
phi=pi/2-theta;

th_opt=theta;
ph_opt=phi;

x_opt=[cos(th_opt) sin(th_opt)];
y_opt=[cos(ph_opt) sin(ph_opt)];



plot([0 x_opt], [0,y_opt],'b*','LineWidth',2)
hold on
plot([-1,1], [0,0],':k')
plot([0,0],[-1,1],':k')
 plot(exp(1j*linspace(0,2*pi)),':k')
hold off

axis equal
axis ([-1 1 -1 1])

xlabel('x')
ylabel('y')

title('Optimal anchor placement')

%minimum of cost 
x1=x_opt(1);
x2=x_opt(2);
y1=y_opt(1);
y2=y_opt(2);

ip=x1.*y1+x2.*y2;
cost_min=(((x1+x2)-(y1+y2).*ip).^2+((y1+y2)-(x1+x2).*ip).^2 +2*(1-ip.^2))./(1-ip.^2).^2;



[th, ph]=meshgrid(linspace(-pi,pi));
    
x1=cos(th);
x2=sin(th);
y1=cos(ph);
y2=sin(ph);
ip=x1.*y1+x2.*y2;
cost=(((x1+x2)-(y1+y2).*ip).^2+((y1+y2)-(x1+x2).*ip).^2 +2*(1-ip.^2))./(1-ip.^2).^2;
for k=1:length(th)
    cost(k,k)=nan;
end
       

figure(4)
mesh(th,ph, cost)

hold on 
contour(th,ph,cost,[  3.73:0.02:3.8 3.85 4 4.5 5 6  8   10 15 30 100 200 1000 2000 5000] )
plot(th_opt,ph_opt, '*r','LineWidth',2)
plot(th_opt -pi,ph_opt-pi, '*r','LineWidth',2)
plot(ph_opt,th_opt, '*r','LineWidth',2)
plot(ph_opt -pi,th_opt-pi, '*r','LineWidth',2)
hold off
axis([-pi pi -pi pi 0 10])
xlabel('\theta')
ylabel('\phi')
title('cost function')

