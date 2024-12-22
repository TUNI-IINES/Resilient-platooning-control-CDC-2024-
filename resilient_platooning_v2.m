function [Y]=resilient_platooning(t,x,n,beta,Mr,d)
a1= 0.1*((sin( t) + cos(2*t) + sin(3*t)))+1;
a2= 0.4*((cos( t) + sin(4*t) + cos(2*t)))+1;
a3= 1*((-cos( t) - sin(3*t) - sin(1*t)))+1;
delta_s=-.5*[a1;0;a2;a3];
delta_c=1*[a2;0;0;a1];

delta=delta_s+delta_c;
delta1=delta_s+3*delta_c;
delta2=2*delta_c;
Y=Mr*x+[zeros(n+1,1);-(1+beta)*d*ones(n-1,1);zeros(n+1,1);-(1+beta)*d*ones(n-1,1)]+[zeros(n+1,1);delta1;zeros(n+1,1);delta2];