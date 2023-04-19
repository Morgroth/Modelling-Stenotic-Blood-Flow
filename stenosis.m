clear all
close all
clc

% CONSTANTS

beta = 1.29;
R_0 = 0.01;
R = 0.15;
R_p = 0.003; 
n = 2;
tao_w = 0.02; 
P = 0.5;
delta = 0.2;
num_v = 10;
num_u = 4;
r = linspace(0,R,num_v);
u = linspace(0.2,0.8,num_u);
markers = ['v','s','X','o'];




% Plug Flow Velocity

v_p_values = ones(size(num_v*num_u));

for vis = 1:num_u
    for i = 1:num_v
        display(r(i))
        v_p_values(num_v*(vis-1)+i) = velocityPlug(n,P,u(vis),r(i),beta);
    end    
end    


figure(1)

for vis = 1:num_u
    plot(r(:),v_p_values(num_v*(vis-1)+1:num_v*vis),strcat('-',markers(vis),'k'));
    legend()
    hold on
end

legend('u = 0.2', 'u = 0.4', 'u = 0.6', 'u = 0.8')
title('Variation of Plug Velocity with Obstructed Radius')
xlabel('Obstructed Radius') 
ylabel('Plug velocity') 


% Flow Velocity Profile

x = linspace(R_0,R,num_v);
 
v_values = ones(size(num_v*num_u));
 
for vis = 1:num_u
    for i = 1:num_v
        display(r(i))
        v_values(num_v*(vis-1)+i) = velocityFlow(n,P,u(vis),R,beta,x(i));
    end    
end    
 

 
figure(2)
 
for vis = 1:num_u
     plot(r(:),v_values(num_v*(vis-1)+1:num_v*vis),strcat('-',markers(vis),'k'));
     legend()
     hold on
end
 
legend('u = 0.2', 'u = 0.4', 'u = 0.6', 'u = 0.8')
title('Variation of Flow Velocity in the Blood Vessel at R = 0.15')
xlabel('Radial Distance') 
ylabel('Flow velocity') 


% Flow Rate

q_values = ones(size(num_v*num_u));

for vis = 1:num_u
    for i = 1:num_v
        display(r(i))
        q_values(num_v*(vis-1)+i) = flowRate(n,P,u(vis),r(i),beta);
    end    
end    


figure(3)

for vis = 1:num_u
    plot(r(:),q_values(num_v*(vis-1)+1:num_v*vis),strcat('-',markers(vis),'k'));
    legend()
    hold on
end

legend('u = 0.2', 'u = 0.4', 'u = 0.6', 'u = 0.8')
title('Variation of Volumetric Flow Rate with Obstructed Radius')
xlabel('Obstructed Radius') 
ylabel('Volumetrix Flow Rate') 


% Apparent Fluidity

phi_values = ones(size(num_v*num_u));

for vis = 1:num_u
    for i = 1:num_v
        display(r(i))
        phi_values(num_v*(vis-1)+i) = apparentFluidity(n,delta,u(vis),r(i),beta);
    end    
end    


figure(4)

for vis = 1:num_u
    plot(r(:),phi_values(num_v*(vis-1)+1:num_v*vis),strcat('-',markers(vis),'k'));
    legend()
    hold on
end

legend('u = 0.2', 'u = 0.4', 'u = 0.6', 'u = 0.8')
title('Variation of Apparent Fluidity with Obstructed Radius')
xlabel('Obstructed Radius') 
ylabel('Apparent Fluidity') 



function [v_p] = velocityPlug(n,P,u,r,beta)

a = n/(n+1);
b = (P/(2*u))^(1/n);
c = r^(1+1/n);
d = (1+beta)^(1+1/n);
e = (2*beta)^(1+1/n);

v_p = a*b*c*(e-d);

display(v_p)
display(r)

end

function [v] = velocityFlow(n,P,u,r,beta,x)

a = n/(n+1);
b = (P/(2*u))^(1/n);
c = r^(1+1/n);
d = (1+beta)^(1+1/n);
e = (x/r+beta)^(1+1/n);

v = a*b*c*(d-e);

display(v)
display(r)

end



function [q] = flowRate(n,P,u,r,beta)

a = n*pi/(3*n+1);
b = (P/(2*u))^(1/n);
c = r^(3*(1+1/n));
d = 1;
e = (3*n+1)*beta/(n*(2*n+1));

q = a*b*c*(d-e);

display(q)
display(r)

end


function [phi] = apparentFluidity(n,delta,u,r,beta)

a = (1/u)^2;
b = (-1*r)^(1+1/n);
c = 1-3*n*delta/b;
d = 1;
e = (3*n+1)*beta/(2*n+1);

phi = a*c*(e-d);

display(phi)
display(r)

end





