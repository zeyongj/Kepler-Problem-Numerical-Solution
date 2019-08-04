a = 0; % Start time
b = 200; % End time
e = 0.6;

%Initial conditions
q10 = 1-e; 
q20 = 0;
p10 = 0;
p20 = sqrt((1+e)/(1-e));


%Stepsize and mesh
h = 0.0005;

t = a:h:b; % Mesh
N = length(t);
q1 = zeros(N,1); % Values for q1
q2 = zeros(N,1); % Values for q2
p1 = zeros(N,1); % Values for p1
p2 = zeros(N,1); % Values for p2

q1(1) = q10; % Initial values
q2(1) = q20; % Initial values
p1(1) = p10; % Initial values
p2(1) = p20; % Initial values

%Euler steps

for i = 1: N-1
    q1(i+1) = q1(i) + h * p1(i);
    q2(i+1) = q2(i) + h * p2(i);
    q = [q1(i+1),q2(i+1)];
    p1(i+1) = p1(i) - h * q1(i+1)/norm(q)^3;
    p2(i+1) = p2(i) - h * q2(i+1)/norm(q)^3;
    % p1(i+1) = p1(i) - h * ((q1(i+1))/(((q1(i+1))^2+(q2(i+1))^2)^(3/2)));
    % p2(i+1) = p2(i) - h * ((q2(i+1))/(((q1(i+1))^2+(q2(i+1))^2)^(3/2)));
end

A = q1.*p2 - q2.*p1;
H = 0.5*(p1.^2+p2.^2)-1./(q1.^2+q2.^2).^(1/2);

figure(1);
plot(q1,q2)
title('Figure 4: Symplectic Euler''s Method','fontsize',16)
xlabel('q1','fontsize',16)
ylabel('q2','fontsize',16)

figure(2);
plot(t,A);
title('Figure 5: New Augular Momentum','fontsize',16)
xlabel('t','fontsize',16)
ylabel('A','fontsize',16)

figure(3);
plot(t,H);
title('Figure 6: New Hamitonian','fontsize',16)
xlabel('t','fontsize',16)
ylabel('H','fontsize',16)

