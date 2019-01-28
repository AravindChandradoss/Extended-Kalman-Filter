clear all; 
close all;
%given data
n=4;
Ap =  0.3;
d = 100;
ti= 0.2;
lamda = 0.9;
Kp = 1088;
Kd = 1525.9;
Tsim=100;
T = 0.01;
T_s = Tsim/T;
%system function
sys = @(pX, sV, sU) [pX(1)+T*( pX(2) -sV); 
                     pX(2)+(T/pX(4))*(-Ap*(pX(2))^2 -d + pX(3));
                     pX(3)+(T/ti)*(-pX(3)+sU);
                     pX(4)] 
%noise covarience
R = [0.0001 0
     0  0.00001]; 
%adding noise
v1 = sqrt(R(1,1)) * randn(1, T_s); 
v2 = sqrt(R(2,2)) * randn(1, T_s); 

x =0*ones(n,Tsim+1);
x(:,1) = [-16.2 18 0 1576]'; 
y=zeros(1,length(Tsim+1));
u=zeros(1,Tsim+1);
x_h =0*ones(n,Tsim+1);
x_h(:,1) = [-14.6 20 0 1500]';
P(:,:,1) = eye(4);

%EKF
for k=1:T_s
 
    %lead vehicle velocity
    if(mod(k,4000)>2000)
        v_l=22;
    elseif(mod(k,4000)<2000)
        v_l=18;
    end
    v_p(k)=v_l;
 z(:,k) =  [x(1,k); x(2,k)]  + [v1(k) ; v2(k)];
    
 %UPDATING
  H(:,:,k) =  [1,T,0,0;
               0,1-( 2*T*Ap*x_h(2,k)/x_h(4,k)),T/x_h(4,k),(-T/x_h(4,k)^2)*(-Ap*x_h(2,k)^2 -d + x_h(3,k))];             
  g(:,:,k) = P(:,:,k) * H(:,:,k)' * (H(:,:,k) * P(:,:,k) * H(:,:,k)' + R)^-1; 
  x_h(:, k) = x_h(:,k) + g(:,:, k) * (z(:,k) - [ x_h(1,k); x_h(2,k)]); 
  P(:,:,k) = (eye(4,4) - g(:,:,k)* H(:,:,k)) * P(:,:,k); 
  y(:,k+1) = -((x_h(1,k) + v1(k)) + lamda *(x_h(2,k) + v2(k)));
  u(:,k+1) = Kp * y(:,k+1) + (Kd/T) * (y(:,k+1) - y(:,k));
 
  F(:,:,k+1) = [1,T,0,0;
                0,1-((2*T*Ap*x_h(2,k))/x_h(4,k)),T/x_h(4,k),(-T/x_h(4,k)^2)*(-Ap*x_h(2,k)^2 -d + x_h(3,k));
                0,0,1-T/ti,0;
                0,0,0,1 ];
     
%PREDICTING
    %Predicted State estimate:
    x_h(1,k+1) = x_h(1,k) + T* (x_h(2,k) - v_l);
    x_h(2,k+1) = x_h(2,k) + (T/x_h(4,k))* (-Ap*(x_h(2,k))^2 - d +x_h(3,k));
    x_h(3,k+1) = x_h(3,k) + (T/ti )* (- x_h(3,k) + u(:,k));
    x_h(4,k+1)= x_h(4,k);
    %x_h(:,k+1)=sys(x_h(:,k),v1(k),u(:,k));    
    P(:,:,k+1) = F(:,:,k+1) * P(:,:,k)* F(:,:,k+1)' ;
              
%Updating State outputs (Actual Model)
    x(1,k+1) = x(1,k) + T* (x(2,k) - v_l);
    x(2,k+1) = x(2,k) + (T/x(4,k))* (-Ap*(x(2,k))^2 - d +x(3,k));
    x(3,k+1) = x(3,k) + (T/ti )* (- x(3,k) + u(:,k));
    x(4,k+1)= x(4,k);
    %x(:,k+1)=sys(x(:,k),v1(k),u(:,k))  
end

time = 0:T_s;
timel = 0:T_s-1;

subplot(5,1,1)
stairs (timel, v_p)
grid on
ylim([16,24])
legend('Velocity (Lead Vehicle','location','best');
xlabel('Time(k)');
ylabel('Velocity(m/s)');
title('Velocity profile of lead vehicle') 


subplot(5,1,2)
plot(time,x(1,:),'b')
hold on
grid on
plot(time,x_h(1,:),'r')
legend('Actual', 'Estimated','location','best');
xlabel('Time(k)');
ylabel('Inter Spacing(m)');
title('InterVehicle spacing X(1)')

subplot(5,1,3)
plot(time,x(2,:),'b')
hold on
grid on
plot(time,x_h(2,:),'r')
legend('Actual', 'Estimated','location','best');
xlabel('Time(k)');
ylabel('Velocity(m/s)');
title('Velocity X(2)')

subplot(5,1,4)
plot(time,x(3,:),'b')
hold on
grid on
plot(time,x_h(3,:),'r')   
legend('Actual', 'Estimated','location','best');
xlabel('Time(k)');
ylabel('Force(N)');
title('Force applied X(3)')

subplot(5,1,5)
stairs(time,x(4,:),'b') 
ylim([1490,1600])
hold on
grid on
stairs(time,x_h(4,:),'r')    
legend('Actual', 'Estimated','location','best');
xlabel('Time(k)');
ylabel('Mass(kg)');                 
title('Mass of the Vehicle')
suptitle('Vehicle Parameter and State Estimation in AHS Using EKF')                   
                   
    
    
    
   