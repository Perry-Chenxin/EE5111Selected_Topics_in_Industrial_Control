clc; clear all;
% GPS can give us x-direction and y-direction position info
% the number of track points and time interval
N = 120; 
dt = 0.2; 
% noise covariance matrix coefficient
sigma_a1 = 0.5; sigma_a2 = 0.5; 
% sensor noise variance 
sigma_1 = 0.25; sigma_2 = 0.25;
% the variable for cumulative error square
sum_error_square = zeros(6,1);

% initialization
% covariance matrix P and Estimated covariance matrix P_
P_=zeros(6,6,N); 
P=zeros(6,6,N); 
P(:,:,1)=eye(6);
% Kamlan gain
K=zeros(6,2,N);
% state value and its initial value
x=zeros(6,N); x(:,1)=[0 0 3 3 2 2]';
% estimated state value and its initial value
x_=zeros(6,N); x_(:,1)=[1 -1 2.5 3.5 1.5 2.5]';
% measurement value and its initial value
z=zeros(2,N); z(:,1)=[0;0];
% state transition matrix F and measurement matrix G
F=[1 0 dt 0 0.5*dt^2 0;
   0 1 0 dt 0 0.5*dt^2;
   0 0 1 0 dt 0;
   0 0 0 1 0 dt;
   0 0 0 0 1 0;
   0 0 0 0 0 1];

G=[1 0 0 0 0 0;
   0 1 0 0 0 0];

% Covariance matrix of process and measurement noise
R=[sigma_1^2 0;0 sigma_2^2];

T = [dt^2/2 0 dt 0 1 0;0 dt^2/2 0 dt 0 1]';
qc = [sigma_a1^2 0;0 sigma_a2^2];
Q = T*qc*T';
% generating noise
q = random('normal',0,sigma_a1,6,N);
r = random('normal',0,sigma_a2,2,N);
r(:,1) = [0;0];

% KF
for n=2:N    
    % update state and measurement
    x(:,n)=F*x(:,n-1)+q(:,n-1);
    z(:,n)=G*x(:,n)+r(:,n);
    x_(:,n)=F*x_(:,n-1);
    P_(:,:,n)=F*P(:,:,n-1)*F'+Q;
    K(:,:,n)=P_(:,:,n)*G'*inv(G*P_(:,:,n)*G'+R);
    x_(:,n)=x_(:,n)+K(:,:,n)*(z(:,n)-G*x_(:,n));
    P(:,:,n)=P_(:,:,n)-K(:,:,n)*(G*P_(:,:,n)*G'+R)*K(:,:,n)';

    % sum of error square
    sum_error_square = sum_error_square+power(abs(x_(:,n)-x(:,n)),2); 
end

% trajectory figure
figure(1)
hold on
p1=plot(x(1,:),x(2,:),'-dr','linewidth',0.5);%real trace
p2=plot(z(1,:),z(2,:),'*g','linewidth',0.5);%measurement trace
p3=plot(x_(1,:),x_(2,:),'-ob','linewidth',0.5);%filtered trace
hold off
legend([p1,p2,p3],'Real Trace','Measurements Trace','Filtered Trace');
%legend([p1,p3],'real trace','Filtered Trace');

% error figure
figure(2)
p_e = plot(1:N,x_(1,:)-x(1,:),'r',1:N,x_(2,:)-x(2,:),'b');
legend([p_e],'x-position error','y-position error');

% the mean squared error and root mean square error
mse = sum_error_square/N
rmse = sqrt(mse)