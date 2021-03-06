close all

Xt= zeros(1,500);
Yt = zeros(1,500);
prev = 2;
Xt(1) = 0;
%X(:,1) = [0;0];
%Q = 0.64*eye(2);
phip = normrnd(0,0.64);
men = zeros(1,1000);
x = Xt(1);
k = 0.5;
T = 1;
fi = exp(-k*T)
V = prev*fi + phip;
Col = zeros(1,1000);
Yt(1) = x^2/20 +  V;
Col(1) = fi*prev + phip;
men(1) = phip; 
for i = 2:1000
    phip = normrnd(0,0.64);
    chip = normrnd(0,0.64);
    Col(i) = fi*Col(i-1) + phip;
    Xt(i) = 0.5*Xt(i-1) + 25*Xt(i-1)/(1+Xt(i-1)^2) + 8*cos(1.2*(i-1)) + chip;
    Yt(i) = Xt(i)^2  + phip;
    %men(i) = phip;
end




q = 0.64; 
um = 0;
R = 0.64*eye(1);
n = 1;
Q = q*eye(1);
P = Q; %Process Covaraince
inidat;% REMOVE THIS, i am Importing my data from this file
y = Yt; %This is a matrix that contains all the ranges coming from the anchors

Mt = zeros(size(um,1),length(y)); %this matrix contains all the measurements we collect
Con = zeros(size(P,1),size(P,2),length(y));%this matrix keeps track of the covariance we calculate in each step


a = 1; %alpha for lambda
kappa = 0; %for later calculations
betta = 1;
lambda = (a^2)*(n+kappa)-n;
N = length(Yt); %for the number of steps in our Experiment
phil = [sqrt(n)*eye(n) -sqrt(n)*eye(n)]
U = chol(P,'lower'); %cholsky matrix for the algorithm
    X = zeros((2*n)); %for the sigma points
   
     
    for i = 1:2*n %this forms a total of n points + the first point equal to m
    
         X(i) = um + U*phil(:,i);

    end
    
for t = 1:N
    
    for i = 1:2*n
        X(i) = 0.5*(X(i)) +25*X(i)/(1+X(i)^2) + 8*cos((1.2*(t-1)));
    end
    %we have p[assed through the dynamic model of the state equation
    
    
        sum = 0;
    for i = 1 : 2*n
        sum = sum + X(i);
    end
    m = (1/(2*n)) + sum;
    sum = 0;
    
    for i = 1:2*n
        sum = sum + (X(i)-m)*(X(i)-m)';
    end
       P = (1/(2*n))*sum + Q;
       P = real(P);
    L = zeros(size(P));
    L = P - Q;
    L = real(L);
    T = chol(L);
    %so far we have schieved the step 3 requirement of creating a
    %decomposition for the positive definite matrices
    %we also generate the error matrix
    Xerr = zeros(size(X));
    for i = 1:2*n %this forms a total of n points + the first point equal to m
    
         Xerr(i) = X(i) - m;

    end
    %we now also have the error matrix like in the paper
%     Xerr
%     L
%     Temp = Xerr(:,3:4);
%     D = [diag(Xerr) diag(Temp)]
%     Resear = Xerr(:,1:2)*D'
    
%     U = chol(P,'lower');  
%     for i = 1:2*n 
%     
%          X(:,i) = m + U*phil(:,i);
% 
%     end
    %we have our new sigma points
    %importing out data
%     K1 = M(t,:);
%     A1 = K1(1:3); %setting the Co-ordinates of the anchors
%     A2 = K1(5:7);
%     A3 = K1(9:11);
%     A4 = K1(13:15);
%     An=[A1; A2; A3;A4];
%     
%     for i = 1 : size(X,2) % this will equal 2n + 1
%     % we wil pass it through the equation for our non-linear function y and
%     % store the points in a matrix
%     Y(1:4,i)=[sqrt((X(1,i)-An(1,1))^2 + (X(2,i)-An(1,2))^2 + (0.84-An(1,3))^2);...
%               sqrt((X(1,i)-An(2,1))^2 + (X(2,i)-An(2,2))^2 + (0.84-An(2,3))^2);...
%               sqrt((X(1,i)-An(3,1))^2 + (X(2,i)-An(3,2))^2 + (0.84-An(3,3))^2);...
%               sqrt((X(1,i)-An(4,1))^2 + (X(2,i)-An(4,2))^2 + (0.84-An(4,3))^2)] ;
%     %where An are the co-ordinates for our Anchor matrix (this is our non-linear function)
%     % you just need to write anything like y(i) = g(X(i)) for each sigma point and
%     % store it in the matrix Y
% end

for i = 1:size(X,2)
    Y(i) = X(i)/20;
end
       
       sum = 0
     for i = 1:2*n
         sum = sum + Y(i);
     end
    uk = (1/(2*n))*sum;
    sum = 0 ;
    for i = 1:2*n
        sum = sum + (Y(i)-uk)*(Y(i)-uk)';
    end
    Sk = (1/(2*n))*sum + R;
    
    sum = 0;
    for i = 1:2*n
        sum = sum + (X(i)-m)*(Y(i)-uk)';
    end
    Ck = (1/(2*n))*sum;
    
    K = Ck/Sk; %filter gain
    um = m + K*(y(t) - uk); %new mean for the X
    um
    P = P - K*Sk*K';
    P = real(P)
% here um and P are the estimates that the paper refers to in step 5
%     W = zeros(1,2*n);
%     for i = 1:2*n
%         W(i) = 1/(2*n);
%     end
    XUerr = zeros(size(X));
    for i = 1:2*n %this forms a total of n points + the first point equal to m
    
         XUerr(i) = X(i) - um;

    end
    
      %D = diag(XUerr);
%     Temp = XUerr(:,3:4);
%     D = [diag(XUerr) diag(Temp)];
%     
%     J = XUerr; %%%REMEMBER YOU MADE A BIG CHANE HERE AND THIS CAN CAUSE PROBLEMS
%     N = chol(J);
      E = 0.5*eye(size(P));
      N = P ;%- E;
      N = real(N);
      N = chol(N); 
      %Now he says that delta E is something error blah blah
      %Updating the sigma points here
    jhi = N/T;
    jhi;
    XUerr = jhi*Xerr;
    %new sigma points for the next filter cycle
    X = um + XUerr;
    %X
    %t
    um;
    t;
    Mt(t) = um;
    
    Con(:,:,t) = P;
    
    
    
       
end

    plot(Mt)
    hold on
    plot(Xt)
    legend('Coloured Noise Cubature','Actual X');
    title('Algorith 2 estimate with E = 0 and 1000 samples');
    
    
