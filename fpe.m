% compute pdf of fpe corr to Levy noise with absorbing BC 
% in a bounded, symmetric domain (-r,r);
% central differencing for 2nd-order deriv 
% and using one-sided differencing near the boundaries;
% 3rd-order WENO for 1st-order deriv
% 3rd-order TVD RK in time
% 

%function U=fpe1(alpha,eps,d,r,J)
clear all;
%clf

alpha=0.5;
eps=1;
d=0;
f=@(x) -x;
r=1; 
fpmax=f(-r);


J=50; 
h=1/J;
dt = 0.4*h^2; 
dtdx = dt/h;
T=1;

Jt=2*J-1;     %Jt=total number of unknowns  
x=-2:h:2;
C=-zeta(alpha-1)*h^(2-alpha);          % correction term u''(x)
cons=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2));
% coeff for the diffusion and the correction term
Chh=( d/(2*r^2) + cons*C*eps/r^(alpha) )/h^2;
c1=eps*cons/(alpha*r^alpha);
c2=eps*cons*h/r^alpha;

b=zeros(Jt,1);
a=zeros(Jt,1);
c=zeros(Jt,1); 

%nonintegral part

% coefficient of U_j
b(2:Jt-1) = (-2*Chh - c1*(1./(1+x(J+3:3*J-1)).^alpha+1./(1-x(J+3:3*J-1)).^alpha))';
% one-sided diff near boundaries
b(1)  = 2*Chh - c1*(1/(1+x(J+2))^alpha+1/(1-x(J+2))^alpha); 
b(Jt) = 2*Chh - c1*(1/(1+x(3*J))^alpha+1/(1-x(3*J))^alpha); % one-sided diff
a= Chh*ones(Jt,1);  % coefficient of U_(j-1)
c= Chh*ones(Jt,1);  % coefficient of U_(j+1) 
c(1)  = -5*Chh; % one-sided diff
a(Jt) = -5*Chh; % one-sided diff
vp2 = zeros(Jt,1); vp2(3) = 4*Chh;  % one-sided diff
vp3 = zeros(Jt,1); vp3(4) =  -Chh; % one-sided diff 
vm2 = zeros(Jt,1); vm2(Jt-2) = 4*Chh;  % one-sided diff
vm3 = zeros(Jt,1); vm3(Jt-3) =  -Chh; % one-sided diff 

% integral part
for j=-J+1:J-1
   b(j+J)= b(j+J) - c2*( sum(1./abs(x(J+2-j:2*J)).^(1+alpha)) ...
                       + sum(1./abs(x(2*J+2:3*J-j)).^(1+alpha)) ...
         + .5/abs(x(J+1-j))^(1+alpha) + .5/abs(x(3*J+1-j))^(1+alpha) );  
end 
A=spdiags([vm3 vm2 [a(2:end); 0] b ...
           [0; c(1:end-1)] vp2 vp3],-3:3,Jt,Jt);

% coefficient of u_(j+k) 
B=zeros(size(A));
for j=-J+1:J-1
  B(J+j,:)=[1./abs(x(J+2-j:2*J)).^(1+alpha)  0  1./abs(x(2*J+2:3*J-j)).^(1+alpha)];
end


%A = diag(ones(Jt,1))+ dt*(A+c2*B); % the iterative matrix for each time-step

X=-r:(h*r):r;
%UU=sqrt(pi)*(r^2-X.^2).^(alpha/2)/(2^alpha*gamma(1+alpha/2)*gamma(1/2+alpha/2)); 
UU=sqrt(40/pi)*exp(-40*X.^2);
U=UU(2:end-1)';
Un=U; U1=U; U2=U;
%U=UU(2:end-1);

nft=round(T/dt);
figure
hold on
nu = length(U);
data = zeros(nu+4,1);
for nt=1:nft

    U1 = U + dt*(A+c2*B)*U;
    % global Lax-Friedrichs(LF) flux splitting
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U/(2*r);
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U/(2*r);
    fx2 = derWENOr2_plus(data,h);
    U1 = U1 - dtdx*(fx1+fx2);
    
    U2 = 0.75*U + U1/4 + (dt/4)*(A+c2*B)*U1;
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U1/(2*r);
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U1/(2*r);
    fx2 = derWENOr2_plus(data,h);
    U2 = U2 - (dtdx/4)*(fx1+fx2);
    
    Un = U/3 + 2*U2/3 + (2*dt/3)*(A+c2*B)*U2;
    data(3:nu+2) = (f(X(2:end-1)) + fpmax)'.*U2/(2*r);
    fx1 = derWENOr2_minus(data,h);
    data(3:nu+2) = (f(X(2:end-1)) - fpmax)'.*U2/(2*r);
    fx2 = derWENOr2_plus(data,h);
    Un = Un - (2*dtdx/3)*(fx1+fx2);
    
    if mod(nt,100)==0
       plot(X,[0; U; 0])    
    end
         %pause
    U=Un;
    t=nt*dt;
end
title(['final time t= ', num2str(t)])
  
%U=A\(-1*ones(Jt,1));
%U=[0; U; 0];
%  rescale:
%X=-r:(h*r):r;


%check the order without true solution
%u(n)=U(J+1);
%end
%order=log((u(2)-u(1))/(u(3)-u(2)))/log(2);


%comparison with analytic solution
%UU=sqrt(pi)*(r^2-X.^2).^(alpha/2)/(2^alpha*gamma(1+alpha/2)*gamma(1/2+alpha/2)); 
%figure
%plot(X,U,'b--')
%hold on
%plot(X,UU,'k-')
%legend('numerical solution','analytic solution ')


%check the order:
 %m=find(abs(X-0.1)<1e-6);
 %Er=abs(U(m)-UU(m))     %error= U_j-u(x_j)    %modified
 %step=h
 %[maxErr,im]=max(abs(U-UU'))

%end
 
