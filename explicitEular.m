
function [x,t] = explicitEular(h)
%parameters
%Initial parameters
%dh= Step size
%t_max=maximum time
%m=mass
%g=gravitational aceleration
dh = h;
t_max=2;
m=15;
g=[0,0,-9.81].';
J=[15.234375, 0       ,0;
    0       , 0.46875 ,0;
    0       , 0       ,15.234375];
Jin= inv(J);
I=[1 0 0;0 1 0;0 0 1];
t=0:dh:t_max;
t_len=length(t);
%Initialization of zero arrays for Rotation,Position and velocity vector
%arrays with zeros
omg=zeros(3,t_len);
omg_div_array=zeros(3,t_len);
R=zeros(3,3,t_len);
X=zeros(3,t_len);
x=zeros(t_len,1);


%values at t=0 ( array index starts at 1)
omg(:,1)=[0, 150, -4.61538].';
X(:,1)=[0,1,0].';

R(:,:,1)=[1,0,0;
          0,1,0; 
          0,0,1];
 %Start of Explicit Euler loop     
for n=2:1:t_len
   %calculation of derivative
  omg_div_array(:,n-1) =  diff(J,R(:,:,(n-1)),m,g,omg(:,(n-1)),X(:,1),Jin);
  %Updating Omega
  omg(:,n)=omg(:,(n-1))+(dh*omg_div_array(:,(n-1)));
  omgt=omg(:,(n-1));
  omgtilda=[0 -omgt(3) omgt(2); omgt(3) 0 -omgt(1); -omgt(2) omgt(1) 0] ;
  expdhxomgtilda= I+((sin(dh*(norm(omgt))))/norm(omgt))*(omgtilda)+((1-cos(dh*(norm(omgt))))/(norm(omgt)^2))*(omgtilda*omgtilda);
 %Updating Rotation matrix
  R(:,:,n)=R(:,:,(n-1))* expdhxomgtilda;
  
end

%Calculating X vector. 
for n=1:1:t_len
  X(:,n)=R(:,:,(n))*X(:,(1));
  x(n,1)=X(3,n);
end
%Function for derivative calculation.
 function [omg_div] = diff(J,R,m,g,omg,X,Jin) 

 Jxomg= J*omg;
 omgxJxomg=cross(omg,Jxomg);
 mg=m*g;
 Rmg=((R.')*mg);
 XxRmg=cross(X,Rmg);
 A=XxRmg-omgxJxomg;
 omg_div=(Jin*A);

 end

end
 