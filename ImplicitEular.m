function [x,t]=ImplicitEular(h)%Parameters Initialization

%Initial parameters
%dh= Step size
%t_max=maximum time
%m=mass
%g=gravitational aceleration
%tolerance=0.000000000001 to stop the iterations.

dh = h;
t_max=2;
m=15;
g=[0,0,-9.81].';
J=[15.234375, 0       ,0;
    0       , 0.46875 ,0;
    0       , 0       ,15.234375];
I=[1 0 0;0 1 0;0 0 1];
tolerance=0.000000000001;

%Initialization of zero arrays for Rotation,Position and velocity vector

t=0:dh:t_max;
t_len=length(t);
%arrays with zeros
omg=zeros(3,t_len);
omg_div_array=zeros(3,t_len);
R=zeros(3,3,t_len);
X=zeros(3,t_len);
x=zeros(t_len,1);
%Initial values stored at index "1"
%values at t=0 ( array index starts at 1)
omg(:,1)=[0, 150, -4.61538].';
X(:,1)=[0,1,0].';
R(:,:,1)=[1,0,0;
          0,1,0;
          0,0,1];
%Start of Implicit Euler loop

for n=2:1:t_len
    
    omg_div_array(:,n)=[0 0 0];
    omg(:,(n))=omg(:,(n-1));
    %Maximum iterations set to 50000. Can be increased.
    for i=1:1:5000
    omgt=omg(:,(n));
    %Calculation of tilda.
    omgtilda=[0 -omgt(3) omgt(2); omgt(3) 0 -omgt(1); -omgt(2) omgt(1) 0] ;
    expdhxomgtilda= I+ ((sin(dh*(norm(omgt))))/norm(omgt))*(omgtilda)+((1-cos(dh*(norm(omgt))))/(norm(omgt)^2))*(omgtilda*omgtilda);
    R(:,:,n)=R(:,:,(n-1))*expdhxomgtilda;
    [Res] = Rescal(J,omg_div_array(:,n),omg(:,(n)),m,g,R(:,:,n),X(:,1));
    if norm(Res)<tolerance
        break 
    end
    %calculation of Third Term in S_t calculation
   [ThirdTerm] = ThirdTermcal(J,omg(:,n));
   Rt=R(:,:,n).';
   %Calculation of Tmatrix
   [Tmatrix] = Tmatrixcal(I,omg(:,1));
   omg(:,n);
   %calculation of last term in S_t equation
   [LTerm] = LTermcal(Rt,m,g,X(:,1),dh,omg(:,n),Tmatrix);
   %S_t calculation
   [St]=(J/dh)+(omgtilda*J)-ThirdTerm-LTerm;
   %Calulation of del.
   Delomg=(-inv(St))*(Res);
   %updating value of derivative and velocity vector.
   omg(:,n)=omg(:,n)+Delomg;
   omg_div_array(:,n)=omg_div_array(:,n)+(Delomg/dh);
 
    end
end
for n=1:1:t_len
    %Calculation ofPosition vector
  X(:,n)=R(:,:,(n))*X(:,(1));
end
% creation of an x vector for graph.
x=X(3,:);
%Function for residual calculation
function [Res] = Rescal(J,omgdiv,omg,m,g,R,X) 

 Jxomgdiv= J*omgdiv;
 JXomg=J*omg;
 omgxJXomg=cross(omg,JXomg);
 mg=m*g;
 Rmg=((R.')*mg);
 XxRmg=cross(X,Rmg);
 Res=Jxomgdiv+omgxJXomg-XxRmg;
end 
%Function for Last term of S_t
function [LTerm] = LTermcal(Rt,m,g,X,t,omg,T)
%Last Term
 Rtxmxg=Rt*m*g;
 Rmgtilda=[0 -Rtxmxg(3) Rtxmxg(2); Rtxmxg(3) 0 -Rtxmxg(1); -Rtxmxg(2) Rtxmxg(1) 0];
 Xt=X(:,1);
 Xtilda=[0 -Xt(3) Xt(2); Xt(3) 0 -Xt(1); -Xt(2) Xt(1) 0];
 hXomg=t*omg;
 LTerm=(t*Xtilda)*(Rmgtilda)*T*hXomg;
 
end

function [ThirdTerm] = ThirdTermcal(J,omg)
%Third Term calculation function
JXomg=J*omg;
ThirdTerm=[0 -JXomg(3) JXomg(2); JXomg(3) 0 -JXomg(1); -JXomg(2) JXomg(1) 0];
 
end
function [Tmatrix] = Tmatrixcal(I,omg)
%Tmatrix function
omgtilda=[0 -omg(3) omg(2); omg(3) 0 -omg(1); -omg(2) omg(1) 0] ;

Tmatrix=I+ ((((cos(norm(omg)))-1)/(norm(omg)^2))*omgtilda)+(1-((sin(norm(omg)))/norm(omg)))*((omgtilda*omgtilda)/(norm(omg))^2);
 
end

end
    



