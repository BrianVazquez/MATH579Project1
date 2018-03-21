% %good example
% A=4;B=5;C=2;
% k=4;
% a1=1;a2=4;
% b1=1;b2=3;
%---random example
A=8;B=15;C=2;
a1=1; a2=5;
b1=1; b2=3;
k=3;
%---
n=25;
rz=7;
h=0.1;

z=@(x,y) C+sqrt((C^2)*(k-(((x-A)^2)/(A^2))-(((y-B)^2)/(B^2))));
gradphi= @(x,y,z) [ 2*x/(A^2), 2*y/(B^2), 2*z/(C^2)];

x=a1:h:a2; y=b1:h:b2;
lengthx=length(x);
lengthy=length(y);
rx=zeros(lengthx,lengthy);
ry=zeros(lengthx,lengthy);
zvalues=zeros(lengthx,lengthy);
% rzvalues=zeros(lengthx,lengthy);
P=zeros(lengthx,lengthy);
%get z values, reflected domain values
for i=1:lengthx
    for j=1:lengthy
        zvalues(i,j)=z(x(i),y(j));
        gradp=gradphi(x(i),y(j),zvalues(i,j));
        
        t=(rz-zvalues(i,j))/gradp(3);
        rx(i,j)=x(i)+gradp(1)*t;
        ry(i,j)=y(j)+gradp(2)*t; 
        
        tempk=(((rx(i,j)-A)^2)/(A^2))+(((ry(i,j)-B)^2)/(B^2))+(((rz-C)^2)/(C^2));
        P(i,j)=(k-tempk);
    end
end

% [X]=meshgrid(x,y);
% figure
% surf(X,P)

alpha=3; beta=5;
I=@(x,y) sin(pi*x/alpha)*sin(pi*y/beta);

mI=zeros(lengthx,lengthy);
for i=1:lengthx
    for j=1:lengthy
        mI(i,j)=I(x(i),y(j));
    end
end
ef=zeros(lengthx,lengthy);
for i= 2:lengthx-1
    for j=2:lengthy-1
       ef(i,j) = ((mI(i,j+1) + mI(i,j))*(P(i,j+1) - P(i,j))...
           - ((mI(i,j) + mI(i,j-1))*(P(i,j)-P(i,j-1)))...
           + (mI(i+1,j)+mI(i,j))*(P(i+1,j)-P(i,j))...
           - (mI(i,j)+mI(i-1,j)*(P(i,j) - P(i-1,j))))/(2*(h^2));
    end
end

approxP=TIENeumann(mI,ef,h,0);
% hold on
% surf(y,x,approxP)
A1=min(min(rx));
A2=max(max(rx));
B1=min(min(ry));
B2=max(max(ry));
approxsurf=phase2surf(P,15,rz,a1,a2,b1,b2,h,A1,A2,B1,B2);
subplot(2,2,1)
% [gridX,gridY]=meshgrid(x,y);
surf(x',y',((approxsurf)+P+rz*ones(lengthx,lengthy))')
% surf(x',y',approxsurf')
title('Approx')

subplot(2,2,2)
surf(x',y',zvalues')
title('Actual')

subplot(2,2,3)
surf(x',y',(P)')
title('P')


error=norm(approxsurf+P+rz*ones(lengthx,lengthy)-zvalues)/norm(zvalues)
error=norm(approxsurf-zvalues)/norm(zvalues)