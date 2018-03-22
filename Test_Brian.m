k = 4;
rz = 5;

z = @(x,y) sqrt(k - (((x-2).^2) + ((y-2).^2))) + 2;
gradphi = @(x,y,z) [2*x,2*y,2*z];
h = 0.1;

% Grid D where gradient of Phi crosses positive z values
senseX = 0:h:2;
senseY = 0:h:2;
lengthx=length(senseX);
lengthy=length(senseY);
rx=zeros(lengthx,lengthy);
ry=zeros(lengthx,lengthy);
% length_gradX = length(senseX);
% length_gradY = length(senseY);
for i=1:lengthx
    for j=1:lengthy
        zvalues(i,j)=z(x(i),y(j));
        gradp=gradphi(x(i),y(j),zvalues(i,j));
        
        t=(rz-zvalues(i,j))/gradp(3);
        rx(i,j)=x(i)+gradp(1)*t;
        ry(i,j)=y(j)+gradp(2)*t; 
        
        tempk=(((rx(i,j)-A)^2)/(A^2))+(((ry(i,j)-B)^2)...
            /(B^2))+(((rz-C)^2)/(C^2));
        P(i,j)=(k-tempk);
    end
end

alpha=3; beta=5;
I=@(x,y) sin(pi*x/alpha)*sin(pi*y/alpha);

mI=zeros(lengthx,lengthy);
for i=1:lengthx
    for j=1:lengthy
        mI(i,j)=I(senseX(i),senseY(j));
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

% approxP=TIENeumann(mI,ef,h,0);
% % hold on
% % surf(y,x,approxP)
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

% Grid A for when gradient of Phi crosses negative z values. In this case,
% the area is smaller than grid D. 
% Dx = zeros(length_gradX,1);
% Dy = zeros(length_gradY,1);
% Dx(1) = 1.3;
% Dy(1) = 1.6;
% Dx(end) = 2.7;
% Dy(end) = 2.2;
% ze = zeros(length_gradX,length_gradY);
% %Bottom/top grids
% for i = 2:length(Dx) - 1
%     % for top grids  -------------------
%     a = gradX(i);
%     b = gradY(1);
%     c = z(a,b);
%     t = (b - Dy(1))/(2*b);
%     Dx(i) = a - a*t;
%     temp_z0 = c - 2*c*t;
%     ze(i,1) = temp_z0;
%     
%     %bottom grid ------------------------
%     
% end

% for j = 1:length(Dy)
%     a1 = senseX(1);
%     b1 = senseY(j);
%     c1 = z(a1,b1);
%     t1 = (a1 - Dx(1))/(2*a1);
%     Dy(j) = a1 - a1*t1;
%     temp_z10 = c1 - 2*c1*t1;
%     ze(1,j) = temp_z10;
% end


% [gX, gY] = meshgrid(senseX, senseY);
% 
% ze = zeros(length_gradX,length_gradY);
% ex = zeros(1,length_gradX);
% wy = zeros(1,length_gradY);
% P = zeros(length_gradX,length_gradY);
% 
% for i = 1:length_gradX
%     for j = 1:length_gradY
%         temp_x = senseX(i);
%         temp_y = senseY(j);
%         temp_z = z(temp_x,temp_y);
%         ze(i,j) = temp_z;
%         t = (-temp_z + sensor_z)/(2*temp_z - 4);
%         ex(i) = temp_x + (2*temp_x - 4)*t;
%         wy(j) = temp_y + (2*temp_y - 4)*t;
%         temp_k = ((ex(i) - 2)^2) + ((wy(j) - 2)^2) + ((sensor_z - 2)^2);
%         P(i,j) = k - temp_k;
% %     end
% end
% 
% ze = ze + P + sensor_z;
% surf(gX,gY,ze');
% 
% 
% mI = zeros(length_gradX,length_gradY);
% ef = zeros(length_gradX,length_gradY);
% 
% for i = 1:length_gradX
%     for j = 1:length_gradY
%         mI(i,j) = I(ex(i),wy(j));
%     end
% end
% 
% for i= 2:length_gradX-1
%     for j=2:length_gradY-1
%        ef(i,j) = ((mI(i,j+1) + mI(i,j))*(P(i,j+1) - P(i,j))...
%            - ((mI(i,j) + mI(i,j-1))*(P(i,j)-P(i,j-1)))...
%            + (mI(i+1,j)+mI(i,j))*(P(i+1,j)-P(i,j))...
%            - (mI(i,j)+mI(i-1,j)*(P(i,j) - P(i-1,j))))/(2*(h^2));
%     end
% end
% 
% figure
% surf(gX',gY',P)
% 
% approxP=TIENeumann(mI,ef,h,0);
% figure
% surf(gX',gY',approxP)

% Calculating Phi(x,y, z(x,y)):




