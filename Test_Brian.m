k = 4;
sensor_z = 3;

z = @(x,y) sqrt(k - (((x-2).^2) + ((y-2).^2))) + 2;
I= @(x,y) sin(pi*x/3)*sin(pi*y/4)+1;
h = 0.05;

% Grid where gradient of Phi crosses positive z values
gradX = 1:h:3;
gradY = 1.5:h:2.5;

% Grid for when gradient of Phi crosses negative z values


length_gradX = length(gradX);
length_gradY = length(gradY);

[gX, gY] = meshgrid(gradX, gradY);

ze = zeros(length_gradX,length_gradY);
ex = zeros(1,length_gradX);
wy = zeros(1,length_gradY);
P = zeros(length_gradX,length_gradY);

for i = 1:length_gradX
    for j = 1:length_gradY
        temp_x = gradX(i);
        temp_y = gradY(j);
        temp_z = z(temp_x,temp_y);
        ze(i,j) = temp_z;
        t = (-temp_z + sensor_z)/(2*temp_z - 4);
        ex(i) = temp_x + (2*temp_x - 4)*t;
        wy(j) = temp_y + (2*temp_y - 4)*t;
        temp_k = ((ex(i) - 2)^2) + ((wy(j) - 2)^2) + ((sensor_z - 2)^2);
        P(i,j) = k - temp_k;
    end
end




mI = zeros(length_gradX,length_gradY);
ef = zeros(length_gradX,length_gradY);

for i = 1:length_gradX
    for j = 1:length_gradY
        mI(i,j) = I(ex(i),wy(j));
    end
end

for i= 2:length_gradX-1
    for j=2:length_gradY-1
       ef(i,j) = ((mI(i,j) + mI(i,j))*(P(i,j+1) - P(i,j))...
           - ((mI(i,j) + mI(i,j-1))*(P(i,j)-P(i,j-1)))...
           + (mI(i+1,j)+mI(i,j))*(P(i+1,j)-P(i,j))...
           - (mI(i,j)+mI(i-1,j)*(P(i,j) - P(i-1,j))))/(2*(h^2));
    end
end
% 
% [pp, lambda] = TIENeumann(mI,ef,h,0);




%size(gradP)