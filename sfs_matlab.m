rgb_matrix = imread('target.png');
figure, imshow(rgb_matrix);
gray_matrix = rgb2gray(rgb_matrix);
% Gauss blur
stru = fspecial('gaussian',[5 5], 2);
gray_matrix_1 = rgb2gray(imfilter(rgb_matrix, stru, 'same'));
[h1 h2] = size(gray_matrix_1);
disp(h1)
disp(h2)
figure, imshow(gray_matrix_1);
rect_cons = 0.9;
% light source direction
light_source_xy1 = pi/6;
light_source_xy2 = pi/3;
% source matrix
R=[cos(light_source_xy1)*cos(light_source_xy2) cos(light_source_xy1)*sin(light_source_xy2) -sin(light_source_xy1);...
    -sin(light_source_xy2) cos(light_source_xy2) 0;...
    sin(light_source_xy1)*cos(light_source_xy2)  sin(light_source_xy1)*sin(light_source_xy2) cos(light_source_xy1)];

[m n] = size(gray_matrix);
gray_matrix_source = ones(m+3, n+3);
gray_matrix_source(4:m+3, 4:n+3) = gray_matrix_1(:, :);
% process edge to 
for i=1:3,
    gray_matrix_source(i, 4:n+3) = gray_matrix_1(4-i,:);
end

for i=1:3,
    gray_matrix_source(:,i) = gray_matrix_source(:, 7-i);
end
% gray_matrix_source along light.

%
diff_x = ones(m + 3, n + 3);
diff_y = ones(m + 3, n + 3);
E_x = ones(m + 3, n + 3);
E_y = ones(m + 3, n + 3);
% first difference 
diff_x(2 : m + 3, :) = diff(gray_matrix_source, 1, 1);
diff_y(:, 2 : n + 3) = diff(gray_matrix_source, 1, 2);

for i = 1:m+3,
    for j = 1:n+3,
        E_x(i, j) = diff_x(i, j)*(cos(light_source_xy1)*cos(light_source_xy2)) +  diff_y(i, j)*(cos(light_source_xy1)*sin(light_source_xy2));
        E_y(i, j) = diff_x(i, j)*(-sin(light_source_xy2)) +  diff_y(i, j)*cos(light_source_xy2);
    end
end

% calculate dip angle and drift angle
dip_angle = ones(m+3, n+3);
drift_angle = ones(m+3, n+3);
max_value = max(max(gray_matrix_source));
for i = 1:m+3,
    for j = 1:n+3,
        dip_angle(i, j)= asin(gray_matrix_source(i, j) / max_value);
        drift_angle(i, j) = atan(E_y(i, j)/E_x(i, j));
    end
end

% normal = ones(m + 3, n + 3);
% for i=1:3,
%     normal_3(:,:,i)=normal;
% end
% 
% % calculate normal 
% for i = 1:m+3,
%     for j = 1:n+3,
%     normal(i,j,1) = sin(dip_angle(i, j)) * cos(drift_angle(i, j));
%     normal(i,j,2) = sin(dip_angle(i, j)) * sin(drift_angle(i, j));
%     normal(i,j,3) = cos(drift_angle(i, j));
%     end
% end

% height along source direction
height = ones(m + 3, n + 3);
for i=1:3,
    height_3(:,:,i)=height;
end
disp('222222222222222222222222')
disp(height_3(180,712,:))

for i = 1:m+3,
    for j = 1:n+3,
    height(i,j,1) = sin(dip_angle(i, j)) * cos(drift_angle(i, j))*rect_cons;
    height(i,j,2) = sin(dip_angle(i, j)) * sin(drift_angle(i, j))*rect_cons;
    height(i,j,3) = cos(drift_angle(i, j))*rect_cons;
    end
end

R_ivs_matrix = [-cos(light_source_xy1)*cos(light_source_xy2) sin(light_source_xy1) sin(light_source_xy1)*cos(light_source_xy2);...
    cos(light_source_xy1)*sin(light_source_xy2) cos(light_source_xy2) -sin(light_source_xy1)*sin(light_source_xy2);...
    sin(light_source_xy1) 0 cos(light_source_xy1)];

% real_h = ones(m+3, n+3);
% for i=1:m+3,
%     for j = 1:n+3,
%         temp = [height(i, j ,1) height(i, j ,2) height(i, j ,3)]';
%         real_h(i, j) = sum((R_ivs_matrix*temp).^2);
%     end
% end

H = ones(m, n);
for i = 4:m +3,
    H(i-3,:) = height(i,4:n+3);
end

for i =2:m-1,
    for j = 2:n-1,
        if H(i, j) == NaN,
            H(i, j) = (H(i,j-1) + H(i-1,j) + H(i,j+1) + H(i+1,j))/4 +(H(i-1,j-1)+ H(i-1,j+1)+H(i+1,j-1)+H(i+1,j+1))/4;
        end
    end
end

for i = 2 : m-1,
    if H(i, 1) == NaN,
        H(i, 1) = (H(i-1,1)+H(i+1,1)+ H(i-1,2)+H(i,2)+H(i+1,2))/5;
    end
end

for i = 2 : m-1,
    if H(i, n) == NaN,
        H(i, n) = (H(i-1,n)+H(i+1,n)+ H(i-1,n-1)+H(i,n-1)+H(i+1,n-1))/5;
    end
end

for j = 2 : n-1,
    if H(1, j) == NaN,
        H(1, j) = (H(1, j - 1)+H(1, j + 1)+ H(2, j)+H(2, j-1)+H(2, j+1))/5;
    end
end

for j = 2 : n-1,
    if H(m, j) == NaN,
        H(m, j) = (H(m, j - 1)+H(m, j + 1)+ H(m-1, j)+H(m-1, j-1)+H(m-1, j+1))/5;
    end
end

H(1, 1) = (H(2, 2) + H(1,2) + H(2, 1))/3;
H(m, n) = (H(m-1, n-1) +H(m, n-1) + H(m-1, n))/3;
H(m, 1) = (H(m-1, 1) + H(m-1, 2) + H(m, 2))/3;
H(1, n) = (H(1, n-1) + H(2, n-1) + H(2, n))/3;

x = 0.1:0.1:n/10;
y = 0.1:0.1:m/10;
z = H(uint16(10*y), uint16(10*x))
% H_last = mat2gray(H);
% draw height pic
% figure,
% h = surf(z, x, y);
% shading flat
figure, mesh(x, y, z);
