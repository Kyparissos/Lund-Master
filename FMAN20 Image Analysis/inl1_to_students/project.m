function r = project(stacks,bases,k,i,j)% r is the norm of the difference
u = stacks{k}(:,:,j);% k is the number of stacks , j is number of images
e_1 = bases{i}(:,:,1);% i is the number of basis , each basis has 4 images
e_2 = bases{i}(:,:,2);
e_3 = bases{i}(:,:,3);
e_4 = bases{i}(:,:,4);
x_1 = sum(dot(u,e_1));% project the image onto a basis so x1 is scalar product
x_2 = sum(dot(u,e_2));% of u and each basis image
x_3 = sum(dot(u,e_3));
x_4 = sum(dot(u,e_4));
u_p = x_1 * e_1 + x_2 * e_2 + x_3 * e_3 + x_4 * e_4;% u_p is projection result
r = sqrt(sum(dot(u-u_p,u-u_p)));% r to the power of 2 =(u-up)·(u-up)
end