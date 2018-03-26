function [new_cell,number_of_points] = Gen_Coord(nx,ny,nn,p,lattice_constant)
digits(10)
a=[1,0];
b=[1/2,sqrt(3)/2];
init_p=p(1)*lattice_constant*a+p(2)*lattice_constant*b;
points=zeros(nn,2);
index_of_point=1;
for i = -ny : nx;
   for j = 0 : (nx + ny*2);
      if((nx*j - ny*i) >= 0 && (nx*j - ny*i) < nn);
          if(((nx + ny)*i + ny*j) >= 0 && ((nx + ny)*i + ny*j) <  nn);
             %Do not shift the lattice
             %{
             points(index_of_point,1)=init_p(1)+i*lattice_constant*a(1)+j*lattice_constant*b(1);%x coordinate
             points(index_of_point,2)=init_p(2)+i*lattice_constant*a(2)+j*lattice_constant*b(2);%y coordinate
             %}
             %Shift the lattice
             
             points(index_of_point,1)=init_p(1)+i*lattice_constant*a(1)+j*lattice_constant*b(1)+0/3*(a(1)+b(1))*lattice_constant;%x coordinate
             points(index_of_point,2)=init_p(2)+i*lattice_constant*a(2)+j*lattice_constant*b(2)+0/3*(a(2)+b(2))*lattice_constant;%y coordinate
             
             
             %plot( lattice_points_1(m,1),lattice_points_1(m,2),'go','MarkerFaceColor','g','MarkerSize',5);
             %axis([0 50 -50 50])
             %axis equal
             %hold on ;
             index_of_point=index_of_point+1;
         end;
      end;
   end;
end;
number_of_points=index_of_point-1;
new_cell=zeros(number_of_points,1);
for i = 1 : number_of_points;
    new_cell(i,1)=points(i,1);
    new_cell(i,2)=points(i,2);
end;
