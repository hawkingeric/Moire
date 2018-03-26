% generate hexagonal unitcells by enlagering the lattice constant
% the new lattice vector can be derived from (nx,ny) which is nx*a+ny*b
format long
window_size=20;
n_min = 0;
n_max = 20;%the maximum of nx or ny
m_min = 0;
m_max = 20;
number_of_unit_cells=(n_max+1)*(n_max+2)/2;%number of different sets (nx,ny)
length_1=2.504;%lattice constant of the 1st layer Ag=4.0853(exp), 4.169(PBE), 4.015(LDA) 
length_2=2.556;%lattice constant of the 2nd layer Pb=4.9508(exp), 5.033(PBE), 4.876(LDA)
struc_1=2; %1=hexagonal, 2=honeycomb
struc_2=1; %1=hexagonal, 2=honeycomb
error_criterion=2;
%new_lattice_n=zeros(number_of_unit_cells,3)
%A records the nx^2+ny^2+nx*ny and the (nx,ny) sets in order
A = zeros(number_of_unit_cells,3);
m = 1;
y = zeros(number_of_unit_cells);
for i=n_min:n_max;
    for j=m_min:m_max;
       if j >= i;
           A(m,1)=(i^2+j^2+i*j);
           A(m,2)=i;
           A(m,3)=j;
           m=m+1;
       end;
    end;
end;

new_lattice_n=A;

%{
plot(length_1*sqrt(new_lattice_n(:,1)),y,'bo','MarkerFaceColor','b','MarkerSize',5);
hold on;
plot(length_2*sqrt(new_lattice_n(:,1)),y,'r*','MarkerFaceColor','r','MarkerSize',5);
%}

%B records the two indexes for mathching new lattice
B=zeros(number_of_unit_cells*number_of_unit_cells,2);
%D records the lattice constant change of layer 1
D=zeros(number_of_unit_cells,1);
p=1;
for i=1:number_of_unit_cells;
    for j=1:number_of_unit_cells;
        if abs((length_2*sqrt(new_lattice_n(j,1))/sqrt(new_lattice_n(i,1))/length_1-1)*100)<error_criterion;
           D(p)=(length_2*sqrt(new_lattice_n(j,1))/sqrt(new_lattice_n(i,1))/length_1-1)*100;
           %plot(multiple(i,1),multiple(j,1),'ko','MarkerFaceColor','k','MarkerSize',5)
           %axis equal
           %hold on;
           B(p,1)=i;%index of layer 1 that matches 
           B(p,2)=j;%index of layer 2 that matches
           p=p+1;
        end
    end
end
number_of_match=p-1;
Errors=zeros(number_of_match,1);
match_index=zeros(number_of_match,2);
match_n=zeros(number_of_match,6);
for i=1:number_of_match;
    Errors(i)=D(i);
end;
[Errors,ix2]=sort(Errors);

for j=1:number_of_match;
    match_index(j,1)=B(ix2(j),1);
    match_index(j,2)=B(ix2(j),2);
end;

total_atom_number=0;
disp(['Index ' 'a_1 ' 'b_1 ' 'phi_1 ' 'a_2 ' 'b_2 ' 'phi_2 ' 'n_1 ' 'n_2 ' 'dphi ' 'total_atom_number ' 'lattice_constant ' 'Error '])
for i=1:number_of_match;
    match_n(i,1)=new_lattice_n(match_index(i,1),1);
    match_n(i,2)=new_lattice_n(match_index(i,1),2);
    match_n(i,3)=new_lattice_n(match_index(i,1),3);
    match_n(i,4)=new_lattice_n(match_index(i,2),1);
    match_n(i,5)=new_lattice_n(match_index(i,2),2);
    match_n(i,6)=new_lattice_n(match_index(i,2),3);
    n_1=sym(sqrt(new_lattice_n(match_index(i,1),1)));
    a_1=(new_lattice_n(match_index(i,1),2));
    b_1=sym(new_lattice_n(match_index(i,1),3));
    n_2=sym(sqrt(new_lattice_n(match_index(i,2),1))); 
    a_2=sym(new_lattice_n(match_index(i,2),2));
    b_2=sym(new_lattice_n(match_index(i,2),3));
    switch struc_1;
       case 1;
          total_atom_number=match_n(i,1);
          switch struc_2;
              case 1;
                  total_atom_number=total_atom_number+match_n(i,4);          
              case 2;
                  total_atom_number=total_atom_number+match_n(i,4)*2;
          end;
       case 2;
          total_atom_number=match_n(i,1)*2;
          switch struc_2;
              case 1;
                  total_atom_number=total_atom_number+match_n(i,4);
              case 2;
                  total_atom_number=total_atom_number+match_n(i,4)*2;
          end;
       case 3;
           total_atom_number=match_n(i,1)*2;
    end;
    
    n1 =match_n(i,1);
    n1x=match_n(i,2);
    n1y=match_n(i,3);
    n2 =match_n(i,4);
    n2x=match_n(i,5);
    n2y=match_n(i,6);
    
    if n1x<n1y
       buffer=n1x;
       n1x=n1y;
       n1y=buffer;
    end;
    if n2x<n2y   
       buffer=n2x;
       n2x=n2y;
       n2y=buffer;
    end;
    
    phi_1=sprintf('%0.2f',acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi);
    phi_2=sprintf('%0.2f',acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi);
    dphi =sprintf('%0.2f',abs(acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi-acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi));
    Error=sprintf('%0.5f',Errors(i)); 
    lattice_constant=sprintf('%0.2f',(length_2*sqrt(new_lattice_n(match_index(i,2),1))));
    X=[i,a_1,b_1,phi_1,a_2,b_2,phi_2,n_1,n_2,dphi,total_atom_number,lattice_constant,Error];
   %if abs(acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi-acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi) < 0.1
   %    disp(X)
   %else
   
   %end;
    
    %disp(X)
    
    if length_2*sqrt(new_lattice_n(match_index(i,2),1)) < 30 % &&...
    %   length_2*sqrt(new_lattice_n(match_index(i,2),1)) < 62
       %{
       abs(acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi-acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi) >= 9&&...%dphi
       abs(acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi-acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi) < 11 &&...%dphi
       acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi > 20     &&...                                                  %phi_1
       acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi < 50                                                            %phi_1
       %acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi < 25  &&...                              %phi_2
       %acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi > 22                                %phi_2
                                                               
       %&&abs(acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)))*180/pi-acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)))*180/pi)>20
        %}
        disp(X)
    else
        
    end;
end;
%for file output
%fid = fopen('d.txt','wt'); 
%fprintf(fid,'%g\n', match_n(:,1));
%fprintf(fid,'%g\n', Errors);
%fclose(fid);
match_n;
plot_index=416;
Gen_NewCell(match_n,length_1,length_2,struc_1,struc_2,plot_index);
