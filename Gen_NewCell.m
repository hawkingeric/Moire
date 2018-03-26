function output = Gen_NewCell(match_n,length_1,length_2,struc_1,struc_2,plot_index)
    output=0;
    n1 =match_n(plot_index,1);
    n1x=match_n(plot_index,2);
    n1y=match_n(plot_index,3);
    n2 =match_n(plot_index,4);
    n2x=match_n(plot_index,5);
    n2y=match_n(plot_index,6);
    
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
   
    
    p_1=[0,0];
    q_1=[1/3,1/3];
    %for WSe2 on Al2O3 c-plane
    %============================
    %r_1=[2/3,1/3];
    %s_1=[0.0,2/3];
    %t_1=[1/3,0.0];
    r_1=[0.63983331,1/3];
    s_1=[0.02683337,0.639833312];
    t_1=[1/3,0.02683337];
    %============================
    p_2=[0,0];
    q_2=[1/3,1/3];
    
    
    %plot the new unitcell based on the exact lattice constant of 1st layer 
    axis equal
    plot([0, length_1*sqrt(n1)],[0,0],'-k','LineWidth',2);%bottom line
    hold on;
    plot([0,length_1*sqrt(n1)*cos(pi/3)],[0,length_1*sqrt(n1)*sin(pi/3)],'-k','LineWidth',2);%left line
    hold on;
    plot([length_1*sqrt(n1),length_1*sqrt(n1)+length_1*sqrt(n1)*cos(pi/3)],[0,length_1*sqrt(n1)*sin(pi/3)],'-k','LineWidth',2);%right line
    hold on
    plot([length_1*sqrt(n1)*cos(pi/3),length_1*sqrt(n1)+length_1*sqrt(n1)*cos(pi/3)],[length_1*sqrt(n1)*sin(pi/3),length_1*sqrt(n1)*sin(pi/3)],'-k','LineWidth',2);%top line
    hold on;
    
     
    phi_1=acos((n1x^2+n1-n1y^2)/(2*n1x*sqrt(n1)));
    phi_2=acos((n2x^2+n2-n2y^2)/(2*n2x*sqrt(n2)));
    rotation_phi_1=[cos(phi_1),sin(phi_1);-sin(phi_1),cos(phi_1)];
    rotation_phi_2=[cos(phi_2),sin(phi_2);-sin(phi_2),cos(phi_2)];
    
    
   
    
    %express the coordinates on the basis of (1,0)L and (0.5,sqrt(3)/2)L
    %where L= length_2*sqrt(n2)
    basis = [length_2*sqrt(n2),1/2*length_2*sqrt(n2);0,sqrt(3)/2*length_2*sqrt(n2)];

    switch struc_1;
        case 1;
            %[new_cell_p_1,number_of_points]=Gen_Coord(n1x,n1y,n1,p_1,length_1);    
            [new_cell_p_1,number_of_points]=Gen_Coord(n1x,n1y,n1,p_1,length_2*sqrt(n2)/sqrt(n1)); 
            new_cell_p_1=rotation_phi_1*new_cell_p_1';
            new_cell_p_1=new_cell_p_1';    
            plot(new_cell_p_1(:,1),new_cell_p_1(:,2),'go','MarkerFaceColor','g','MarkerSize',10);
            axis([0 50 -50 50]);
            axis equal;
            hold on;  
            disp('p_1(btm p) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_p_1')') %first layer basis 1
        case 2;
            %[new_cell_p_1,number_of_points]=Gen_Coord(n1x,n1y,n1,p_1,length_1);    
            [new_cell_p_1,number_of_points]=Gen_Coord(n1x,n1y,n1,p_1,length_2*sqrt(n2)/sqrt(n1));
            new_cell_p_1=rotation_phi_1*new_cell_p_1';
            new_cell_p_1=new_cell_p_1';    
            plot(new_cell_p_1(:,1),new_cell_p_1(:,2),'go','MarkerFaceColor','g','MarkerSize',10);
            axis([0 50 -50 50]);
            axis equal;
            hold on;  
            %[new_cell_q_1,number_of_points]=Gen_Coord(n1x,n1y,n1,q_1,length_1);
            [new_cell_q_1,number_of_points]=Gen_Coord(n1x,n1y,n1,q_1,length_2*sqrt(n2)/sqrt(n1));
            new_cell_q_1=rotation_phi_1*new_cell_q_1';
            new_cell_q_1=new_cell_q_1';
            plot(new_cell_q_1(:,1),new_cell_q_1(:,2),'go','MarkerFaceColor','g','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            disp('p_1(btm p) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_p_1')') %first layer basis 1
            disp('q_1(btm q) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_q_1')') %first layer basis 2
        case 3;
            %[new_cell_p_1,number_of_points]=Gen_Coord(n1x,n1y,n1,p_1,length_1);
            [new_cell_p_1,number_of_points]=Gen_Coord(n1x,n1y,n1,p_1,length_2*sqrt(n2)/sqrt(n1));
            new_cell_p_1=rotation_phi_1*new_cell_p_1';
            new_cell_p_1=new_cell_p_1';    
            plot(new_cell_p_1(:,1),new_cell_p_1(:,2),'go','MarkerFaceColor','g','MarkerSize',10);
            axis([0 50 -50 50]);
            axis equal;
            hold on;  
            %[new_cell_q_1,number_of_points]=Gen_Coord(n1x,n1y,n1,q_1,length_1);
            [new_cell_q_1,number_of_points]=Gen_Coord(n1x,n1y,n1,q_1,length_2*sqrt(n2)/sqrt(n1));
            new_cell_q_1=rotation_phi_1*new_cell_q_1';
            new_cell_q_1=new_cell_q_1';
            plot(new_cell_q_1(:,1),new_cell_q_1(:,2),'go','MarkerFaceColor','g','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            %[new_cell_r_1,number_of_points]=Gen_Coord(n1x,n1y,n1,r_1,length_1);
            [new_cell_r_1,number_of_points]=Gen_Coord(n1x,n1y,n1,r_1,length_2*sqrt(n2)/sqrt(n1));
            new_cell_r_1=rotation_phi_1*new_cell_r_1';
            new_cell_r_1=new_cell_r_1';
            plot(new_cell_r_1(:,1),new_cell_r_1(:,2),'go','MarkerFaceColor','g','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            %[new_cell_s_1,number_of_points]=Gen_Coord(n1x,n1y,n1,s_1,length_1);
            [new_cell_s_1,number_of_points]=Gen_Coord(n1x,n1y,n1,s_1,length_2*sqrt(n2)/sqrt(n1));
            new_cell_s_1=rotation_phi_1*new_cell_s_1';
            new_cell_s_1=new_cell_s_1';
            plot(new_cell_s_1(:,1),new_cell_s_1(:,2),'go','MarkerFaceColor','g','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            %[new_cell_t_1,number_of_points]=Gen_Coord(n1x,n1y,n1,t_1,length_1);
            [new_cell_t_1,number_of_points]=Gen_Coord(n1x,n1y,n1,t_1,length_2*sqrt(n2)/sqrt(n1));
            new_cell_t_1=rotation_phi_1*new_cell_t_1';
            new_cell_t_1=new_cell_t_1';
            plot(new_cell_t_1(:,1),new_cell_t_1(:,2),'go','MarkerFaceColor','g','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            
            
            
            disp('p_1(btm p) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_p_1')') %first layer basis 1
            disp('q_1(btm q) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_q_1')') %first layer basis 2
            disp('r_1(btm r) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_r_1')') %first layer basis 3
            disp('s_1(btm s) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_s_1')') %first layer basis 4
            disp('t_1(btm t) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_t_1')') %first layer basis 5
            
    end;
    switch struc_2;
        case 1;
            [new_cell_p_2,number_of_points]=Gen_Coord(n2x,n2y,n2,p_2,length_2);   
            new_cell_p_2=rotation_phi_2*new_cell_p_2';
            new_cell_p_2=new_cell_p_2';
            plot(new_cell_p_2(:,1),new_cell_p_2(:,2),'ko','MarkerFaceColor','k','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            disp('p_2(top p) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_p_2')') %second layer basis 1
        case 2;
            [new_cell_p_2,number_of_points]=Gen_Coord(n2x,n2y,n2,p_2,length_2);   
            new_cell_p_2=rotation_phi_2*new_cell_p_2';
            new_cell_p_2=new_cell_p_2';
            plot(new_cell_p_2(:,1),new_cell_p_2(:,2),'ro','MarkerFaceColor','r','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            [new_cell_q_2,number_of_points]=Gen_Coord(n2x,n2y,n2,q_2,length_2);
            new_cell_q_2=rotation_phi_2*new_cell_q_2';
            new_cell_q_2=new_cell_q_2';
            plot(new_cell_q_2(:,1),new_cell_q_2(:,2),'ro','MarkerFaceColor','r','MarkerSize',5);
            axis([0 50 -50 50]);
            axis equal;
            hold on;
            disp('p_2(top p) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_p_2')') %second layer basis 1
            disp('q_2(top q) coordinates on the hexagonal basis')
            disp((inv(basis)*new_cell_q_2')') %second layer basis 2
    end;
   %new_cell_p_1
   %new_cell_q_1
   %new_cell_p_2
   %new_cell_q_2
   disp(['Plotting index: ' num2str(plot_index)])
   disp(['Green circles : 1st(btm) layer(a=' num2str(length_1) ')'])
   disp(['Black circles : 2nd(top) layer(a=' num2str(length_2) ')'])






