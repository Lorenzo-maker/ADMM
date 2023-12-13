function dt_J = dt_floatingBaseBodyJacZYX(J,dt_q)

dt_E = dt_q(4:6);

u_x = dt_E(3);
u_y = J(5,5)*dt_E(2);
u_z = J(6,5)*dt_E(2);
w_x = J(1,3)*dt_E(1)+u_x;
w_y = J(2,3)*dt_E(1)+u_y;
w_z = J(3,3)*dt_E(1)+u_z;

dt_J = [ J(2,1)*w_z-J(3,1)*w_y, J(2,2)*w_z-J(3,2)*w_y, J(2,3)*w_z-J(3,3)*w_y, 0, 0, 0;
         J(3,1)*w_x-J(1,1)*w_z, J(3,2)*w_x-J(1,2)*w_z, J(3,3)*w_x-J(1,3)*w_z, 0, 0, 0;
         J(1,1)*w_y-J(2,1)*w_x, J(1,2)*w_y-J(2,2)*w_x, J(1,3)*w_y-J(2,3)*w_x, 0, 0, 0;
         0, 0, 0, J(2,3)*u_z-J(3,3)*u_y,          0, 0;
         0, 0, 0, J(3,3)*u_x-J(1,3)*u_z, J(6,5)*u_x, 0;
         0, 0, 0, J(1,3)*u_y-J(2,3)*u_x,-J(5,5)*u_x, 0];

end