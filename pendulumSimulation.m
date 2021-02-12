function [x,y,z] = pendulumSimulation(pendulum_length,step_size,end_time,init_cond) 
    g = 9.81;
    fy=@(x,y,z) z;
    fz=@(x,y,z) -3*g/(2*pendulum_length)*sin(y);
    x(1) = 0;
    y(1) = init_cond(1);
    z(1) = init_cond(2);
    h = step_size;
    xfinal=end_time;
    N=ceil((xfinal-x(1))/h);
    for j=1:N
        x(j+1)=x(j)+h;
        k1y=fy(x(j),y(j),z(j));
        k1z=fz(x(j),y(j),z(j));
        k2y=fy(x(j)+h/2,y(j)+h/2*k1y,z(j)+h/2*k1z);
        k2z=fz(x(j)+h/2,y(j)+h/2*k1y,z(j)+h/2*k1z);
        k3y=fy(x(j)+h/2,y(j)+h/2*k2y,z(j)+h/2*k2z);
        k3z=fz(x(j)+h/2,y(j)+h/2*k2y,z(j)+h/2*k2z);
        k4y=fy(x(j)+h,y(j)+h*k3y,z(j)+h*k3z);
        k4z=fz(x(j)+h,y(j)+h*k3y,z(j)+h*k3z);
        y(j+1)=y(j)+h/6*(k1y+2*k2y+2*k3y+k4y);
        z(j+1)=z(j)+h/6*(k1z+2*k2z+2*k3z+k4z);
    end
end