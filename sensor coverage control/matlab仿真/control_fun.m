function [hx,hy] = control_fun(sx,sy,agentposition,delta,V,D)
    R0 = 3.0;
    beita = 0.1;
    p0 = 1.0;
    lamda = 1.0;
    Bx=1;
    position = zeros(1,2);
    rangeposition = zeros(6,3);
    hx = 0;
    hy = 0;
    
    %agentposition = [s1x,s1y;s2x,s2y;s3x,s3y;s4x,s4y;s5x,s5y;s6x,s6y];


    for k = 1:6
        distance = sqrt((agentposition(k,1)-sx)*(agentposition(k,1)-sx)+(agentposition(k,2)-sy)*(agentposition(k,2)-sy));
        if (distance ~= 0) &&(distance<2*D)
            rangeposition(k,1)=1;
            rangeposition(k,2:3)=agentposition(k,:);
        end
    end
    for u = -1*V:V
        for v = -1*V:V
            if u~=0 &&v~=0
                [position(:,1),position(:,2)] = coortrans_narr_to_gen(u,v,sx,sy,delta);
                for m = 1:6
                    if rangeposition(m,1) == 1
                        Bx = Bx*(1-p0*exp(-1*lamda*sqrt((position(:,1)-rangeposition(m,2))*(position(:,1)-rangeposition(m,2))+(position(:,2)-rangeposition(m,3))*(position(:,2)-rangeposition(m,3)))));
                    end
                end
                RX = R0-beita*sqrt(position(:,1)*position(:,1)+(position(:,2)-20)*(position(:,2)-20));
                PX = p0*exp(-1*lamda*sqrt((position(:,1)-sx)*(position(:,1)-sx)+(position(:,2)-sy)*(position(:,2)-sy)));
                dPX = -1*lamda*PX;
                hx = hx+RX*Bx*dPX*u/sqrt(u*u+v*v);
                hy = hy+RX*Bx*dPX*v/sqrt(u*u+v*v);
                Bx = 1;
            end
        end
    end
    hx = hx*delta*delta;
    hy = hy*delta*delta;