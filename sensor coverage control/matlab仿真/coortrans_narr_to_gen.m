function [x,y] = coortrans_narr_to_gen(u,v,sx,sy,delta)
    x = sx + u*delta;
    y = sy + v*delta;