function [u,v] = coortrans_gen_to_narr(x,y,sx,sy,delta)
    u = (x-sx)/delta;
    v = (y-sy)/delta;