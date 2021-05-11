function D_func = sol_D_star(g,R,nu,v_s,D)

R_ep = ((g*R*D)^0.5)*D/nu;
R_f = exp(-2.891394+0.95296*log(R_ep)-0.056835*(log(R_ep))^2-0.002892*(log(R_ep))^3+0.000245*(log(R_ep))^4);
D_func = v_s-R_f*(g*R*D)^0.5;
