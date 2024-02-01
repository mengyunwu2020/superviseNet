
Gamma_int<-function(v_0, v_1, p_2, p, k,Theta_l){

  P_l=list();
  for(i in 1:k) {
    P_l[[i]] = (1.0/(1.0 + v_1/v_0*exp(-abs((Theta_l[[i]]))/v_0 + abs((Theta_l[[i]]))/v_1)*(1-p_2)/p_2));
  }


  return(P_l)

}
