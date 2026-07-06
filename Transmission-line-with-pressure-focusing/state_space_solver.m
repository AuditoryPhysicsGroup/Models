function yp=state_space_solver(t,y,model)


 
yp=((model.A*y)+model.B*compute_q(model,t));
