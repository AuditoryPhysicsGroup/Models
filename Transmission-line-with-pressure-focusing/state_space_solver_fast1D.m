function yp=state_space_solver_fast1D(t,y,model)

 
q=compute_q(model,t);
Ay=model.Ae*y;
CAy=model.Ce*Ay;
yp=Ay+model.Be*((model.Fcb\(CAy+q)));

