function yp=state_space_solver_fast2D(t,y,model)

q=compute_q(model,t);
Ay=model.Ae*y;
CAy=model.Ce*Ay;

y2=zeros(2*model.N,1);
y2(1:model.N)=CAy+q;

sol=model.U \ (model.L \y2(model.pp2));
yp=Ay+model.Be1d*(model.Ee\(model.bigD*sol));