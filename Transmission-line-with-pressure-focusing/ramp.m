function out=ramp(inp,ton)
out=inp;
out(1:ton)=inp(1:ton).*sin(pi/2*linspace(0,1,ton)).^2;
out(end-ton+1:end)=inp(end-ton+1:end).*sin(pi/2*linspace(1,0,ton)).^2;