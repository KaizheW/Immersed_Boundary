% Delta Function Y
function w=phi2(r)
w=zeros(4,4);
q=sqrt(1+4*r*(1-r));
w(:,4)=(1+2*r-q)/8;
w(:,3)=(1+2*r+q)/8;
w(:,2)=(3-2*r+q)/8;
w(:,1)=(3-2*r-q)/8;
