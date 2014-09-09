function a = full(b)
% FULL makes an adiff jacobian full

a = adiff(b.x,full(b.dx),b.root);