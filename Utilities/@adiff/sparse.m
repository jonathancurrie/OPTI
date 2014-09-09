function a = sparse(b)
% SPARSE makes an adiff jacobian sparse

a = adiff( b.x, sparse(b.dx), b.root);