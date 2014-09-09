function [id]=dmumps(id,mat)
%
% [id]=dmumps(id,mat)
% id is a structure (see details in initmumps.m and MUMPS documentation)
% mat is optional if the job is -1 or -2
% mat is a square sparse matrice
% information are return in id fields
%
% Use help mumps_help for detailed information
%

errmsg = nargoutchk(1,1,nargout);
if(~isempty(errmsg))
     disp(errmsg);
     return;
end

arithtype = 1;

if(id.JOB == -2)     
     if(id.INST==-9999)
         disp('Uninitialized instance');
         return;
     end
     if(id.TYPE ~= arithtype)
       disp('You are trying to call z/d version on a d/z instance');
       return;
     end
     mumps(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS);
     id = [];
     return;
end


if(id.JOB == -1)
     if(id.INST~=-9999)
         disp('Allready initialized instance');
         return;
     end
     [inform,rinform,sol,inst,schur,redrhs,pivnul_list,sym_perm,uns_perm,icntl,cntl] = mumps(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS);
     id.INFOG = inform;
     id.RINFOG = rinform;
     id.SOL = sol;
     id.INST = inst;
     id.SCHUR = schur;
     id.REDRHS = redrhs;
     id.PIVNUL_LIST = pivnul_list;
     id.SYM_PERM = sym_perm;
     id.UNS_PERM = uns_perm;
     id.TYPE = arithtype;
     id.ICNTL=icntl;
     id.CNTL=cntl;
     return;
end

if(id.INST==-9999)
         disp('Uninitialized instance');
         return;
end

if(id.TYPE ~= arithtype)
   disp('You are trying to call z/d version on a d/z instance');
   return;
end

[inform,rinform,sol,inst,schur,redrhs,pivnul_list,sym_perm,uns_perm,icntl,cntl] = mumps(id.SYM,id.JOB,id.ICNTL,id.CNTL,id.PERM_IN,id.COLSCA,id.ROWSCA,id.RHS,id.VAR_SCHUR,id.INST,id.REDRHS,mat);
id.INFOG = inform;
id.RINFOG = rinform;
id.SOL = sol;
id.INST = inst;
if(id.JOB == 2 | id.JOB == 4 | id.JOB == 6)
  if(id.SYM == 0)
	id.SCHUR = schur';
  else
        id.SCHUR = triu(schur)+tril(schur',-1);
  end
end
id.REDRHS = redrhs;
id.PIVNUL_LIST = pivnul_list;
id.SYM_PERM(sym_perm) = [1:size(mat,1)];
id.UNS_PERM = uns_perm;
id.ICNTL=icntl;
id.CNTL=cntl;
