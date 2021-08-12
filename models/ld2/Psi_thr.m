function [Psi] = Psi_thr(pus,pds,param)
Pi = pds./pus;
   
Psi = param.psi_thr.c1.*sqrt( 1 - (Pi.*param.psi_thr.c2).^param.psi_thr.c3 ) + param.psi_thr.c4;

end