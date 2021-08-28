function [Psi] = Psi_wg(pus,pds,param)
Pi = pds./pus;

Psi = param.psi_wg.c1.*sqrt( 1 - (Pi.*param.psi_wg.c2).^param.psi_wg.c3 ) + param.psi_wg.c4;

end