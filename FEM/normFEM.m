function normsq = normFEM(v,fmdl)

normsq = v'*fmdl.Carea*v;