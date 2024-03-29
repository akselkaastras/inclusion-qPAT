function cmdsub = submitToDTUCluster(jobid,cmd,email)

id = num2str(jobid);
jobname = ['job_', id];

memcore = 6000;
maxmem  = 10000;
ncores = 1;

% update jobscript string
str = '#!/bin/sh\n\n#BSUB -q hpc\n';
str = append(str,['#BSUB -J ', jobname, '\n']);
str = append(str,['#BSUB -n ', num2str(ncores), '\n']);
str = append(str,'#BSUB -R "span[hosts=1]"\n');
str = append(str,['#BSUB -R "rusage[mem=', num2str(memcore),'MB]"\n']);
str = append(str,['#BSUB -M ', num2str(maxmem), 'MB\n']);
str = append(str,'#BSUB -W 64:00\n');
str = append(str,['#BSUB -u ', email,'\n']);
str = append(str,'#BSUB -N \n');
str = append(str,['#BSUB -o Driverscripts/Output/OutputCluster_',id,'_.out\n']);
str = append(str,['#BSUB -e Driverscripts/Error/ErrorCluster_',id,'_.err\n']);
%str = append(str,'export OMP_NUM_THREADS=$LSB_DJOB_NUMPROC\n');
%str = append(str,['module load gcc\n']);
str = append(str,cmd);

% name
jobscript = ['Driverscripts/Jobscripts/submitToCluster_',jobid,'.sh'];
fileid = fopen(jobscript,'w');
fprintf(fileid,str);
fclose(fileid);

% Submit job
cmdsub = ['bsub < ', jobscript];
[~, cmdout] = system(cmdsub);
delete(jobscript)
