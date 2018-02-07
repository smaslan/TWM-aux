%% Starts multicore slave in desired folder.

function [] = run_mc_slave(shr_fld)

  %% get script folder path
  rpth = mfilename('fullpath');
  rpth = rpth(1:strchr(rpth,filesep(),1,'last')-1);
  
  %% set current path
  cd(rpth);
  addpath(rpth);
    
  % try to set process priority to "below normal"
  if(ispc())
    syscmd = ['cmd /Q /C "wmic process where handle=' int2str(getpid()) ' CALL SetPriority "Below Normal" "'];
    [sout,stxt] = system(syscmd);            
  endif
  
  %% set cwd to script folder (since version 3.6.4 the Octave won't find data at '--path')
  mfld = mfilename('fullpath'); mfld = mfld(1:strchr(mfld,filesep(),1,'last'));
  addpath([mfld 'multicore']);
  
  %% run slave in this process as well
  startmulticoreslave(shr_fld);

endfunction