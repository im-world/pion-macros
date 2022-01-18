#!/usr/local/bin/tcsh
 setenv HOME /eic/u/$LOGNAME
 source /etc/csh.login
 foreach i (/etc/profile.d/*.csh)
   source $i
 end
 source $HOME/.login
 source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/eic_setup.csh -n

set FILE = condor.out
set STRING = condorjob done
set nJobs = 1000
set j = 1

#folder not created if pre-existing
while ($j < $nJobs)
    mkdir macros$j
    @ j++
end

echo Fin

