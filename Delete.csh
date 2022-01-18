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

while ($j < $nJobs)
   # echo inloop
    cd macros$j
    grep -q "$STRING" $FILE
    if ( $status != 0 ) then
	 # echo then
	 cd ../
         # echo cd
 	 mv macros$j ../brownJobs/macros$j
	 echo Deleted macros$j
    else
         # echo found
	 grep -n "$STRING" $FILE
         cd ../
    endif
    @ j++
end

echo Fin
