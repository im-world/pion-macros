#!/usr/local/bin/tcsh
 setenv HOME /eic/u/$LOGNAME
 source /etc/csh.login
 foreach i (/etc/profile.d/*.csh)
   source $i
 end
 source $HOME/.login
 source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/eic_setup.csh -n

set nJobs = 1000
set j = 0
set count = 1

mkdir EvalFiles
mv macros macros0

while ($j < $nJobs)
   # echo inloop
    cd macros$j
    if ( -f Eval_CEMC.root  ) then
	 cp Eval_CEMC.root ../EvalFiles/Eval_CEMC_$count.root
	 cp Eval_FEMC.root ../EvalFiles/Eval_FEMC_$count.root
	 cp Eval_EEMC.root ../EvalFiles/Eval_EEMC_$count.root
	 cp Eval_FHCAL.root ../EvalFiles/Eval_FHCAL_$count.root
	 cp Eval_HCALIN.root ../EvalFiles/Eval_HCALIN_$count.root
	 cp Eval_HCALOUT.root ../EvalFiles/Eval_HCALOUT_$count.root
	 echo extracted from macros$j
	 @ count++
    else
         echo Empty Directory
    endif
    cd ../
    @ j++
end

sed -i "s/successfulJobs/$count/g" hadd.C
echo $count
echo Fin
