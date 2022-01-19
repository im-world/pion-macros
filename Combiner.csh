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

echo Deletion done


set j = 1

#folder not created if pre-existing
while ($j < $nJobs)
    mkdir macros$j
    @ j++
end

echo Creation done

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
echo Done Moving


 root.exe -q -b hadd.C\(\"EEMC\"\)
 root.exe -q -b hadd.C\(\"CEMC\"\)
 root.exe -q -b hadd.C\(\"FEMC\"\)
 root.exe -q -b hadd.C\(\"HCALIN\"\)
 root.exe -q -b hadd.C\(\"HCALOUT\"\)
 root.exe -q -b hadd.C\(\"FHCAL\"\)

echo Statistics Combined


