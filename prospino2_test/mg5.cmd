#!/bin/sh

MG5=/Users/misho/HEPcodes/MG5_aMC/bin/mg5_aMC

$MG5 <<_EOC_
import model mssm
define sucR = ur ur~ cr cr~
define q = d d~ u u~ s s~ c c~ b b~
generate    g g > sucR sucR / go
add process q q > sucR sucR / go
output uRcRnogo -f
launch uRcRnogo --laststep=parton
set lhc 13
set pdlabel cteq6l1
SLHA/mg100_sq1.slha
done
_EOC_

for process in "ur ur" "ur ur~" "cr cr" "cr cr~" "ur cr" "ur cr~" "cr ur~"; do
  name=`echo uRcR_$process | sed -e 's/ //' -e 's/~/x/'`
  $MG5 <<_EOC_
import model mssm
define sucR = ur ur~ cr cr~
generate p p > $process
output $name -f
launch $name --laststep=parton
set lhc 13
set pdlabel cteq6l1
SLHA/mg100_sq1.slha
done
launch $name --laststep=parton
set lhc 13
set pdlabel cteq6l1
SLHA/mg7_sq1.slha
done
_EOC_
done
