#!/bin/sh
set -eu

prospino_dir=`git rev-parse --show-toplevel`/vendor/prospino2
input=prospino.in.les_houches
output=result.dat
run_avg=./pros_squark_avg.run 
run_uc=./pros_squark_uc.run 

make

if [ ! -d Pro2_subroutines ]; then
  ln -s ${prospino_dir}/Pro2_subroutines .
fi
rm -f $input $output

function run(){
  sed -e "s/<<<%MGLU>>>/$MGLU/g" \
      -e "s/<<<%MN1>>>/$MN1/g" \
      -e "s/<<<%MST>>>/$MST/g" \
      -e "s/<<<%MSB>>>/$MSB/g" \
      -e "s/<<<%MSQ2>>>/$MSQ2/g" \
      -e "s/<<<%MSQ4>>>/$MSQ4/g" \
      -e "s/<<<%MSQ8>>>/$MSQ8/g" \
      colored.slha.template | tee $input > /dev/null
  $1 >&2
  cat prospino.dat | grep '^\w'
}

function default(){
  MN1=100
  MGLU=7000
  MSQ2=1000
  MSQ4=1000
  MSQ8=1000
  MST=1000
  MSB=1000
}

function execute(){
  echo "# Gluino 7.0TeV, All squark 1.0TeV"
  default; run $run_avg; run $run_uc
  cp $input SLHA/mg7_sq1.slha

  echo "# Gluino 7.0TeV, All squark 1.0TeV but Stop 10.0TeV"
  MST=10000; run $run_avg; run $run_uc
  cp $input SLHA/mg7_sq1_st10.slha

  echo "# Gluino 7.0TeV, All squark 1.0TeV but Stop/Sbottom 10.0TeV"
  MSB=10000; run $run_avg; run $run_uc
  cp $input SLHA/mg7_sq1_st10_sb10.slha

  echo "# Gluino 100TeV, All squark 1.0TeV"
  default; MGLU=100000; run $run_avg; run $run_uc
  cp $input SLHA/mg100_sq1.slha

  echo "# Gluino 100.0TeV, All squark 1.0TeV but Stop 10.0TeV"
  MST=10000; run $run_avg; run $run_uc
  cp $input SLHA/mg100_sq1_st10.slha

  echo "# Gluino 100.0TeV, All squark 1.0TeV but Stop/Sbottom 10.0TeV"
  MSB=10000; run $run_avg; run $run_uc
  cp $input SLHA/mg100_sq1_st10_sb10.slha

  echo "# Gluino 100.0TeV, right-handed 1.0TeV and others 10.0TeV"
  MSQ8=10000; run $run_uc
  cp $input SLHA/mg100_sqR1_sq10.slha

  echo "# Gluino 100.0TeV, uR+cR 1.0TeV and others 10.0TeV"
  MSQ4=10000; run $run_uc
  cp $input SLHA/mg100_uRcR1_sq10.slha
}

execute | tee $output
