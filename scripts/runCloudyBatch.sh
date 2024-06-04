#!/bin/sh
echo $PATH
echo $PYTHONPATH
echo $CLOUDY_EXE
dir=$1
suffix='.in'
FILES=$2*$suffix
nofiles=$(cd $dir && ls $FILES | wc -l)
for i in $(seq 1 $nofiles)
do
  cd $dir && $CLOUDY_EXE -r $2$i
  echo "CLOUDY finished $2$i"
  outfix='.out'
  outfile="$dir$2$i$outfix"
  echo $outfile
  if [ -f $outfile ];
  then
    echo "File $outfile exists."
    echo "Running python..."
    python /Users/kgarofal/software/cloudyfsps/scripts/runCloudy.py $outfile
    echo "Python finished."
  else
    echo "ERROR: file $outfile does NOT exist."
  fi
  echo "Done."
done
