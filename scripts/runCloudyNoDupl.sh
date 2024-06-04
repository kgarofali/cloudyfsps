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
  linefix='.lineflux'
  linefile="$dir$2$i$linefix"
  if [ -f $linefile ];
  then
    echo "File $linefile exists."
  else
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
  fi
  echo "Done."
done







#!/bin/sh
echo $PATH
echo $PYTHONPATH
echo $CLOUDY_EXE
suffix='.in'
mod_prefix=$2$3
infile=$mod_prefix$suffix
dir=$1
echo $dir
echo $infile
suffix='.out'
outfile="$dir/$mod_prefix$suffix"
if [ -f $outfile ];
then
  echo "File $outfile exists."
  echo $outfile
else
  echo "Running CLOUDY on $infile"
  cd $dir && $CLOUDY_EXE -r $mod_prefix
  echo "CLOUDY finished $infile"
  echo "Running python..."
  python /Users/kgarofal/software/cloudyfsps/scripts/runCloudy.py $outfile
  echo "Python finished."
fi
echo "Done."
