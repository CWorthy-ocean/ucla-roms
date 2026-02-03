#!/bin/bash

# need this here as well, in case example run on its own:
case "$1" in
    expanse|maya|laptop|github_gnu|github_ifx)
	echo "running test for $1"
	;;
    *)
	echo "Script must have argument 'expanse' or 'maya'! E.g.: './do_test_all.sh maya'. Try again!"
	exit 1
	;;
esac

bm_file="benchmark.result_$1"                    # set benchmark specific to machine (maya/expanse)
echo "$bm_file"

# 1) Compile test case:
echo "##############################"
echo "  test compiling..."
echo "##############################"

make compile_clean &> /dev/null
make BUILD_MODE=test


# 2) Run test case:
echo "##############################"
echo "  test running..."
echo "##############################"

if [ "$1" = "expanse" ]
then
    srun --mpi=pmi2 -n 6 ./roms benchmark.in 2>&1 | tee -a test.log
else
    mpirun -n 6 ./roms benchmark.in 2>&1 | tee -a test.log
fi

# 2) Python - confirm values:
cp $ROMS_ROOT/tests/scripts/test_roms.py .
python3 test_roms.py $bm_file
retval=$?

rm test_roms.py

# 3) Rename results logs so they can't be mistakenly read by the
#    python script even if new simulation doesn't run
mv test.log test_old.log

exit $retval                                     # pass success value onto do_test_all script

# Notes:
# Compile directories are left in the various examples for fast re-compilation of all
#  test cases.






