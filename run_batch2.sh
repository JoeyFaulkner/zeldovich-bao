limits=( 40 10 25 50 100)
for i in "${limits[@]}"
do
for j in $(seq 1 100)
do

    python runZA.py nice_run $j $i
done
done