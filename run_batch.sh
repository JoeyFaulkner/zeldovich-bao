limits=( 10 25 50 100)
for i in "${limits[@]}"
do
    python runZA.py production $i
done