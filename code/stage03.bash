for f in *.json
do
    echo $f
    strun $JYP/drizzle.image3 $f
done
