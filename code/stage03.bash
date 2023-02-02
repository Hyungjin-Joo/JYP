asn_make_pool st3.ecsv *cal.fits
asn_generate st3.ecsv
rm *image2*  
for f in *.json
do
    echo $f
    strun $JYP/drizzle.image3 $f
done
