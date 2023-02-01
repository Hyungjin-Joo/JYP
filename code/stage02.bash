rm *json
asn_make_pool st2.ecsv *rate.fits
asn_generate st2.ecsv
rm *image3*
for f in *json
do
	echo $f
	strun calwebb_image2 $f --steps.flat_field.skip True --steps.resample.skip  True
done
