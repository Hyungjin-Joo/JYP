#for f in *cal.fits
#do
#    echo "$f"
#    IFS='.'
#    read -ra starr<<<"$f"
#    conda run -n astroconda sex ${starr[0]}.${starr[1]} -c /Volumes/Metal_Empire_001/JWST/code/default.sex -CATALOG_NAME ${starr[0]}.cat -CHECKIMAGE_NAME ${starr[0]}_check.fits
#    echo "SExtractor Done"
##    rm *cat
#done

#echo "mask expansion"
#conda run -n astroconda python /Volumes/Metal_Empire_001/JWST/code/mask_expansion.py
echo "sky estimation"
conda run -n astroconda python /Volumes/Metal_Empire_001/JWST/code/sky_estimation.py
