for f in *cal.fits
do
    echo "$f"
    IFS='.'
    read -ra starr<<<"$f"
    sex ${starr[0]}.${starr[1]} -c $JYP/code/default.sex -CATALOG_NAME ${starr[0]}.cat -CHECKIMAGE_NAME ${starr[0]}_check.fits
    echo "SExtractor Done"
    rm *cat
done

echo "mask expansion"
python $JYP/code/mask_expansion.py
echo "sky estimation"
python $JYP/code/sky_estimation.py
