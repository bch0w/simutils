# Download a large file from google drive based on the file id

# FID="188VDHS1VdMUYzJkzITZ80hZFGpL90lAe"  # Carl 1
# FID="1dBMWhK9njipZvhegtdeYvNLHIND9zjEF"
FID="1qLuxyyNaZoD1hLeZrYpykAFlgfsDuY4"
FIDOUT="nz_north_2019.tar"

curl -L -c cookies.txt 'https://docs.google.com/uc?export=download&id='$FID \
     | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p' > confirm.txt

curl -L -b cookies.txt -o $FIDOUT \
     'https://docs.google.com/uc?export=download&id='$FID'&confirm='$(<confirm.txt)

rm -f confirm.txt cookies.txt
