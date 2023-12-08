# Use Imagemagick to crop whitespace from ALL .png files in a directory
# Acts in place
for filename in *.png; do
	echo ${filename}
	convert ${filename} -fuzz 1% -trim +repage ${filename}
done
