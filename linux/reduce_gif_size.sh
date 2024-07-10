# Uses Imagemagcik to reduce file size on gifs, fuzz blends similar colors 
# together. Can also be used on images I guess
mogrify -layers 'optimize' -fuzz 5% ${1}
