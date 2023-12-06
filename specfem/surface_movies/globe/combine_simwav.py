"""
Utility function to combine two sets of .png frames into a single tiled .png,
and then combine everything into a .gif. Uses ImageMagick under the hood.
"""
import os
from glob import glob
from subprocess import run

movies = "simulation/*.png"
waveforms = "waveform/*.png"

mov_files = sorted(glob(movies))
wav_files = sorted(glob(waveforms))
assert(mov_files)
# assert(len(mov_files) == len(wav_files))

for i, (m, w) in enumerate(zip(mov_files, wav_files)):
    callcmd = f"montage -tile 2x1 {m} {w} -geometry +0+0 ./combine/{i:0>3}_mov.png"
    run(callcmd.split(" "))

os.chdir("combine")
callcmd = f"convert -delay 0 -loop 1 *.png ../output/simwav.gif"
run(callcmd.split(" "))



