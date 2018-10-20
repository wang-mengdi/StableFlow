import os

os.system("ffmpeg -r 20 -i images/%08d.png -crf 0 -y fluid.mov")
