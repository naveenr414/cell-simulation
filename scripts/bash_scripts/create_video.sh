rm $1/video.mp4
cd $1; ffmpeg -f image2 -pattern_type glob -i '*.png' -an -c:v libx264 -r 10 video.mp4
