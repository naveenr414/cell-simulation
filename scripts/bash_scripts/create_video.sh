rm ../src/data_film2/video.mp4
cd ../src/data_film2; ffmpeg -f image2 -pattern_type glob -i '*.png' -an -c:v libx264 -r 10 video.mp4
