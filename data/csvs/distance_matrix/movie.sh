rm out.mp4

ffmpeg -framerate 60 -pattern_type glob -i 'fr_*.png' \
  -c:a copy -shortest -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
