echo "Title of movie?"
read title

rm out.mp4
rm final.mp4

convert -size 500x270 -background transparent -fill black -pointsize 24 -font Helvetica label:"\n\n\n\n\n\n${title}\n" logo.png

ffmpeg -framerate 60 -pattern_type glob -i 'fr_*.png' \
  -c:a copy -shortest -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4

ffmpeg -i out.mp4 -i logo.png -filter_complex "overlay=x=(10):y=(main_h-overlay_h-10)" final.mp4
