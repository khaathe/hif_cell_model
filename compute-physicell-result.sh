#/usr/bin/bash

last_povwriter_index=$1
result_dir="./result"
log_file="$result_dir/log.txt"

if [ -d "$result_dir" ]; then 
    rm -R $result_dir
fi

echo ">>>>> mkdir $result_dir"
echo ">>>>> mkdir $result_dir/pov $result_dir/png $result_dir/gif"

mkdir $result_dir
mkdir $result_dir/pov $result_dir/png $result_dir/gif

if [$# lt 1]; then
    povwriter
else
    echo ">>>>> povwriter 0:1:$last_povwriter_index"
    povwriter 0:1:$last_povwriter_index
fi

echo ">>> mv *.pov $result_dir/pov"
mv *.pov "$result_dir/pov"

for  pov_file in `ls $result_dir/pov`; do
    png_file=$(echo "$pov_file"|sed 's/pov/png/')
    echo ">>> creating png : $result_dir/png/$png_file"
    povray +O$result_dir/png/$png_file  $result_dir/pov/$pov_file
done

echo "magick convert $result_dir/png/*.png $result_dir/gif/simu.gif"
magick convert $result_dir/png/*.png $result_dir/gif/simu.gif
