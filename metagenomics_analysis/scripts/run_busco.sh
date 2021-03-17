#Wrapper script used to run BUSCO4 on all bins

for file in ../binning/TS1974/concoct_bins/*.fa
do
  tmp=${file##*/}
  bin=${tmp%.*}
  busco -i "$file" -o "$bin" --out_path ./TS1974 -m genome -c 2 -f
  #echo "$bin"
done