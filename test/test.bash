



exe=$1
otherArgs=$(echo "$@" | awk '{$1=""; print $0}')


args=$( echo "-N 10000 -Nsnapshots 100 -nbins 100 -rcut 1 -L 10 -outputDecimals 6 -fixBIAS $otherArgs")
echo "Testing with a random distribution..."

echo "3 snapshots of 1e4 particles"

echo "A q2D distribution"
echo "$exe $args  -dim q2D > q2D.gdr"
cat inipos |
    awk '{print $1, $2, $3*0.1}' |
    { /usr/bin/time  -p sh -c "$exe $args  -dim q2D > q2D.gdr"; }  2>&1 | paste -sd" " | awk '{print "q2D :", $0}' >> timing.dat

echo "A true 2D random distribution"
cat inipos |
    { /usr/bin/time  -p sh -c "$exe $args -dim 2D > 2D.gdr" ; }  2>&1 | paste -sd" " | awk '{print "2D :", $0}' >> timing.dat
echo "$exe $args  -dim 2D > 2D.gdr"
echo "A 3D random distribution"
cat inipos |
    { /usr/bin/time  -p sh -c "$exe $args  -dim 3D > 3D.gdr"; }  2>&1 | paste -sd" " | awk '{print "3D :", $0}' >> timing.dat

echo "$exe $args  -dim 3D > 3D.gdr"
