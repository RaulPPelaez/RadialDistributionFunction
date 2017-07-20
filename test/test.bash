



exe=$1
echo "Testing with a random distribution..."

echo "3 snapshots of 2e5 particles"

echo "A q2D distribution"
octave -qf random_positions.m 2>/dev/null | awk '{print $1, $2, $3*0.1}' | $exe -N 200000 -Nsnapshots 3 -nbins 1000 -rcut 0.5 -L 1  -dim q2D > q2D.gdr

echo "A true 2D random distribution"
octave -qf random_positions.m 2>/dev/null| $exe -N 200000 -Nsnapshots 3 -nbins 1000 -rcut 0.5 -L 1  -dim 2D > 2D.gdr

echo "A 3D random distribution"
octave -qf random_positions.m 2>/dev/null | $exe -N 200000 -Nsnapshots 3 -nbins 1000 -rcut 0.5 -L 1  -dim 3D > 3D.gdr

