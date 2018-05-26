#this script runs the DataJK matlab script on the 12 datafiles specified below
# Matlab R2017b or later is recommended
echo "DataJK.m will be called for the following particles in order of appearance:"
names="kaonstar_us pion_ud proton_12_uud kaon_us delta_32_uuu omega_32_sss phi_ss rho_ud sigma_12_uus sigmastar_32_uus xi_12_uss"
i=0;
for name in $names
do
	i=$((i+1));
	printf "%i\t%s\n" $i "$name"
done

printf "\n"

if [ ! -d results ]; then
	mkdir results
fi

echo "Proceeding to run the script."
echo "MATLAB will abort job when Ctrl+C is pressed"

printf "\n"

n=0
for name in $names
do
	n=$((n+1));
	printf "Running\t%i of %i \t %s\n" $n $i "$name"
	if [ ! -d results/$name ]; then
		mkdir results/$name	
	fi

	datafile=data/$name'_correl_all_data_ascii.dat'
	matlab -r "DataJK '$datafile' '$name' ; exit"
done
