#!/bin/bash

Setup(){
echo "Sending the process name $process_name"
{
echo "$process_name" 
cat
} | ./setup.x
}
## VMC part

file_indexer=$(find . -type f -name "*.vmc*" -exec grep -Fl "$process_name" {} + | wc -l)
weight_col=-1

Input_QMC(){
	echo "What VMC arguments do you want? To deliver the arguments press Ctrl+D. (it may take sometime)"
vmc_args=$(cat)
{
echo "$vmc_args"
} | cat > "$process_name.in"
if [[ $vmc_args =~ "restart" ]];then
	echo "$process_name" | ./qmc.x >> "${process_name}$file_indexer.vmc"
else
	(( file_indexer+=1 ))
	echo "$process_name"|./qmc.x > "$process_name$file_indexer.vmc"

fi

if [[ $vmc_args =~ "optimize" ]];then

	echo "Convergence code (repeat sims until info != -1)"
	echo $(grep info "$process_name$file_indexer.vmc")
fi
}

Analysis(){
	echo "We now do the analysis of our vmc/dmc sim"
	local chosen_index=${1:-$file_indexer}
	echo "Energy\n"
	if [[ $weight_col == -1 ]];then
		Find_Weight "$process_name$chosen_index.vmc"
	fi
	grep elocal "$process_name$chosen_index.vmc" | awk "{print \$1, \$$weight_col}"  | ./statforw.x | tee "${process_name}${chosen_index}_elocal.out"
	echo "Acceptance rate\n"
	grep "acc.rate" "$process_name$chosen_index.vmc" | awk "{print \$1, \$$weight_col}" | ./statforw.x| tee "${process_name}${chosen_index}_accrate.out"
	gnuplot -persist <<-EOF 
	plot '< grep elocal $process_name${chosen_index}.vmc' u 0:1 with lines
	EOF
}

Find_Weight(){
	local file=$1
	for n in {2..4}; do
		output1=$(grep elocal "$file" | awk "{print \$$n}")
		output2=$(grep "acc.rate" "$file" | awk "{print \$$n}")
		if [[ "$output1" == "$output2" && -n "$output1"  ]];then
			weight_col=$n
		fi
	done
	echo "Collumn found for weights: $weight_col"
}

Delete(){
	find . -type f -name "*$process_name*"
	echo "These are the files that will be deleted. Are you sure?(1:yes,anything else:no)"
	local answer
	read answer
	find . -name "*$process_name*" -type f -exec rm {} +
	echo "Deleted!"
}

choice=-1

while true ; do
	
	if [[ -z "$process_name" ]]; then
		echo "Type the name of your program"
		read process_name
		file_indexer=$(find . -type f -name "*.vmc*" -exec grep -Fl "$process_name" {} + | wc -l)

	fi


	read -p "What to do? Setup:0 Simulate:1 Analyse:2 Delete_Sims:3 Quit: anything else\n" choice
	if [[ $choice == 0 ]]; then
		Setup
	elif [[ $choice == 1 ]]; then
		Input_QMC
	elif [[ $choice == 2 ]]; then
		read -p "Choose an index to analyse: (max: $file_indexer)\n" ch_ind
		Analysis $ch_ind
	elif [[ $choice == 3 ]]; then
		Delete
	else
		break
	fi
done
