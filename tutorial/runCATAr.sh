for a in Imidazole 
do
	cd $a
	for i in 1 2 3
	do
		cd $i
		echo "1
		0" > choices1.txt
		gmx trjconv -f prod1.xtc -s prod1.tpr -o dynamic_clust.xtc -pbc cluster  < choices1.txt
		gmx trjconv -f dynamic_clust.xtc -s prod1.tpr -o dynamic_njump.xtc -pbc mol -center < choices1.txt
echo "1 | 20
q
" > index.txt
		gmx make_ndx -f prod1.gro -o index.ndx < index.txt
		echo "24" > 24.txt
		gmx trjconv -s prod1.tpr -f dynamic_njump.xtc -o protein.xtc -n index.ndx < 24.txt
		echo "0" > 0.txt
		gmx convert-tpr -s prod1.tpr -o protein.tpr -n index.ndx < 24.txt
		gmx trjconv -s protein.tpr -f protein.xtc -o protein.gro -dump 0 < 0.txt
		echo "3
		1
		0" > choices.txt
		gmx trjconv -s protein.tpr -f protein.xtc -o protein-fit.xtc -fit rot+trans -center < choices.txt
		rm dynamic_njump.xtc
		rm protein.xtc
		rm dynamic_clust.xtc
		rm choices.txt
		rm 0.txt
		rm 24.txt
		rm choices1.txt
		rm index.txt
		cp *GMX.itp $a.itp
		cd ..
	done
	for i in 1 2 3
	do
		mkdir snapsCAT
		cd $i
		echo "0" > 0.txt
		gmx trjconv -f protein-fit.xtc -s protein.tpr -o ../snapsCAT/protein_$i.pdb -sep -skip 2 < 0.txt
		cd ..
	done
	cd snapsCAT
	for e in protein*pdb;
	do 
		echo $e;
	done > list.dat
	mkdir 0.4 0.5 0.8
	for i in 0.4 0.5 0.8 1 2 3 4 5 6
	do
	/home/quicosabanes/CAT_git/CAT -i list.dat -o output.dat -c 10 -n 5000 -t ../1/topol.top -e $i -s ../1/$a.itp -v $i -g 5
	mv all* $i
	done
cd ../..
done

