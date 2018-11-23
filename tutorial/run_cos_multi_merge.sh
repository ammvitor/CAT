for i in 1 2 3
do
cd $i
echo "6
1" > choices.txt
gmx pdb2gmx -f ar_mod.pdb -o prot.gro -merge all < choices.txt
rm choices.txt 
gmx editconf -f prot.gro -o box.gro -c -d 1.0 -bt cubic
gmx insert-molecules -f box.gro -ci acetamide_GMX.gro -o boxLig.gro -nmol 469
sed -i '/forcefield.itp/ a\  #include "acetamide_GMX.itp"' topol.top
sed -i '/forcefield.itp/ a\  ; Include acetamide_GMX.itp topology' topol.top
sed -i '$ a\acetamide          469' topol.top
gmx solvate -cp boxLig.gro -cs spc216.gro -o solvbox.gro -p topol.top
gmx grompp -f ions.mdp -c solvbox.gro -p topol.top -o ions.tpr
echo "15" > choices.txt
gmx genion -s ions.tpr -o sbi.gro -p topol.top -pname NA -nname CL -conc 0.1 -neutral < choices.txt
rm choices.txt
cd ..
done

