#include "receptor_pdb.h"

receptor_pdb::receptor_pdb()
{

/*
Br          35      79.90    0.0000  A   0.00000e+00  0.00000e+00
C            6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CA           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CB           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CC           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CK           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CM           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CN           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CQ           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CR           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CT           6      12.01    0.0000  A   3.39967e-01  4.57730e-01
CV           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
CW           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
C*           6      12.01    0.0000  A   3.39967e-01  3.59824e-01
C0          20      40.08    0.0000  A   3.05240e-01  1.92376e+00
F            9      19.00    0.0000  A   3.11815e-01  2.55224e-01
H            1       1.008   0.0000  A   1.06908e-01  6.56888e-02
HC           1       1.008   0.0000  A   2.64953e-01  6.56888e-02
H1           1       1.008   0.0000  A   2.47135e-01  6.56888e-02
H2           1       1.008   0.0000  A   2.29317e-01  6.56888e-02
H3           1       1.008   0.0000  A   2.11499e-01  6.56888e-02
HA           1       1.008   0.0000  A   2.59964e-01  6.27600e-02
H4           1       1.008   0.0000  A   2.51055e-01  6.27600e-02
H5           1       1.008   0.0000  A   2.42146e-01  6.27600e-02
HO           1       1.008   0.0000  A   0.00000e+00  0.00000e+00
HS           1       1.008   0.0000  A   1.06908e-01  6.56888e-02
HW           1       1.008   0.0000  A   0.00000e+00  0.00000e+00
HP           1       1.008   0.0000  A   1.95998e-01  6.56888e-02
I           53     126.9     0.0000  A   4.18722e-01  1.67360e+00
Cl          17      35.45    0.0000  A   4.40104e-01  4.18400e-01
Na          11      22.99    0.0000  A   3.32840e-01  1.15897e-02
IB           0     131.0     0.0000  A   8.90899e-01  4.18400e-01
MG          12      24.305   0.0000  A   1.41225e-01  3.74342e+00
N            7      14.01    0.0000  A   3.25000e-01  7.11280e-01
NA           7      14.01    0.0000  A   3.25000e-01  7.11280e-01
NB           7      14.01    0.0000  A   3.25000e-01  7.11280e-01
NC           7      14.01    0.0000  A   3.25000e-01  7.11280e-01
N2           7      14.01    0.0000  A   3.25000e-01  7.11280e-01
N3           7      14.01    0.0000  A   3.25000e-01  7.11280e-01
N*           7      14.01    0.0000  A   3.25000e-01  7.11280e-01
O            8      16.00    0.0000  A   2.95992e-01  8.78640e-01
OW           8      16.00    0.0000  A   3.15061e-01  6.36386e-01
OH           8      16.00    0.0000  A   3.06647e-01  8.80314e-01
OS           8      16.00    0.0000  A   3.00001e-01  7.11280e-01
O2           8      16.00    0.0000  A   2.95992e-01  8.78640e-01
P           15      30.97    0.0000  A   3.74177e-01  8.36800e-01
S           16      32.06    0.0000  A   3.56359e-01  1.04600e+00
SH          16      32.06    0.0000  A   3.56359e-01  1.04600e+00
CU          29      63.55    0.0000  A   3.39967e-01  3.59824e-01
FE          26      55.00    0.0000  A   0.00000e+00  0.00000e+00
K           19      39.10    0.0000  A   4.73602e-01  1.37235e-03
Rb          37      85.47    0.0000  A   5.26699e-01  7.11280e-04
Cs          55     132.91    0.0000  A   6.04920e-01  3.37230e-04
*/



    this->atom_names.push_back("Br");
    this->atom_names.push_back("C");
    this->atom_names.push_back("CA");
    this->atom_names.push_back("CB");
    this->atom_names.push_back("CC");
    this->atom_names.push_back("CK");
    this->atom_names.push_back("CM");
    this->atom_names.push_back("CN");
    this->atom_names.push_back("CQ");
    this->atom_names.push_back("CR");
    this->atom_names.push_back("CT");
    this->atom_names.push_back("CV");
    this->atom_names.push_back("CW");
    this->atom_names.push_back("C*");
    this->atom_names.push_back("C0");
    this->atom_names.push_back("F");
    this->atom_names.push_back("H");
    this->atom_names.push_back("HC");
    this->atom_names.push_back("H1");
    this->atom_names.push_back("H2");
    this->atom_names.push_back("H3");
    this->atom_names.push_back("HA");
    this->atom_names.push_back("H4");
    this->atom_names.push_back("H5");
    this->atom_names.push_back("HO");
    this->atom_names.push_back("HS");
    this->atom_names.push_back("HW");
    this->atom_names.push_back("HP");
    this->atom_names.push_back("I");
    this->atom_names.push_back("Cl");
    this->atom_names.push_back("Na");
    this->atom_names.push_back("IB");
    this->atom_names.push_back("MG");
    this->atom_names.push_back("N");
    this->atom_names.push_back("NA");
    this->atom_names.push_back("NB");
    this->atom_names.push_back("NC");
    this->atom_names.push_back("N2");
    this->atom_names.push_back("N3");
    this->atom_names.push_back("O");
    this->atom_names.push_back("OW");
    this->atom_names.push_back("OH");
    this->atom_names.push_back("OS");
    this->atom_names.push_back("O2");
    this->atom_names.push_back("P");
    this->atom_names.push_back("S");
    this->atom_names.push_back("SH");
    this->atom_names.push_back("CU");
    this->atom_names.push_back("FE");
    this->atom_names.push_back("K");
    this->atom_names.push_back("Rb");
    this->atom_names.push_back("Cs");


    this->atom_masses.push_back(79.90);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(12.01);
    this->atom_masses.push_back(40.08);
    this->atom_masses.push_back(19.00);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(1.008);
    this->atom_masses.push_back(126.9);
    this->atom_masses.push_back(35.45);
    this->atom_masses.push_back(22.99);
    this->atom_masses.push_back(131.0);
    this->atom_masses.push_back(24.305);
    this->atom_masses.push_back(14.01);
    this->atom_masses.push_back(14.01);
    this->atom_masses.push_back(14.01);
    this->atom_masses.push_back(14.01);
    this->atom_masses.push_back(14.01);
    this->atom_masses.push_back(14.01);
    this->atom_masses.push_back(14.01);
    this->atom_masses.push_back(16.00);
    this->atom_masses.push_back(16.00);
    this->atom_masses.push_back(16.00);
    this->atom_masses.push_back(16.00);
    this->atom_masses.push_back(16.00);
    this->atom_masses.push_back(30.97);
    this->atom_masses.push_back(32.06);
    this->atom_masses.push_back(32.06);
    this->atom_masses.push_back(63.55);
    this->atom_masses.push_back(55.00);
    this->atom_masses.push_back(39.10);
    this->atom_masses.push_back(85.47);
    this->atom_masses.push_back(132.91);


    this->atom_sigma.push_back(0.00000e+00);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(3.05240e-01);
    this->atom_sigma.push_back(3.11815e-01);
    this->atom_sigma.push_back(1.06908e-01);
    this->atom_sigma.push_back(2.64953e-01);
    this->atom_sigma.push_back(2.47135e-01);
    this->atom_sigma.push_back(2.29317e-01);
    this->atom_sigma.push_back(2.11499e-01);
    this->atom_sigma.push_back(2.59964e-01);
    this->atom_sigma.push_back(2.51055e-01);
    this->atom_sigma.push_back(2.42146e-01);
    this->atom_sigma.push_back(0.00000e+00);
    this->atom_sigma.push_back(1.06908e-01);
    this->atom_sigma.push_back(0.00000e+00);
    this->atom_sigma.push_back(1.95998e-01);
    this->atom_sigma.push_back(4.18722e-01);
    this->atom_sigma.push_back(4.40104e-01);
    this->atom_sigma.push_back(3.32840e-01);
    this->atom_sigma.push_back(8.90899e-01);
    this->atom_sigma.push_back(1.41225e-01);
    this->atom_sigma.push_back(3.25000e-01);
    this->atom_sigma.push_back(3.25000e-01);
    this->atom_sigma.push_back(3.25000e-01);
    this->atom_sigma.push_back(3.25000e-01);
    this->atom_sigma.push_back(3.25000e-01);
    this->atom_sigma.push_back(3.25000e-01);
    this->atom_sigma.push_back(3.25000e-01);
    this->atom_sigma.push_back(2.95992e-01);
    this->atom_sigma.push_back(3.15061e-01);
    this->atom_sigma.push_back(3.06647e-01);
    this->atom_sigma.push_back(3.00001e-01);
    this->atom_sigma.push_back(2.95992e-01);
    this->atom_sigma.push_back(3.74177e-01);
    this->atom_sigma.push_back(3.56359e-01);
    this->atom_sigma.push_back(3.56359e-01);
    this->atom_sigma.push_back(3.39967e-01);
    this->atom_sigma.push_back(0.00000e+00);
    this->atom_sigma.push_back(4.73602e-01);
    this->atom_sigma.push_back(5.26699e-01);
    this->atom_sigma.push_back(6.04920e-01);



    this->atom_epsilon.push_back(0.00000e+00);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(4.57730e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(1.92376e+00);
    this->atom_epsilon.push_back(2.55224e-01);
    this->atom_epsilon.push_back(6.56888e-02);
    this->atom_epsilon.push_back(6.56888e-02);
    this->atom_epsilon.push_back(6.56888e-02);
    this->atom_epsilon.push_back(6.56888e-02);
    this->atom_epsilon.push_back(6.56888e-02);
    this->atom_epsilon.push_back(6.27600e-02);
    this->atom_epsilon.push_back(6.27600e-02);
    this->atom_epsilon.push_back(6.27600e-02);
    this->atom_epsilon.push_back(0.00000e+00);
    this->atom_epsilon.push_back(6.56888e-02);
    this->atom_epsilon.push_back(0.00000e+00);
    this->atom_epsilon.push_back(6.56888e-02);
    this->atom_epsilon.push_back(1.67360e+00);
    this->atom_epsilon.push_back(4.18400e-01);
    this->atom_epsilon.push_back(1.15897e-02);
    this->atom_epsilon.push_back(4.18400e-01);
    this->atom_epsilon.push_back(3.74342e+00);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(8.78640e-01);
    this->atom_epsilon.push_back(6.36386e-01);
    this->atom_epsilon.push_back(8.80314e-01);
    this->atom_epsilon.push_back(7.11280e-01);
    this->atom_epsilon.push_back(8.78640e-01);
    this->atom_epsilon.push_back(8.36800e-01);
    this->atom_epsilon.push_back(1.04600e+00);
    this->atom_epsilon.push_back(1.04600e+00);
    this->atom_epsilon.push_back(3.59824e-01);
    this->atom_epsilon.push_back(0.00000e+00);
    this->atom_epsilon.push_back(1.37235e-03);
    this->atom_epsilon.push_back(7.11280e-04);
    this->atom_epsilon.push_back(3.37230e-04);



}

void receptor_pdb::topology_parser(string topol_file){

    vector< int> natoms_per_residue;
    int natoms = 0;
    char buffer[900];
    FILE * ftopol_rec;
    ftopol_rec = fopen (topol_file.c_str(),"r");
    if (ftopol_rec == NULL) {
        cout << "Can't open file Rec topol: " << topol_file << endl;
        exit(0);
    }

    fgets(buffer, 900, ftopol_rec);

    string p(buffer);

    while(p.substr(0,1) != "["){
        fgets(buffer, 900, ftopol_rec);
         p =  string(buffer);

    }

    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);

    while(p.substr(0,1) != "["){

        fgets(buffer, 900, ftopol_rec);
        p =  string(buffer);

    }
    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);

    int nresidues = 0;

    while(p.substr(0,1) != "["){

    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);

    if(p.substr(0,1) == ";"){
        nresidues++;
        if(nresidues != 1){
        natoms_per_residue.push_back(natoms-1);
        }
        natoms=0;

    }
    if(nresidues > 0 ){
        natoms++;
    }

    }
    natoms_per_residue.push_back(natoms-3);


    rewind (ftopol_rec);

    while(p.substr(0,1) != "["){
        fgets(buffer, 900, ftopol_rec);
        p =  string(buffer);

    }

    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);

    while(p.substr(0,1) != "["){

    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);

    }
    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);

    while(p.substr(0,1) != "["){

    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);


    }

    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);


    vector< string> atomic_name;
    vector< string> residue_name;

    vector< string> atomic_type;

    vector<  int> t_res_index;
    vector<  int> t_a_index;

    vector < double> tcharges;

    int atom_index, res_index;
    char a_name[80], a_type[80], res_name[80] ,tstring[80];
    double charge, tdouble;
    int cnt =0;

    for(int i =0; i < nresidues; i++){

    fgets(buffer, 900, ftopol_rec);
    p =  string(buffer);

           for(int atom =0; atom < natoms_per_residue[i] ;atom++){
            fgets(buffer, 900, ftopol_rec);
            p =  string(buffer);
            sscanf (p.substr(0, p.size() - 2).c_str(),"%d %s %d %s %s %d %lf %s %s %s ",
                 &atom_index,&a_type,&res_index,&res_name,&a_name,&atom_index,&charge,&tstring,&tstring,&tstring,&tstring);
            atomic_name.push_back(a_name);
            atomic_type.push_back(a_type);
            t_res_index.push_back(res_index);

            residue_name.push_back(res_name);
            tcharges.push_back(charge);
            t_a_index.push_back(atom_index);
        }




            this->atom_name_per_residue.push_back(atomic_name);
            this->atom_type_per_residue.push_back(atomic_type);
            this->parameters_charges.push_back(tcharges);
            this->residue_id_per_atom.push_back(t_res_index);
            this->residue_name_per_atom.push_back(residue_name);
            this->atomic_index.push_back(t_a_index);
            atomic_name.clear();
            atomic_type.clear();
            t_res_index.clear();
            residue_name.clear();
            tcharges.clear();
            t_a_index.clear();

    }
    
     fclose(ftopol_rec);
 //sigma epsilon definer

    vector < double> tsigma;
    vector < double> tepsilon;
    for(int i =0; i< this->atom_type_per_residue.size(); i++){
        for(int k =0; k< this->atom_type_per_residue[i].size(); k++) {
             for(int j =0; j < this->atom_names.size(); j++){
                 if(this->atom_type_per_residue[i][k] == this->atom_names[j]){
                     tsigma.push_back(this->atom_sigma[j]);
                     tepsilon.push_back(this->atom_epsilon[j]);

                 }
             }



        }


        this->rec_sigma.push_back(tsigma);
        this->rec_epsilon.push_back(tepsilon);
        tsigma.clear();
        tepsilon.clear();

    }


}

void receptor_pdb::xyz_parser(string pdb_file){
    char buffer[900];


    FILE * fpdb_rec;
    fpdb_rec = fopen (pdb_file.c_str(),"r");
    if (fpdb_rec == NULL) {
        cout << "Can't open file Rec Structure: " << pdb_file  << endl;
        exit(0);
    }
    fgets(buffer, 100, fpdb_rec);

    while(buffer[0] != 'M'){
        fgets(buffer, 100, fpdb_rec);
    }
    int atom_index, res_index;
    char atom_name[40], resname[40];
    double x, y ,z , buf ;
    double x_aver = 0;
    double y_aver = 0;
    double z_aver = 0;
    vector < double> cog;

    vector < double> txyz;
    vector < vector< double>> all_atoms_1resid_xyz;
    vector< vector < vector< double>>> resid_xyz;

    int counter = 0;

        for(int i =0; i< int(this->atom_type_per_residue.size()); i++){
            for(int k =0; k< this->atom_type_per_residue[i].size(); k++) {
                counter++;
                fgets(buffer, 100, fpdb_rec);
                string s(buffer);
                 //sscanf (s.substr(0, s.size() - 2).c_str(),"%s %d %s %s %d %lf %lf %lf %d %d ",
                    //&atom_name,&atom_index,&atom_name,&resname,&res_index,&x,&y,&z,&buf,&buf);
                    x = stof(s.substr(30, 8));
                    y = stof(s.substr(38, 8));
                    z = stof(s.substr(46, 8));
                    atom_index = stoi(s.substr(4, 8));
                    strcpy(atom_name,s.substr(11, 6).c_str());
                    *std::remove_copy_if( s.substr(11, 6).begin(), s.substr(11, 6).end(), atom_name, (int (*)(int))std::isspace ) = 0;
                    strcpy(resname,s.substr(16, 4).c_str());
                    *std::remove_copy_if( s.substr(16, 4).begin(), s.substr(16, 4).end(), resname, (int (*)(int))std::isspace ) = 0;
                    res_index = stoi(s.substr(22, 8));
                    txyz.push_back(x);
                    txyz.push_back(y);
                    txyz.push_back(z);
                    all_atoms_1resid_xyz.push_back(txyz);
                    txyz.clear();
                    x_aver = x_aver +x;
                    y_aver = y_aver +y;
                    z_aver = z_aver +z;

                 //cout << s.substr(0, 4) << "a "  <<  s.substr(4, 8) << "b " << s.substr(11, 6) << "c " << s.substr(16, 4) << "d " << s.substr(20, 2) << "e " <<  s.substr(22, 8) <<
                       //"f " << s.substr(30, 8) << "g " << s.substr(38, 8) << "h " <<  s.substr(46, 8) << endl;



            }
            cog.push_back(x_aver/this->atom_type_per_residue[i].size());
            cog.push_back(y_aver/this->atom_type_per_residue[i].size());
            cog.push_back(z_aver/this->atom_type_per_residue[i].size());
            x_aver=0;
            y_aver=0;
            z_aver=0;
            this->residue_cog.push_back(cog);
            cog.clear();
            this->atom_per_resid_xyz.push_back(all_atoms_1resid_xyz);
            all_atoms_1resid_xyz.clear();
        }
        this->total_atoms = counter;
	fclose(fpdb_rec);

}
