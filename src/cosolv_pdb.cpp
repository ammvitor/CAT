#include "cosolv_pdb.h"


cosolv_pdb::cosolv_pdb()
{

}


void cosolv_pdb::topology_parser(string itp_file){

    char buffer[90];
    FILE * ftopol_rec;

    ftopol_rec = fopen (itp_file.c_str(),"r");

    if (ftopol_rec == NULL) {
        cout << "Can't open file Cosolv itp: " << itp_file  << endl;
        exit(0);
    }
    fgets(buffer, 90, ftopol_rec);

    string p(buffer);

    while(p.substr(0,1) != "["){
        fgets(buffer, 900, ftopol_rec);
         p =  string(buffer);
    }

        char atom_name[80], atom_type[80], temp_string[80];
        double temp_double, sigma, epsilon;
        vector < string> vec_atom_type;
        vector < double> vec_sigma;
        vector < double> vec_epsilon;
        string string_Atom_type;

        fgets(buffer, 900, ftopol_rec);
        fgets(buffer, 900, ftopol_rec);
        p =  string(buffer);

            while( p.size() != 1){


                    if(p.size() != 1){

                        sscanf (p.substr(0, p.size() - 2).c_str(),"%s %s %lf %lf %s %lf %lf %s %s %s ",
                             &atom_name,&atom_type,&temp_double,&temp_double,&temp_string,&sigma,&epsilon,&temp_string,&temp_string,&temp_string);
                        fgets(buffer, 900, ftopol_rec);
                        p =  string(buffer);
                        string_Atom_type = string(atom_type);
                        vec_epsilon.push_back(epsilon);
                        vec_sigma.push_back(sigma);
                        vec_atom_type.push_back(string_Atom_type);
                        }

                     }
            fgets(buffer, 900, ftopol_rec);
            fgets(buffer, 900, ftopol_rec);
            fgets(buffer, 900, ftopol_rec);
            fgets(buffer, 900, ftopol_rec);
            fgets(buffer, 900, ftopol_rec);
            fgets(buffer, 900, ftopol_rec);
            fgets(buffer, 900, ftopol_rec);

            p =  string(buffer);


            int atom_index, res_index;
            char a_name[80], a_type[80], res_name[80] ,tstring[80];
            double charge, mass;
            int cnt =0;


            while( p.size() != 1){


                    if(p.size() != 1){

                        sscanf (p.substr(0, p.size() - 2).c_str(),"%d %s %d %s %s %d %lf %lf %s %s %s ",
                             &atom_index,&a_type,&res_index,&res_name,&a_name,&atom_index,&charge,&mass,&tstring,&tstring,&tstring);
                        fgets(buffer, 900, ftopol_rec);
                        p =  string(buffer);
                        this->atom_names.push_back(a_name);
                        this->atom_types.push_back(a_type);
                        this->charges.push_back(charge);
                        this->atom_masses.push_back(mass);

                        }

        }

            for(int i =0; i< vec_atom_type.size(); i++){
                for(int j =0; j< this->atom_types.size(); j++) {


                    if(vec_atom_type[i] == this->atom_types[j]){
                        //cout << vec_atom_type[i] << " " << this->atom_types[j] << " " << vec_sigma[i] << " " << vec_epsilon[i] << endl;
                        this->atom_sigma.push_back(vec_sigma[i]);
                        this->atom_epsilon.push_back(vec_epsilon[i]);


                }



           }


       }

	fclose(ftopol_rec);
}


void cosolv_pdb::xyz_parser(string pdb_file, int protein_n){
char buffer[900];
FILE * fpdb_rec;
fpdb_rec = fopen (pdb_file.c_str(),"r");
if (fpdb_rec == NULL) {
    cout << "Can't open file Cosolv Rec: " << pdb_file  << endl;
    exit(0);
}
fgets(buffer, 100, fpdb_rec);

    while(buffer[0] != 'M'){
        fgets(buffer, 100, fpdb_rec);
    }
    for(int i=0; i < protein_n+1; i++){
        fgets(buffer, 100, fpdb_rec);

    }

   fgets(buffer, 100, fpdb_rec);

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

         while(buffer[0] != 'T'){


            for(int k =0; k< this->atom_names.size(); k++) {
                string s(buffer);
                if(s.substr(0,3) != "TER"){


                 /*
                 sscanf (s.substr(0, s.size() - 2).c_str(),"%s %d %s %s %d %lf %lf %lf %d %d ",
                    &atom_name,&atom_index,&atom_name,&resname,&res_index,&x,&y,&z,&buf,&buf);
                 //cout << x << " " << y << " " << resname <<  " " << this->atom_names.size() << endl;
                 //cout << s << endl;

                    txyz.push_back(x);
                    txyz.push_back(y);
                    txyz.push_back(z);
                    all_atoms_1resid_xyz.push_back(txyz);
                    txyz.clear();
                    x_aver = x_aver +x;
                    y_aver = y_aver +y;
                    z_aver = z_aver +z;
                    */


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
                    fgets(buffer, 100, fpdb_rec);
                    s = string(buffer);
                 //cout << s.substr(0, 4) << "a "  <<  s.substr(4, 8) << "b " << s.substr(11, 6) << "c " << s.substr(16, 4) << "d " << s.substr(20, 2) << "e " <<  s.substr(22, 8) <<
                        //"f " << s.substr(30, 8) << "g " << s.substr(38, 8) << "h " <<  s.substr(46, 8) << endl;

                }

            }
            cog.push_back(x_aver/this->atom_names.size());
            cog.push_back(y_aver/this->atom_names.size());
            cog.push_back(z_aver/this->atom_names.size());
            x_aver=0;
            y_aver=0;
            z_aver=0;

            this->residue_cog.push_back(cog);
            cog.clear();
            this->atom_per_resid_xyz.push_back(all_atoms_1resid_xyz);
            all_atoms_1resid_xyz.clear();
        }
	fclose(fpdb_rec);
}

