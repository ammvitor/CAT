#include <iostream>
#include<vector>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<time.h>
#include <iostream>
#include <iomanip>
#include<algorithm>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <getopt.h>
#include "receptor_pdb.h"
#include "cosolv_pdb.h"

using namespace std;

/*
( \
 \ \
 / /                |\\
/ /     .-`````-.   / ^`-.
\ \    /         \_/  {|} `o
 \ \  /   .---.   \\ _  ,--'
  \ \/   /     \,  \( `^^^
   \   \/\      (\  )
    \   ) \     ) \ \
jgs  ) /__ \__  ) (\ \___
    (___)))__))(__))(__)))
*/

//@@%%%%%%%%%%%%%%%%%%%%%%%%@@@@@@@@@@@@@@@@@@@@@@@@%%%%%%%%%%%%%%%%%%%%%%%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%@@@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%@@@@@@@@@@@@@@@@@@@@@@@@
//@@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%@@@&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&@@@@@@@@@@@@@@@@@@@@@@@@
//@@%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&@@@&&&&&&&///////////////////////////////////////&&&&&&&&&@@@@@@@@@@@@@@@@@@@@@
//@@************************%%%%%%%%%%%%%%%%%%%%%%%%************************#%%&@@@&&&&////&/(/%##(%&/@/@#&/@/&/(@/@##%/////////////&&&&&&@@@@@@@@@@@@@@@@@@@@@
//@@************************%%%%%%%%%%%%%%%%%%%%%%%%************************#%%&@@@&&&&//////&/%##%#&(@/@##/@/&/(&/@##&//(//////////&&&&&&@@@@@@@@@@@@@@@@@@@@@
//@@***************************************************************************&@@@&&&///////////////////////////////////////////////////&&&@@@@@@@@@@@@@@@@@@@@@
//@@***************************************************************************&@@@&&&//////////////////////////////((((((///////////////&&&@@@@@@@@@@@@@@@@@@@@@
//@@***************************************************************************&@@@&&&//////////////////////////////@@@@@@///////////////&&&@@@@@@@@@@@@@@@@@@@@@
//@@........................************************........................,**&@@@&&&///////(##/////##(/////////@@@/*****@@@////////////&&&@@@@@@******#@@@@@@@@
//@@........................................................./%%%%%%%%%%%,.....%@@@&&&//////////(###(////////////@@@/*****(((%%%/////////&&&@@@#((******#@@@@@@@@
//@@.........................................................(@@@@@@@@@@@*.....%@@@&&&////////((#(/(#(///////////@@@/*******/&@@/////////&&&@@@/********#@@@@@@@@
//@@.........................................................(@@(****#@@@@@..%@@@&&&///////###/////##(/////////@@@/***********@@@@@@@@@@@@(***********#@@@@@@@@
//@@,,,,,,,,,,,,,,,,,,,,,,,,........................,,,,,,,,,(@@%##/**(%%@@@(//&@@@&&&//////////####(////////////@@@/***********%%%%%%%%%%%%/***********#@@@@@@@@
//@@************************........................*********#@@@@@(*****&@@@@@@@@@&&&////////((#(/##(///////////@@@/***********************************#@@@@@@@@
//@@************************************************************&@@@@@&*****&@@@@@@&&&///////(##////(#%///////&@@****************************,,*,,*,,******@@@@@@
//@@************************************************************///&@@@&&(*//(&@@@&&&//////////(###(/////////@@@********  .%&&**********,,***.  &&&*,,***@@@@@@
//@@***************************************************************&@@@@@(****&@@@&&&////////((#(/##(////////@@@********   %@@***************,,,@@@(*****@@@@@@
//@@((((((((((((((((((((((((***********************((((((((((((((((((#@@@@@@@@@@@@&&&///////###/////#%///////@@@*********@@@@@@*******,/@@(*,*@@@(*****@@@@@@
//@@((((((((((((((((((((((((////////////////////////(((((((((((((((((((#####@@@@@@@&&&////////////////////////@@@*********((((((******,*((/**(((,,***,*@@@@@@
//@@((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((&@@@@@@&&&////////////////////////@@@***,,,,,,**************,*************,,,,,@@@@@@
//@@(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((&@@@&&&&&&/////////////////////@@@***,,,,,,**/&@@(****@%*******@@@***,,*,,@@@@@@
//@@(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((&@@@&&&&&&((///////////////////(((###********/&@@%#####@%,*,,,##@@@****,**##@@@@@@
//@@########################((((((((((((((((((((((((########################(((&@@@&&&&&&&/////////////////////@@@/*******/&@@@@@@/,,*,,*,,*,,*,,,,***#@@@@@@@@
//@@########################################################################@@@@@@@@@@&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&@@***********,,,,,*,,***********&@@@@@@@@@@@
//@@#####################################################################%%%&&&&&&&&&&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&@@@((((((((//*,*****/((((((((((((&@@@@@@@@@@@
//@@#####################################################################&@@(********@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@,,,,,,,*@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@########################@@@@@@@@@@@@@@@@@@@@@@@@(****(@@@@@@##%@@@******@@@@@@@@@@@@@@@@@@@@@**,,*,,,,/@@@@******@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&%%%%%&@@@@@@@@@@@@%%%%%%@@@@@@@@@@@@@@@@@@@@@%%/,,*,,@@@@@@@%%%%%%@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%#%%&&@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




int main(int argc, char *argv[])
{

    string input_name, output_name, itp_input, topology_input;
    double cutoff,delta_elec_input,delta_vdw_input,cutoff_clustering;
    int nframes;
    int c;

    delta_elec_input = 0;
    delta_vdw_input = 0;
    // example of inputs:
    //-i /media/jvscunha/My\ Passport/snaps/list_snaps1.dat -c 10 -n 10 -t /media/jvscunha/My\ Passport/snaps/list_snaps1.dat/topol.top -e 0 -s /media/jvscunha/My\ Passport/snaps/list_snaps1.dat/acetaldehyde_GMX.itp -v 0.5 -e 0.5 -o /home/jvscunha/snaps1.txt
    // /home/jvscunha/workspace/Solv_Docking/Solv_docking/Solv_docking -i list.dat -c 10 -t topol.top -e 6 -v 6 -n 1500 -s Br5_GMX.itp -o output.dat -g 5
    while ((c = getopt(argc, argv, "i:o:c:n:t:v:e:s:h:g:")) != -1)



            switch (c){
            case 'i':
                input_name = string(optarg);
                break;
            case 'h':
                printf("Usage %s -i <inputlist> -o <output_distribution> -c <cutoff(angstrom)>  -t <protein_topol_top> -s <cosolv_itp> -e <delta_electrostatic(angstrom)> -v <delta_vdw(angstrom)> \n");
                exit(1);
            case 'o':
                output_name = string(optarg);
                break;
            //case 'n':
              //  nframes = stoi(string(optarg));
                //break;
            case 'c':
                cutoff = stod(string(optarg));
                break;
            case 'v':
                delta_vdw_input = stod(string(optarg));
                break;
            case 'e':
                delta_elec_input = stod(string(optarg));
                break;
            case 's':
                itp_input = string(optarg);
                break;
            case 't':
                topology_input = string(optarg);
                break;
            case 'g':
                cutoff_clustering = stod(string(optarg));
                break;
           }


    //-i /home/jvscunha/cosolvent/Isopropylamine/snaps100/list.log -c 10 -n 99 -t /home/jvscunha/cosolvent/Isopropylamine/topol.top -e 0 -s /home/jvscunha/cosolvent/Isopropylamine/snaps100/imidazole_GMX.itp -v 0.5 -e 0.5 -o /home/jvscunha/cosolvent/Isopropylamine/snaps100.txt

    cout << "******************************************************************" << endl;
    cout << "                            CAT 0.92                                " << endl;
    cout << "                   Cosolvent Analysis  Tool           " << endl;
    cout << "                    Newcastle University                         " << endl;
    cout << "             Cunha, JVS; Sabanes FZ;  Bronowska, AK              " << endl;
    cout << "******************************************************************" << endl;
    //-i /home/jvscunha/cosolvent/Isopropylamine/snaps100/list.log -c 10 -n 80 -t /home/jvscunha/cosolvent/Isopropylamine/topol.top -e 0 -s /home/jvscunha/cosolvent/Isopropylamine/isopropylamine_GMX.itp -v 0 -o /home/jvscunha/cosolvent/Isopropylamine/output.txt





    string pdb_input;
    string topol_input_string = topology_input;

    char buffer[100];
    vector < string> filenames;
    string inpfile = input_name;

    FILE * input;
    input = fopen (inpfile.c_str(),"r");
    if (input == NULL) {
        cout <<  "Can't open list of structures: " << input_name << endl;
        exit(0);
    }
    int initial_residue;

    fgets(buffer, 100, input);
    int nfiles = nframes;
    double vdw_cutoff = 10;
    double d_elec = delta_elec_input;
    double d_vdw = delta_vdw_input;

    vector < double> vdw_final;
    vector < double> elec_final;
    vector < double> count_final;
    vector < double> vdw_final_var;
    vector < double> elec_final_var;
    vector < double> count_final_var;
    vector < vector <vector < double> > >average_structure_perresidue;
    vector <vector < double> > average_residual_xyz;
    vector < double> average_xyz;
    vector < vector <vector < double> > > atomic_spheres_average;
    vector < vector < double> > atomic_per_residue;
    vector < double> atom_sphere_i;


    ifstream counter(input_name.c_str());
    string line;
    int number_of_structures=0;

    if(!counter) {
      cout << "Cannot open input file.\n";
      return 1;
    }

    while(getline(counter, line)) {
          number_of_structures++;
  }

    counter.close();


  nfiles = number_of_structures;

   cout << "File: " << input_name << endl;
   cout << "Output: " << output_name << endl;
   cout << "Number Of Frames: " << nfiles << endl;
   cout << "Cutoff (A) : " << cutoff << endl;
   cout << "Cosolv itp: " << itp_input << endl;
   cout << "Protein topol: " << topology_input << endl;
   cout << "Eletrostatic delta: " << delta_elec_input << endl;
   cout << "Eletrostatic VdW: " << delta_vdw_input << endl;




    //reader for input list
    for(int i =0; i < nfiles; i++ ){
        string s(buffer);
        filenames.push_back(s.substr(0, s.size()-1));
        fgets(buffer, 100, input);

    }



    //reading open file from the list and analyzing it
    int filenames_size = filenames.size();
    //"for" to read the files from the list



    for(int file =0; file < filenames_size; file++){

        receptor_pdb* rec_pdb = new receptor_pdb();
        cosolv_pdb* cov_pdb = new cosolv_pdb();


    pdb_input = filenames[file];

     cout << "Reading file " << file << " " << filenames[file] << '\xd';
    //call methods to read topol (could be read only once
    rec_pdb->topology_parser(topol_input_string);
    rec_pdb->xyz_parser(filenames[file]);
    cov_pdb->topology_parser(itp_input);
    cov_pdb->xyz_parser(filenames[file],rec_pdb->total_atoms-1);

    if(file == 0){

         for(int resid = 0; resid < int(rec_pdb->atom_per_resid_xyz.size()); resid++){
             vdw_final.push_back(0);
             elec_final.push_back(0);
             count_final.push_back(0);
         }
         for(int resid = 0; resid < int(rec_pdb->atom_per_resid_xyz.size()); resid++){
             for(int resid_atom = 0; resid_atom < int(rec_pdb->atom_per_resid_xyz[resid].size()); resid_atom++){
                 average_xyz.push_back(0);
                 average_xyz.push_back(0);
                 average_xyz.push_back(0);

                 average_residual_xyz.push_back(average_xyz);
                 average_xyz.clear();
              }
             atom_sphere_i.push_back(0);
             atom_sphere_i.push_back(0);
             atom_sphere_i.push_back(0);
             atomic_per_residue.push_back(atom_sphere_i);
             atom_sphere_i.clear();
             average_structure_perresidue.push_back(average_residual_xyz);
             average_residual_xyz.clear();
         }
    }


    /*
     * Lennard Jones
     *
     * V_lj = 4sqrt(epsilon_j*epsilon_i)( ( (sigma_j+sigma_i)/(2r_ij))^12 - (sigma_j+sigma_i)/(2r_ij))^6)
     *
     * V_c = 138.935458*q_j*q_i/rij
     *
     *
     *
       */

    vector < double> vdw_per_residue;
    vector < double> elec_per_residue;
    vector < double> count_per_residue;

    int count =0;
    int n_inside;
    double dr , dcog, dx, dy, dz ,epsilon_ij,sigma_ij, tvdw, telec;
    double vdw = 0;
    double elec = 0;





    for(int resid = 0; resid < int(rec_pdb->atom_per_resid_xyz.size()); resid++){
        vdw = 0;
        elec = 0;
        count = 0;
        n_inside=0;
        //protein average structure calculator
        for(int resid_atom = 0; resid_atom < int(rec_pdb->atom_per_resid_xyz[resid].size()); resid_atom++){
            average_structure_perresidue[resid][resid_atom][0] = average_structure_perresidue[resid][resid_atom][0] +rec_pdb->atom_per_resid_xyz[resid][resid_atom][0];
            average_structure_perresidue[resid][resid_atom][1] = average_structure_perresidue[resid][resid_atom][1] +rec_pdb->atom_per_resid_xyz[resid][resid_atom][1];
            average_structure_perresidue[resid][resid_atom][2] = average_structure_perresidue[resid][resid_atom][2] +rec_pdb->atom_per_resid_xyz[resid][resid_atom][2];
         }
        
        
        //Calculate distances between COGs (centers of geometry)
        for(int cosolv_id =0; cosolv_id < int(cov_pdb->atom_per_resid_xyz.size()); cosolv_id++){
            dcog = sqrt((rec_pdb->residue_cog[resid][0]-cov_pdb->residue_cog[cosolv_id][0])*(rec_pdb->residue_cog[resid][0]-cov_pdb->residue_cog[cosolv_id][0])+
                        (rec_pdb->residue_cog[resid][1]-cov_pdb->residue_cog[cosolv_id][1])*(rec_pdb->residue_cog[resid][1]-cov_pdb->residue_cog[cosolv_id][1])+
                        (rec_pdb->residue_cog[resid][2]-cov_pdb->residue_cog[cosolv_id][2])*(rec_pdb->residue_cog[resid][2]-cov_pdb->residue_cog[cosolv_id][2]));
            //cutoff = sphere radius assigned in the beginning
            if(dcog < cutoff){

                    atomic_per_residue[resid][0] = atomic_per_residue[resid][0] + cov_pdb->residue_cog[cosolv_id][0];
                    atomic_per_residue[resid][1] = atomic_per_residue[resid][1] + cov_pdb->residue_cog[cosolv_id][1];
                    atomic_per_residue[resid][2] = atomic_per_residue[resid][2] + cov_pdb->residue_cog[cosolv_id][2];
                    n_inside++;
                    count++;
                    
                    //interaction calculator using the hybrid score function
                    for(int resid_atom = 0; resid_atom < int(rec_pdb->atom_per_resid_xyz[resid].size()); resid_atom++){

                        for(int cosolv_atom = 0; cosolv_atom < int(cov_pdb->atom_per_resid_xyz[cosolv_id].size()); cosolv_atom++){

                            dx = (rec_pdb->atom_per_resid_xyz[resid][resid_atom][0]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][0])*(rec_pdb->atom_per_resid_xyz[resid][resid_atom][0]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][0]);
                            dy = (rec_pdb->atom_per_resid_xyz[resid][resid_atom][1]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][1])*(rec_pdb->atom_per_resid_xyz[resid][resid_atom][1]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][1]);
                            dz = (rec_pdb->atom_per_resid_xyz[resid][resid_atom][2]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][2])*(rec_pdb->atom_per_resid_xyz[resid][resid_atom][2]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][2]);
                            dr = (sqrt(dx + dy + dz)/10);
                            epsilon_ij = sqrt(rec_pdb->rec_epsilon[resid][resid_atom]*cov_pdb->atom_epsilon[cosolv_atom]);
                            sigma_ij = (rec_pdb->rec_epsilon[resid][resid_atom]+cov_pdb->atom_epsilon[cosolv_atom])/2;
                            tvdw = 4*epsilon_ij*(pow((sigma_ij/(dr+d_vdw)),12) - pow((sigma_ij/(dr+d_vdw)),6));
                            telec = (138.935458*rec_pdb->parameters_charges[resid][resid_atom]*cov_pdb->charges[cosolv_atom])/(dr+d_elec);
                            
                            //cuttoff for extreme atomic classes
                            if(vdw_cutoff > vdw){
                                vdw = vdw + tvdw;
                                elec = elec + telec;

                            }

                        }

                    }


              }


        }

        if(n_inside > 0){
            atomic_per_residue[resid][0] = atomic_per_residue[resid][0]/n_inside;
            atomic_per_residue[resid][1] = atomic_per_residue[resid][1]/n_inside;
            atomic_per_residue[resid][2] = atomic_per_residue[resid][2]/n_inside;
        }



        vdw_per_residue.push_back(vdw);
        elec_per_residue.push_back(elec);
        count_per_residue.push_back(count);

    }

    //average vdw elec and count colculator
    for(int i =0; i < int(vdw_per_residue.size()); i++){
        vdw_final[i] = vdw_final[i] + vdw_per_residue[i]/nfiles;
        elec_final[i] = elec_final[i] + elec_per_residue[i]/nfiles;
        count_final[i] = count_final[i] + count_per_residue[i]/nfiles;
    }

    
    
    //average protein structure printer

    if(file == int(filenames.size()-1)){
        FILE * output_average;
        output_average = fopen ("aver_withB.pdb","w");

        int atom_id=0;
        for(int resid = 0; resid < int(rec_pdb->atom_per_resid_xyz.size()); resid++){
            for(int resid_atom = 0; resid_atom < int(rec_pdb->atom_per_resid_xyz[resid].size()); resid_atom++){
                atom_id++;
                if(count_final[resid] >= 1){
                fprintf(output_average, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                         "ATOM",
                         atom_id,
                         rec_pdb->atom_name_per_residue[resid][resid_atom].c_str(),
                         " ",
                         rec_pdb->residue_name_per_atom[resid][resid_atom].c_str(),
                         "A",
                         resid+rec_pdb->residue_id_per_atom[0][0],
                         " ",
                         average_structure_perresidue[resid][resid_atom][0]/nfiles,
                         average_structure_perresidue[resid][resid_atom][1]/nfiles,
                         average_structure_perresidue[resid][resid_atom][2]/nfiles,
                         1.0,
                         (elec_final[resid]+vdw_final[resid])/count_final[resid]);
                  }
                
                
                if(count_final[resid] < 1){
                    fprintf(output_average, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                             "ATOM",
                             atom_id,
                             rec_pdb->atom_name_per_residue[resid][resid_atom].c_str(),
                             " ",
                             rec_pdb->residue_name_per_atom[resid][resid_atom].c_str(),
                             "A",
                             resid+rec_pdb->residue_id_per_atom[0][0],
                             " ",
                             average_structure_perresidue[resid][resid_atom][0]/nfiles,
                             average_structure_perresidue[resid][resid_atom][1]/nfiles,
                             average_structure_perresidue[resid][resid_atom][2]/nfiles,
                             1.0,
                             0.0);
                      }
             }

        }
            fclose(output_average);
    }
    //printed_average



    vdw_per_residue.clear();
    elec_per_residue.clear();
    count_per_residue.clear();
    atomic_spheres_average.push_back(atomic_per_residue);
    atomic_per_residue.clear();


    for(int resid = 0; resid < int(rec_pdb->atom_per_resid_xyz.size()); resid++){
        atom_sphere_i.push_back(0);
        atom_sphere_i.push_back(0);
        atom_sphere_i.push_back(0);
        atomic_per_residue.push_back(atom_sphere_i);
        atom_sphere_i.clear();
    }


    delete rec_pdb;
    delete cov_pdb;

   }



   cout << "Averages calculated" << endl;

   //First Cycle of calculations done - to reduce ram memory usage - since the traj size is huge, it rereads the traj to calculate the variance


   cout << "Calculating Variances" << endl;


   //variance calculator
    for(int file =0; file < int(filenames.size()); file++){

        receptor_pdb* rec_pdb = new receptor_pdb();
        cosolv_pdb* cov_pdb = new cosolv_pdb();

        cout << "Reading file " << file << '\xd';

    pdb_input = filenames[file];


    rec_pdb->topology_parser(topol_input_string);
    rec_pdb->xyz_parser(filenames[file]);

    cov_pdb->topology_parser(itp_input);
    initial_residue= rec_pdb->residue_id_per_atom[0][0];

    cov_pdb->xyz_parser(filenames[file],rec_pdb->total_atoms-1);
    if(file == 0){
         for(int resid = 0; resid < int(rec_pdb->atom_per_resid_xyz.size()); resid++){
             vdw_final_var.push_back(0);
             elec_final_var.push_back(0);
             count_final_var.push_back(0);
         }
    }


    /*
     * Lennard Jones
     *
     * V_lj = 4sqrt(epsilon_j*epsilon_i)( ( (sigma_j+sigma_i)/(2r_ij))^12 - (sigma_j+sigma_i)/(2r_ij))^6)
     *
     * V_c = 138.935458*q_j*q_i/rij
     *
     *
     *
       */

    vector < double> vdw_per_residue;
    vector < double> elec_per_residue;
    vector < double> count_per_residue;


    int count =0;
    double dr , dcog, dx, dy, dz ,epsilon_ij,sigma_ij, tvdw, telec;
    double vdw = 0;
    double elec = 0;
    for(int resid = 0; resid < int(rec_pdb->atom_per_resid_xyz.size()); resid++){
        vdw = 0;
        elec = 0;
        count = 0;
        for(int cosolv_id =0; cosolv_id < int(cov_pdb->atom_per_resid_xyz.size()); cosolv_id++){
            dcog = sqrt((rec_pdb->residue_cog[resid][0]-cov_pdb->residue_cog[cosolv_id][0])*(rec_pdb->residue_cog[resid][0]-cov_pdb->residue_cog[cosolv_id][0])+
                        (rec_pdb->residue_cog[resid][1]-cov_pdb->residue_cog[cosolv_id][1])*(rec_pdb->residue_cog[resid][1]-cov_pdb->residue_cog[cosolv_id][1])+
                        (rec_pdb->residue_cog[resid][2]-cov_pdb->residue_cog[cosolv_id][2])*(rec_pdb->residue_cog[resid][2]-cov_pdb->residue_cog[cosolv_id][2]));

            if(dcog < cutoff){
                    count++;
                    for(int resid_atom = 0; resid_atom < int(rec_pdb->atom_per_resid_xyz[resid].size()); resid_atom++){
                        for(int cosolv_atom = 0; cosolv_atom < int(cov_pdb->atom_per_resid_xyz[cosolv_id].size()); cosolv_atom++){

                            dx = (rec_pdb->atom_per_resid_xyz[resid][resid_atom][0]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][0])*(rec_pdb->atom_per_resid_xyz[resid][resid_atom][0]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][0]);
                            dy = (rec_pdb->atom_per_resid_xyz[resid][resid_atom][1]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][1])*(rec_pdb->atom_per_resid_xyz[resid][resid_atom][1]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][1]);
                            dz = (rec_pdb->atom_per_resid_xyz[resid][resid_atom][2]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][2])*(rec_pdb->atom_per_resid_xyz[resid][resid_atom][2]-cov_pdb->atom_per_resid_xyz[cosolv_id][cosolv_atom][2]);
                            dr = (sqrt(dx + dy + dz)/10);
                            epsilon_ij = sqrt(rec_pdb->rec_epsilon[resid][resid_atom]*cov_pdb->atom_epsilon[cosolv_atom]);
                            sigma_ij = (rec_pdb->rec_epsilon[resid][resid_atom]+cov_pdb->atom_epsilon[cosolv_atom])/2;
                            tvdw = 4*epsilon_ij*(pow((sigma_ij/(dr+d_vdw)),12) - pow((sigma_ij/(dr+d_vdw)),6));
                            telec = (138.935458*rec_pdb->parameters_charges[resid][resid_atom]*cov_pdb->charges[cosolv_atom])/(dr+d_elec);
                            elec = elec + telec;
                            if(vdw_cutoff > vdw){
                                vdw = vdw + tvdw;
                            }

                        }

                    }




              }


        }

        vdw_per_residue.push_back(vdw);
        elec_per_residue.push_back(elec);
        count_per_residue.push_back(count);
     }

    for(int i =0; i < int(vdw_per_residue.size()); i++){

        vdw_final_var[i] =  vdw_final_var[i] + ((vdw_final[i] - vdw_per_residue[i])*(vdw_final[i] - vdw_per_residue[i]))/nfiles;
        elec_final_var[i] = elec_final_var[i] + ((elec_final[i] - elec_per_residue[i])*(elec_final[i] - elec_per_residue[i]))/nfiles;
        count_final_var[i] =count_final_var[i] + ((count_final[i] - count_per_residue[i])*(count_final[i] - count_per_residue[i]))/nfiles;



    }

    vdw_per_residue.clear();
    elec_per_residue.clear();
    count_per_residue.clear();
    delete rec_pdb;
    delete cov_pdb;

   }



    for(int i =0; i < int(vdw_final_var.size()); i++){
        vdw_final_var[i] =  sqrt(vdw_final_var[i]);
        elec_final_var[i] = sqrt(elec_final_var[i]);
        count_final_var[i] = sqrt(count_final_var[i]);


    }

   //printer of frag score count vdw and elec

    FILE * output2;
    output2 = fopen (output_name.c_str(),"w");

    fprintf (output2, "%10.10s %10.5s %10.5s %10.5s %10.5s %10.5s %10.5s %10.5s %10.5s \n","resid"," ele ","e_VAR"," vdw ","v_var","count","c_var","vol_c","vol_");


    for(int i =0; i < int(elec_final.size()); i++){
        if(count_final[i] != 0){
       fprintf (output2, "%10.5d %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n",i+initial_residue,elec_final[i]/count_final[i],elec_final_var[i]/count_final[i],vdw_final[i]/count_final[i],vdw_final_var[i]/count_final[i],count_final[i],count_final_var[i],(cutoff*cutoff*cutoff)/count_final[i],((cutoff*cutoff*cutoff)/count_final[i])-((cutoff*cutoff*cutoff)/count_final_var[i]));
        }
        if(count_final[i] == 0){
       fprintf (output2, "%10.5d %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n",i+initial_residue,0,0,0,0,count_final[i],count_final_var[i],0,0);
        }
    }


    cout << "Destribution calculated:" << output_name << '\xd';

    vector < vector < double> > final_ball_residue;
    vector < double > xyz_ball;
    xyz_ball.push_back(0);
    xyz_ball.push_back(0);
    xyz_ball.push_back(0);

    for(int frame = 0; frame < int(atomic_spheres_average[0].size()); frame++){

                final_ball_residue.push_back(xyz_ball);

      }


    vector < int> normalizer;
    for(int resid = 0; resid < int(atomic_spheres_average[0].size()); resid++){
        normalizer.push_back(0);
    }


    for(int frame = 0; frame < int(atomic_spheres_average.size()); frame++){
        for(int resid = 0; resid < int(atomic_spheres_average[frame].size()); resid++){
            if(atomic_spheres_average[frame][resid][0] != 0){
                final_ball_residue[resid][0] =  final_ball_residue[resid][0] + atomic_spheres_average[frame][resid][0];
                final_ball_residue[resid][1] =  final_ball_residue[resid][1] + atomic_spheres_average[frame][resid][1];
                final_ball_residue[resid][2] =  final_ball_residue[resid][2] + atomic_spheres_average[frame][resid][2];
                normalizer[resid]++;

            }

        }
      }

    for(int resid = 0; resid < int(atomic_spheres_average[0].size()); resid++){
        if(normalizer[resid] != 0){
            final_ball_residue[resid][0] =  final_ball_residue[resid][0]/normalizer[resid];
            final_ball_residue[resid][1] =  final_ball_residue[resid][1]/normalizer[resid];
            final_ball_residue[resid][2] =  final_ball_residue[resid][2]/normalizer[resid];

        }

    }



    FILE * output_balls;
    output_balls = fopen ("all_spheres_count_variance.pdb","w");
    vector < vector < double> > sphere_bkp;

    for(int i =0; i < int(final_ball_residue.size()); i++){
        if(int(final_ball_residue[i][0]) != 0){
        fprintf(output_balls, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                 "ATOM",
                 i+1+initial_residue,
                 "C",
                 " ",
                 "C",
                 "A",
                 i,
                 " ",
                 final_ball_residue[i][0],
                 final_ball_residue[i][1],
                 final_ball_residue[i][2],
                 1.0,
                 count_final_var[i]);
        fprintf(output_balls, "TER\n");

            }

        sphere_bkp.push_back(final_ball_residue[i]);
    }





    double dcog;
    double x_cluster_temp =0;
    double y_cluster_temp =0;
    double z_cluster_temp =0;
    int cluster_int = 0;

    double energy_per_ball =0;
    double variance_per_ball =0;
    double occupancy_per_ball =0;

    vector < vector < double > > cluster_end;
    vector < double > cluster_i;
    vector < double > cluster_size;
    vector < double > energy_per_cluster;
    vector < double > occupancy_perball_vec;
    vector < double > variance_occupancy_perball_vec;


    for(int i =0; i < int(final_ball_residue.size()); i++){
        x_cluster_temp = final_ball_residue[i][0];
        y_cluster_temp = final_ball_residue[i][1];
        z_cluster_temp = final_ball_residue[i][2];
        energy_per_ball = elec_final[i]+vdw_final[i];
        for(int k =i+1; k < int(final_ball_residue.size()); k++){

            dcog=   (final_ball_residue[i][0]-final_ball_residue[k][0])*(final_ball_residue[i][0]-final_ball_residue[k][0])+
                    (final_ball_residue[i][1]-final_ball_residue[k][1])*(final_ball_residue[i][1]-final_ball_residue[k][1])+
                    (final_ball_residue[i][2]-final_ball_residue[k][2])*(final_ball_residue[i][2]-final_ball_residue[k][2]);
            dcog =sqrt(dcog);

        if(dcog < cutoff_clustering && dcog != 0){

            x_cluster_temp = x_cluster_temp + final_ball_residue[k][0];
            y_cluster_temp = y_cluster_temp + final_ball_residue[k][1];
            z_cluster_temp = z_cluster_temp + final_ball_residue[k][2];
            //final_ball_residue.erase(final_ball_residue.begin()+k+1);
            final_ball_residue[k][0] = 9999;
            final_ball_residue[k][1] = 9999;
            final_ball_residue[k][2] = 9999;
            energy_per_ball = energy_per_ball + elec_final[k]+vdw_final[k];
            variance_per_ball = variance_per_ball + count_final_var[k];
            occupancy_per_ball = occupancy_per_ball + count_final[k];
            cluster_int++;
            }
          }

            if(cluster_int > 0){
                x_cluster_temp = x_cluster_temp/(cluster_int+1);
                y_cluster_temp = y_cluster_temp/(cluster_int+1);
                z_cluster_temp = z_cluster_temp/(cluster_int+1);
                cluster_i.push_back(x_cluster_temp);
                cluster_i.push_back(y_cluster_temp);
                cluster_i.push_back(z_cluster_temp);
                cluster_end.push_back(cluster_i);
                cluster_i.clear();
                  x_cluster_temp =0;
                  y_cluster_temp =0;
                  z_cluster_temp =0;
                  energy_per_cluster.push_back(energy_per_ball/(cluster_int+1));
                  occupancy_perball_vec.push_back(occupancy_per_ball/(cluster_int+1));
                  variance_occupancy_perball_vec.push_back(variance_per_ball/(cluster_int+1));
                  cluster_size.push_back(cluster_int);
                  variance_per_ball=0;
                  occupancy_per_ball=0;
                  cluster_int=0;
                  energy_per_ball =0;

              }


    }

    int contact_number = 0;
    int max_contact_number = 0;

    vector < vector < double > > sorted_clusters;
    vector < double > t_sorted_clusters;
    vector < double> n_contact_vec;





    for(int i =0; i < int(cluster_end.size()); i++){
        for(int resid = 0; resid < int(average_structure_perresidue.size()); resid++){
                    for(int resid_atom = 0; resid_atom < int(average_structure_perresidue[resid].size()); resid_atom++){
                    dcog = (average_structure_perresidue[resid][resid_atom][0]/nfiles-cluster_end[i][0])*(average_structure_perresidue[resid][resid_atom][0]/nfiles-cluster_end[i][0])+
                           (average_structure_perresidue[resid][resid_atom][1]/nfiles-cluster_end[i][1])*(average_structure_perresidue[resid][resid_atom][1]/nfiles-cluster_end[i][1])+
                           (average_structure_perresidue[resid][resid_atom][2]/nfiles-cluster_end[i][2])*(average_structure_perresidue[resid][resid_atom][2]/nfiles-cluster_end[i][2]);
                        if(sqrt(dcog) < cutoff){
                             contact_number++;
                        }

                    }
                    if(max_contact_number < contact_number){
                        max_contact_number = contact_number;
                    }
        }
        n_contact_vec.push_back(contact_number);
        contact_number = 0;
    }

    for(int i =0; i < int(cluster_end.size()); i++){
        if(cluster_end[i][0] != 0){
                t_sorted_clusters.push_back(cluster_end[i][0]);
                t_sorted_clusters.push_back(cluster_end[i][1]);
                t_sorted_clusters.push_back(cluster_end[i][2]);
                t_sorted_clusters.push_back((energy_per_cluster[i])*(n_contact_vec[i]/max_contact_number)*(1-variance_occupancy_perball_vec[i]/(variance_occupancy_perball_vec[i]+occupancy_perball_vec[i])));
                sorted_clusters.push_back(t_sorted_clusters);
                t_sorted_clusters.clear();
            }

    }

    std::sort(sorted_clusters.begin(), sorted_clusters.end(),[](const std::vector<double>& a, const std::vector<double>& b) {
      return a[3] < b[3];
    });

    //print energy
    FILE * output_cluster_top;
    output_cluster_top = fopen ("all_clusters_top.pdb","w");

    if(cluster_end.size()>10){
        for(int i =0; i < 10; i++){
            if(cluster_end[i][0] != 0){
                fprintf(output_cluster_top, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                         "ATOM",
                         i+initial_residue,
                         "C",
                         " ",
                         "C",
                         "A",
                         i+initial_residue,
                         " ",
                         sorted_clusters[i][0],
                         sorted_clusters[i][1],
                         sorted_clusters[i][2],
                         1.0,
                         sorted_clusters[i][3]);
                fprintf(output_cluster_top, "TER\n");
                }

        }
        }
    else{
        for(int i =0; i < int(cluster_end.size()); i++){
            if(cluster_end[i][0] != 0){
                fprintf(output_cluster_top, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                         "ATOM",
                         i+initial_residue,
                         "C",
                         " ",
                         "C",
                         "A",
                         i,
                         " ",
                         cluster_end[i][0],
                         cluster_end[i][1],
                         cluster_end[i][2],
                         1.0,
                         n_contact_vec[i]/max_contact_number);
                fprintf(output_cluster_top, "TER\n");
                }

        }
    }

    FILE * output_cluster;
    output_cluster = fopen ("all_clusters.pdb","w");


    for(int i =0; i < int(cluster_end.size()); i++){
        if(cluster_end[i][0] != 0){
            fprintf(output_cluster, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                     "ATOM",
                     i+initial_residue,
                     "C",
                     " ",
                     "C",
                     "A",
                     i,
                     " ",
                     cluster_end[i][0],
                     cluster_end[i][1],
                     cluster_end[i][2],
                     1.0,
                     n_contact_vec[i]/max_contact_number);
            fprintf(output_cluster, "TER\n");
            }

    }

    //print energy
    FILE * output_energy;
    output_energy = fopen ("all_energy.pdb","w");

    for(int i =0; i < int(cluster_end.size()); i++){
        if(cluster_end[i][0] != 0){
            fprintf(output_energy, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                     "ATOM",
                     i+initial_residue,
                     "C",
                     " ",
                     "C",
                     "A",
                     i,
                     " ",
                     cluster_end[i][0],
                     cluster_end[i][1],
                     cluster_end[i][2],
                     1.0,
                     n_contact_vec[i]/max_contact_number);
            fprintf(output_energy, "TER\n");
            }

    }



    //print contact
    FILE * output_contact;
    output_contact = fopen ("all_contact.pdb","w");

    for(int i =0; i < int(cluster_end.size()); i++){
        if(cluster_end[i][0] != 0){
            fprintf(output_contact, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                     "ATOM",
                     i+initial_residue,
                     "C",
                     " ",
                     "C",
                     "A",
                     i,
                     " ",
                     cluster_end[i][0],
                     cluster_end[i][1],
                     cluster_end[i][2],
                     1.0,
                     n_contact_vec[i]/max_contact_number);
            fprintf(output_contact, "TER\n");
            }

    }

    //print energy
    FILE * output_cluster_spheres_regions;
    output_cluster_spheres_regions = fopen ("all_clusters_top_spheresregions.pdb","w");
    long double sphere_distance;
    if(cluster_end.size()>10){
        for(int i =0; i < 10; i++){
            if(cluster_end[i][0] != 0){
                for(int k =0; k < int(sphere_bkp.size()); k++){
                    if(sphere_bkp[i][0] != 0){

                        sphere_distance = (((sorted_clusters[i][0]-sphere_bkp[k][0])*(sorted_clusters[i][0]-sphere_bkp[k][0]))+
                                           ((sorted_clusters[i][1]-sphere_bkp[k][1])*(sorted_clusters[i][1]-sphere_bkp[k][1]))+
                                           ((sorted_clusters[i][2]-sphere_bkp[k][2])*(sorted_clusters[i][2]-sphere_bkp[k][2])));


                        if(sqrt(sphere_distance) < cutoff_clustering){

                         fprintf(output_cluster_spheres_regions, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                         "ATOM",
                         i+initial_residue,
                         "C",
                         " ",
                         "C",
                         "A",
                         k,
                         " ",
                         sphere_bkp[k][0],
                         sphere_bkp[k][1],
                         sphere_bkp[k][2],
                         1.0,
                         sorted_clusters[i][3]);
                         //cout << k << " " << i << " " << sphere_bkp[k][0] << endl;

                        }
                    }
                    }
                }
            fprintf(output_cluster_spheres_regions, "TER\n");


        }
        }



    else{
        for(int i =0; i < int(cluster_end.size()); i++){
            if(cluster_end[i][0] != 0){
                for(int k =0; k < int(final_ball_residue.size()); k++){
                        sphere_distance = (sqrt((sorted_clusters[i][0]-final_ball_residue[k][0])*(sorted_clusters[i][0]-final_ball_residue[k][0])+
                                                (sorted_clusters[i][1]-final_ball_residue[k][1])*(sorted_clusters[i][1]-final_ball_residue[k][1])+
                                                (sorted_clusters[i][2]-final_ball_residue[k][2])*(sorted_clusters[i][2]-final_ball_residue[k][2])));
                        if(sphere_distance<cutoff_clustering){
                           // cout << i << " " << sphere_distance << endl;

                         fprintf(output_cluster_spheres_regions, "%4s%7d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                         "ATOM",
                         i+initial_residue,
                         "C",
                         " ",
                         "C",
                         "A",
                         i,
                         " ",
                         final_ball_residue[k][0],
                         final_ball_residue[k][1],
                         final_ball_residue[k][2],
                         1.0,
                         sorted_clusters[i][3]);
                        cout << i << " " << sphere_distance << endl;

                        }
                    }
                }
            fprintf(output_cluster_spheres_regions, "TER\n");


        }

        }




cout << "Thanks for using CAT" << endl;

}

