#ifndef RECEPTOR_PDB_H
#define RECEPTOR_PDB_H

#include <iostream>
#include<vector>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<cstring>
#include<time.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <getopt.h>
#include <stdio.h>

using namespace std;


class receptor_pdb
{
public:
    receptor_pdb();
    vector< string> atom_names;
    vector< double> atom_masses;
    vector< double> atom_sigma;
    vector< double> atom_epsilon;
    void topology_parser(string topol_file);

    vector < vector< double> > rec_sigma;
    vector < vector < string> > atom_name_per_residue;
    vector < vector < string> > atom_type_per_residue;
    vector < vector< double> > parameters_charges;
    vector < vector< int> > residue_id_per_atom;
    vector < vector< string> > residue_name_per_atom;
    vector < vector< int> > atomic_index;
     vector < vector < vector< double> > > atom_per_resid_xyz;
    vector < vector< double > > rec_epsilon;
    void xyz_parser(string pdb_file);
    int total_atoms;
    vector <vector< double > > residue_cog;


};

#endif // RECEPTOR_PDB_H
