#ifndef COSOLV_PDB_H
#define COSOLV_PDB_H
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
#include<algorithm>
#include <stdio.h>
#include <string>
#include <getopt.h>
#include <stdio.h>

using namespace std;


class cosolv_pdb
{
public:
    cosolv_pdb();
    vector< string> atom_names;
    vector< string> atom_types;
    vector< double> atom_masses;
    vector< double> atom_sigma;
    vector< double> atom_epsilon;
    void topology_parser(string itp_file);

     vector< double> sigma;
     vector< double> charges;
     vector< double> epsilon;
     void xyz_parser(string pdb_file, int protein_n);
     vector < vector < vector< double > > > atom_per_resid_xyz;
     vector <vector< double > > residue_cog;


};

#endif // COSOLV_PDB_H
