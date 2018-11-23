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
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "clusterer.h"

using namespace std;



int main(int argc, char *argv[])
{



    string inputfile_set, output_name,traj_name;
    int nclusters;
    int c;
    traj_name = "";
    bool printtraj = false;
    while ((c = getopt(argc, argv, "i:n:h:o:t:")) != -1)



            switch (c){
            case 'i':
                inputfile_set = string(optarg);
                break;
            case 'h':
                printf("Usage %s -i <2dproj_xvg> -o <output_prefix> -n n_clusters -t trjname \n");
                exit(1);
            case 'o':
                output_name = string(optarg);
                break;
           case 'n':
                nclusters = stoi(string(optarg));
                break;
           case 't':
                traj_name = string(optarg);
                break;
             }



    if(traj_name != ""){
        printtraj = true;
    }

    clusterer* PCA_object = new clusterer(inputfile_set,nclusters);


    vector < vector < double > > centers_prev = PCA_object->centers;

    PCA_object->classifier_knn();
    PCA_object->center_redefinition();

    while(PCA_object->center_comparator(PCA_object->centers,centers_prev)){
        centers_prev = PCA_object->centers;
        PCA_object->classifier_knn();
        PCA_object->center_redefinition();
    }

    PCA_object->print_clusters(output_name);

    if(printtraj){
    PCA_object->print_divided_traj(traj_name);
    }
}

