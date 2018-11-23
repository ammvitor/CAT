#ifndef CLUSTERER_H
#define CLUSTERER_H
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

using namespace std;

class clusterer
{
public:

    clusterer(string input_2dproj_name, int n_centroids);
    string inputfile_set;
    vector< vector <double> > dproj;
    int total_clusters;
    vector < vector < double > > centers;
    void define_starting_centers(int number_of_centroids);
    vector < int > knn_class;
    vector <double> occupancy;
    vector < vector < double > > distances_to_Centers();
    bool center_comparator(vector< vector <double > > center, vector< vector <double > > prev_center);
    void classifier_knn();
    void center_redefinition();
    void print_clusters(string outname);
    void print_divided_traj(string outname);



};

#endif // CLUSTERER_H
