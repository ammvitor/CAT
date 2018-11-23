#include "clusterer.h"

clusterer::clusterer(string input_2dproj_name, int n_centroids)
{

        this->inputfile_set = input_2dproj_name;
        ifstream counter(inputfile_set.c_str());
        string line;
        vector < double> vec_temp;
        double X1, X2;
        if(!counter) {
          cout << "Cannot open input file.\n";
           exit (EXIT_FAILURE);
        }
        string a  ;
        while(getline(counter, line)) {
            a = line.at(0);
           if ( a != "#" && a != "@" ){
           //cout << line << endl;
            stringstream ss(line.c_str());
            ss >>  X1 >> X2 ;
            vec_temp.push_back(X1);
            vec_temp.push_back(X2);
            this->dproj.push_back(vec_temp);
            vec_temp.clear();
            }
        }
        counter.close();

        this->define_starting_centers(n_centroids);
}



bool clusterer::center_comparator(vector< vector <double > > center, vector< vector <double > > prev_center){
    int score =0;
    for(int i  = 0; i < center.size() ; i++){
        if(center[i][0] == prev_center[i][0]){
            score++;
        }
        if(center[i][1] == prev_center[i][1]){
            score++;
        }
    }
    if(score == center.size()*2){
        return false;

    }
}


void clusterer::define_starting_centers(int number_of_centroids){
    this->total_clusters = number_of_centroids;
    vector < double >temp_xyz;

        for(int i = 0; i< this->total_clusters; i++){

                temp_xyz.push_back(this->dproj[double(this->dproj.size())*(double(i)/double(this->total_clusters))][0]);
                temp_xyz.push_back(this->dproj[double(dproj.size())*(double(i)/double(this->total_clusters))][1]);
                cout << "Starting center coordinates for c" << i+1 << " is " << this->dproj[double(this->dproj.size())*(double(i)/double(this->total_clusters))][0] << " " << this->dproj[double(this->dproj.size())*(double(i)/double(this->total_clusters))][1] << endl;
                this->centers.push_back(temp_xyz);
                temp_xyz.clear();
        }

}



vector < vector < double > > clusterer::distances_to_Centers(){

    vector < vector < double > > distances_to_Cn;
    vector < double > index_distance_to_C;



    for(int k=0; k < this->centers.size(); k++){
        for(int i=0; i < this->dproj.size(); i++){
                index_distance_to_C.push_back(sqrt((this->dproj[i][0]-this->centers[k][0])*(dproj[i][0]-this->centers[k][0])+(this->dproj[i][1]-this->centers[k][1])*(this->dproj[i][1]-this->centers[k][1])));
        }
        distances_to_Cn.push_back(index_distance_to_C);
        index_distance_to_C.clear();
    }
    return distances_to_Cn;
}


void clusterer::classifier_knn(){

       vector < vector < double > > distances_to_Cn = this->distances_to_Centers();
      //define to which group it is
      vector <int > knn_class(distances_to_Cn[0].size());


      double min_distance = 9999;
      double belong_to_group;

      for(int i=0; i < distances_to_Cn[0].size(); i++){

          for(int k=0; k < distances_to_Cn.size(); k++){

           if(distances_to_Cn[k][i] < min_distance ){
               min_distance = distances_to_Cn[k][i];
               belong_to_group = k;
           }



          }

          knn_class[i] = belong_to_group;

          min_distance = 9999;

      }
      this->knn_class = knn_class;


}

void clusterer::center_redefinition(){


    vector < double > occup_local(this->centers.size());

    std::fill (occup_local.begin(),occup_local.begin()+this->centers.size(),1);
    this->occupancy = occup_local;
    vector < vector < double > > distances_to_Cn = this->distances_to_Centers();

    for(int k=0; k < distances_to_Cn.size(); k++){

      for(int i=0; i < distances_to_Cn[0].size(); i++){


         if(this->knn_class[i] == k ){
             this->centers[k][0] =  this->centers[k][0] + this->dproj[i][0];
             this->centers[k][1] =  this->centers[k][1] + this->dproj[i][1];
             this->occupancy[k]++;

         }

        }


    }
    for(int k=0; k < this->centers.size(); k++){
        this->centers[k][0] =  this->centers[k][0]/this->occupancy[k];
        this->centers[k][1] = this->centers[k][1]/this->occupancy[k];
        //cout << "centers " << this->centers[k][0] << " " << this->centers[k][1] << endl;
    }

  }


void clusterer::print_clusters(string outname){

    for(int blob =0; blob < this->centers.size() ; blob++){

          string out_string;
          std::stringstream ss;
          ss << blob;
          out_string = ss.str();
           string out_name = outname+out_string+".xvg";
           //cout << out_name << endl;
          ofstream myfile;

          myfile.open (out_name.c_str());

      for(int i =0; i < this->knn_class.size(); i++){
           if(knn_class[i] == blob){
                myfile << dproj[i][0] << " " << dproj[i][1] << " " << knn_class[i] << endl;
           }
       }
        myfile.close();
   }
    std::cout << "Population per Cluster: " ;
      for (std::vector<double>::iterator it=occupancy.begin(); it!=occupancy.end(); ++it){

          std::cout << ' ' << *it;
      }
        std::cout << '\n';


}

void clusterer::print_divided_traj(string input_traj){

    string input_name = input_traj;
        string line;

        ifstream counter(input_name.c_str());

        int number_of_structures =0 ;
        int number_of_atoms =0 ;


        if(!counter) {
          cout << "Cannot open trajectory input file.\n";
          exit (EXIT_FAILURE);
        }
        vector < string > lines_vc;
        while(getline(counter, line)) {
              std::string firstWord = line.substr(0, line.find(" "));
              if(firstWord == "MODEL"){
                number_of_structures++;
              }
              if(number_of_structures == 1 && firstWord == "ATOM"){
                number_of_atoms++;
              }
            }


        for(int i = 0; i < lines_vc.size(); i++){
            cout << lines_vc[i] << endl;
        }
        counter.close();
        cout << "Number of structures " << number_of_structures << endl;
        cout << "Number of atoms per structure " << number_of_atoms << endl;




        ifstream trajectory(input_name.c_str());
        getline(trajectory, line);
        vector < ofstream > file_vectors;
        ofstream file_temp[this->centers.size()];
        string output_clusters_traj;

        for(int i =0; i < this->centers.size(); i++){
            output_clusters_traj="CLUSTER_traj_"+to_string(i)+".pdb";
            file_temp[i].open(output_clusters_traj.c_str());


        }

        for(int mol =0; mol < number_of_structures ; mol++){
            std::string firstWord = line.substr(0, line.find(" "));

            while(firstWord != "MODEL"){
                getline(trajectory, line);
                firstWord = line.substr(0, line.find(" "));

            }



           file_temp[this->knn_class[mol]] << "MODEL        1\n" ;

            for(int atom = 0; atom < number_of_atoms+2 ; atom++){
                getline(trajectory, line);
                std::string firstWord = line.substr(0, line.find(" "));
                file_temp[this->knn_class[mol]] << line << "\n" ;
            }
            getline(trajectory, line);
            firstWord = line.substr(0, line.find(" "));


        }
        counter.close();



}



