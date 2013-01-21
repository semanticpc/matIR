/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */
#include "DivMeasures.hpp"
#include <dirent.h>
#include <errno.h>
#include <string.h>
using namespace std;

static void printHeader(){
    cout << "topic,profileID,runid";
    cout << ",srecall@5,srecall@10,srecall@20";
    cout << ",P-IA@5,P-IA@10,P-IA@20";
    cout << ",andcg-IA@5,andcg-IA@10,andcg-IA@20";
    cout << ",erria@5,erria@10,erria@20";
    cout << ",prf_ave_none@5,prf_ave_none@10,prf_ave_none@20";
    cout << ",prf_ave_RBP@5,prf_ave_RBP@10,prf_ave_RBP@20";
    cout << ",prf_ave_RR@5,prf_ave_RR@10,prf_ave_RR@20";
    cout << ",prf_ave_DCG@5,prf_ave_DCG@10,prf_ave_DCG@20";

    cout << ",prf_min_none@5,prf_min_none@10,prf_min_none@20";
    cout << ",prf_min_RBP@5,prf_min_RBP@10,prf_min_RBP@20";
    cout << ",prf_min_RR@5,prf_min_RR@10,prf_min_RR@20";
    cout << ",prf_min_DCG@5,prf_min_DCG@10,prf_min_DCG@20";
    cout << endl;
}

int getdir (string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if(dirp->d_name[0] != '.')
            files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

static void printResultsFolder(string runFolderPath, vector<string> runFiles, map<int, Qrels> qrels,
        map<int, map<int, Profiles* > > &profiles){


    // Iterate through the files and store the results in a vector
    vector<map<int, vector<Document> > > runs;
    for(int i=0;i<runFiles.size();i++){
        runs.push_back( readRunFile(runFolderPath + "/" +runFiles.at(i)));
    }
    double numOfQ = qrels.size();
    map<int, Qrels>::iterator it;
    map<string, arma::vec>::iterator prefScore_iter;
    int rank = 20;

    // Initialize variable to keep track of the sum of scores
    arma::vec allScores = arma::zeros(33);

    printHeader();

    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        int i =0;
        int query = it->first;
        Qrels qrels = it->second;
        PrefSimulation utility_scores(qrels, vector<Qrels>());



        // Iterate user profiles
        map<int, map<int, Profiles* > >::iterator iter;
        double srecallSum_query = 0, andcgSum_query = 0, erriaSum_query = 0;
        int numOfProfiles = 0;
        vector<Qrels> new_qrels_vector;
        for(iter=profiles.begin();iter!=profiles.end();iter++){

            Qrels new_qrels;
            map<int, Profiles* >::iterator find_query = iter->second.find(query);
            if(find_query != iter->second.end())
                new_qrels = update_SubtopicImportance(qrels, find_query->second->subtopic_importance);
            else
                continue;

            if(arma::sum(new_qrels.subtopicImportance) <= 0)
                continue;
            new_qrels_vector.push_back(new_qrels);
            numOfProfiles++;
        }

        for(int run_index=0;run_index<runFiles.size();run_index++){

            map<string, arma::vec> prfScore = pref_measure(runs.at(i).find(query)->second, qrels, rank,new_qrels_vector);
            map<string, arma::vec>::iterator prefScore_iter;
            int i =0;
            for(prefScore_iter = prfScore.begin();prefScore_iter != prfScore.end(); prefScore_iter++ ){
                    cout << "," << prefScore_iter->second(rank - 1);
            }
            cout << endl;
        }
    }
}

static void printResults(map<int, vector<Document> > &run, map<int, Qrels> &qrels,
                                        map<int, map<int, Profiles* > > &profiles){
    arma::vec prfScoreSum = arma::zeros(8);
    double numOfQ = qrels.size();
    map<int, Qrels>::iterator it;

    cout << "topic,profileID,runid";
    cout << ",P-IA@5,P-IA@10,P-IA@20";
    cout << ",andcg-IA@5,andcg-IA@10,andcg-IA@20";
    cout << ",erria@5,erria@10,erria@20";
    cout << ",prf_ave_none@5,prf_ave_none@10,prf_ave_none@20";
    cout << ",prf_ave_RBP@5,prf_ave_RBP@10,prf_ave_RBP@20";
    cout << ",prf_ave_RR@5,prf_ave_RR@10,prf_ave_RR@20";
    cout << ",prf_ave_DCG@5,prf_ave_DCG@10,prf_ave_DCG@20";

    cout << ",prf_min_none@5,prf_min_none@10,prf_min_none@20";
    cout << ",prf_min_RBP@5,prf_min_RBP@10,prf_min_RBP@20";
    cout << ",prf_min_RR@5,prf_min_RR@10,prf_min_RR@20";
    cout << ",prf_min_DCG@5,prf_min_DCG@10,prf_min_DCG@20";
    cout << endl;
    // Iterate through each query
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        int query = it->first;
        Qrels qrels = it->second;
        int rank = 20;
        // Iterate user profiles
        map<int, map<int, Profiles* > >::iterator iter;
        double srecallSum_query = 0, andcgSum_query = 0, erriaSum_query = 0;
        int numOfProfiles = 0;
        for(iter=profiles.begin();iter!=profiles.end();iter++){

            Qrels new_qrels;
            map<int, Profiles* >::iterator find_query = iter->second.find(query);
            if(find_query != iter->second.end())
                new_qrels = update_SubtopicImportance(qrels, find_query->second->subtopic_importance);
            else
                continue;

            if(arma::sum(new_qrels.subtopicImportance) <= 0)
                continue;
            arma::mat run_matrix = judge_diversity(run.find(query)->second, new_qrels, rank);
            new_qrels.matrix.print("qrels");
            new_qrels.subtopicImportance.print("stimp");
            run_matrix.print("run");
            cout << query << "," << numOfProfiles;
            cout << ",RUN";


            arma::vec preciaScore = precia(run_matrix, new_qrels, rank);
            cout << "," << preciaScore(2) << "," << preciaScore(9) << "," << preciaScore(19);

            arma::vec ndcgiaScore = ndcgia(run_matrix, new_qrels, rank);
            cout << "," << ndcgiaScore(2) << "," << ndcgiaScore(9) << "," << ndcgiaScore(19);

            arma::vec erriaScore = erria(run_matrix, new_qrels, rank);
            cout << "," << erriaScore(2) << "," << erriaScore(9) << "," << erriaScore(19);


            /*map<string, arma::vec> prfScore = pref_measure(run.find(query)->second, qrels, rank,qrels_vector);
            map<string, arma::vec>::iterator prefScore_iter;
            int i =0;
            for(prefScore_iter = prfScore.begin();prefScore_iter != prfScore.end(); prefScore_iter++ ){
                cout << "," << prefScore_iter->second(rank - 1);
            prfScoreSum(i++) += prefScore_iter->second(rank - 1);*/
            numOfProfiles++;
            cout << endl;
        }
    }
    map<string, arma::vec>::iterator prefScore_iter;
    for(int i =0; i< 8; i ++ )
            cout << "," << prfScoreSum(i)/numOfQ;
    cout << endl;
}

int main(int argc, char** argv) {
    string usage = "Usage:\n\t divEval <qrels_file> <run_file>\n";
    if (argc < 3){
        cout << usage;
        exit(0);
    }
    string qrelsFile, runFile = "", runFolder;
    string profilesFile = "";
    double error, missing;
    //int nextIndex = 3;

    for(int nextIndex= 1; nextIndex < argc; nextIndex++){
        if(strcmp(argv[nextIndex], "-q") == 0){
            qrelsFile = argv[nextIndex + 1];
            nextIndex += 1;
        }
        else if(strcmp(argv[nextIndex], "-rf") == 0){
            runFolder = argv[nextIndex + 1];
            nextIndex += 1;
        }
        else if(strcmp(argv[nextIndex], "-r") == 0){
            runFile = argv[nextIndex + 1];
            nextIndex += 1;
        }
        else if(strcmp(argv[nextIndex], "-e") == 0){
            error = atof(argv[nextIndex + 1]);
            nextIndex += 1;
        }
        else if(strcmp(argv[nextIndex], "-m") == 0 ){
            missing = atof(argv[nextIndex + 1]);
            nextIndex += 1;
        }
        else if(strcmp(argv[nextIndex], "-profile") == 0 ){
            profilesFile = argv[nextIndex + 1];
            nextIndex += 1;
        }
    }
   if ((runFolder == "" && runFile == "") || qrelsFile == "" || profilesFile == "" ){
        cout << usage;
        exit(0);
   }
    // Get Qrels File
    map<int, Qrels> qrels = readDiversityQrelsFile(qrelsFile);




    // Judging User Profiles vs Traditional TREC Assessments
    // Judging User Profiles vs Traditional TREC Assessments
    if(runFolder == ""){
        map<int, map<int, Profiles* > > profiles = read_userProfiles(profilesFile);
        map<int, vector<Document> > run = readRunFile(runFile);
        printResults(run, qrels, profiles);
    }else{
        vector<string> runFiles;
        getdir(runFolder, runFiles);
        map<int, map<int, Profiles* > > profiles = read_userProfiles(profilesFile);
        printResultsFolder(runFolder, runFiles, qrels, profiles);
    }



    return 0;
}



