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
    cout << "topic,runid,srecall@5,srecall@10,srecall@20";
    cout << ",andcg@5,andcg@10,andcg@20";
    cout << ",erria@5,erria@10,erria@20";
    cout << ",prf_ave_none@5,prf_ave_none@10,prf_ave_none@20";
    cout << ",prf_ave_RBP@5,prf_ave_RBP@10,prf_ave_RBP@20";
    cout << ",prf_ave_RR@5,prf_ave_RR@10,prf_ave_RR@20";
    cout << ",prf_ave_DCG@5,prf_ave_DCG@10,prf_ave_DCG@20";

    cout << ",prf_min_none@5,prf_min_none@10,prf_min_none@20";
    cout << ",prf_min_RBP@5,prf_min_RBP@10,prf_min_RBP@20";
    cout << ",prf_min_RR@5,prf_min_RR@10,prf_min_RR@20";
    cout << ",prf_min_DCG@5,prf_min_DCG@10,prf_min_DCG@20";

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

static void printResultsFolder(string runFolderPath, vector<string> runFiles, map<int, Qrels> qrels, int e=0, int m= 0){


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
    if(e > 0 || m > 0)
        cout << ",TotalPairs, TotalRelDocs";

    cout << endl;
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        int query = it->first;
        Qrels qrels = it->second;

        PrefSimulation utility_scores(qrels, vector<Qrels>(), e, m);
        for(int run_index=0;run_index<runFiles.size();run_index++){
            if(e > 0 || m > 0){
                for(int i=0; i< qrels.matrix.n_rows; i++)
                utility_scores.get_UtilityScore(0, i);
            }
            arma::mat run_matrix = judge_diversity(runs.at(run_index).find(query)->second, qrels, rank);
            cout << query;
            cout << "," << runFiles.at(run_index);

            // Get Subtopic-Recall Score
            arma::vec srecallScore = s_recall(run_matrix, qrels, rank);
            cout << "," << srecallScore(4) << "," << srecallScore(9) << "," << srecallScore(19);

            // Get Alpha-nDCG Score
            arma::vec andcgScore = andcg(run_matrix, qrels, rank);
            cout << "," << andcgScore(4) << "," << andcgScore(9) << "," << andcgScore(19);

            // Get ERR-IA Score
            arma::vec erriaScore = erria(run_matrix, qrels, rank);
            cout << "," << erriaScore(4) << "," << erriaScore(9) << "," << erriaScore(19);

            // Print all preference measure scores
            map<string, arma::vec> prfScore = pref_measure(runs.at(run_index).find(query)->second, qrels, rank, utility_scores);
            for(prefScore_iter = prfScore.begin();prefScore_iter != prfScore.end(); prefScore_iter++ ){
                    cout << "," << prefScore_iter->second(4) << "," << prefScore_iter->second(9)
                         << "," << prefScore_iter->second(19);
            }
            if(e > 0 || m > 0)
                cout << "," << utility_scores.getTotalPairs() << "," << utility_scores.getTotalRelDocs() << endl;
            else
                cout << endl;

        }
    }
}

static void printResults(map<int, vector<Document> > run, map<int, Qrels> qrels, double e=0, double m=0){


    double numOfQ = qrels.size();
    map<int, Qrels>::iterator it;
    map<string, arma::vec>::iterator prefScore_iter;
    int rank = 20;

    // Initialize variable to keep track of the sum of scores
    arma::vec allScores = arma::zeros(33);

    printHeader();
    if(e > 0 || m > 0)
        cout << ",TotalPairs, TotalRelDocs";

    cout << endl;
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        int i =0;
        int query = it->first;
        Qrels qrels = it->second;
        arma::mat run_matrix = judge_diversity(run.find(query)->second, qrels, rank);
        PrefSimulation utility_scores(qrels, vector<Qrels>(), e, m);
        if(e > 0 || m > 0){
            for(int i=0; i< qrels.matrix.n_rows; i++)
                utility_scores.get_UtilityScore(0, i);
        }
        cout << query << ",sysrun";
        //run_matrix.print("runs");
        // Get Subtopic-Recall Score
        arma::vec srecallScore = s_recall(run_matrix, qrels, rank);
        cout << "," << srecallScore(4) << "," << srecallScore(9) << "," << srecallScore(19);

        allScores(i++) += srecallScore(4);
        allScores(i++) += srecallScore(9);
        allScores(i++) += srecallScore(19);

        // Get Alpha-nDCG Score
        arma::vec andcgScore = andcg(run_matrix, qrels, rank);
        cout << "," << andcgScore(4) << "," << andcgScore(9) << "," << andcgScore(19);
        allScores(i++) += andcgScore(4);
        allScores(i++) += andcgScore(9);
        allScores(i++) += andcgScore(19);

        // Get ERR-IA Score
        arma::vec erriaScore = erria(run_matrix, qrels, rank);
        cout << "," << erriaScore(4) << "," << erriaScore(9) << "," << erriaScore(19);
        allScores(i++) += erriaScore(4);
        allScores(i++) += erriaScore(9);
        allScores(i++) += erriaScore(19);


        // Print all preference measure scores
        map<string, arma::vec> prfScore = pref_measure(run.find(query)->second, qrels, rank, utility_scores);
        for(prefScore_iter = prfScore.begin();prefScore_iter != prfScore.end(); prefScore_iter++ ){
                cout << "," << prefScore_iter->second(4) << "," << prefScore_iter->second(9)
                     << "," << prefScore_iter->second(19);
                allScores(i++) += prefScore_iter->second(4);
                allScores(i++) += prefScore_iter->second(9);
                allScores(i++) += prefScore_iter->second(19);
        }
        if(e > 0 || m > 0)
            cout << "," << utility_scores.getTotalPairs() << "," << utility_scores.getTotalRelDocs() << endl;
        else
            cout << endl;

    }

    // Print the mean scores for all measures
    cout << "all,sysrun";
    for(int i =0; i< 33; i ++ ){
        cout << "," << allScores(i)/numOfQ;
    }
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
    }
   if ((runFolder == "" && runFile == "") || qrelsFile == ""){
        cout << usage;
        exit(0);
   }
    // Get Qrels File
    map<int, Qrels> qrels = readDiversityQrelsFile(qrelsFile);

    // Judging User Profiles vs Traditional TREC Assessments
    if(runFolder == ""){
        // Read the run file
        map<int, vector<Document> > run = readRunFile(runFile);
        printResults(run, qrels, error, missing);

    }else{
        vector<string> runFiles;
        getdir(runFolder, runFiles);
        printResultsFolder(runFolder, runFiles, qrels,error, missing);
    }



    return 0;
}

