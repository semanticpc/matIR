/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */


#include "Utils.hpp"
#include "DivMeasures.hpp"
using namespace std;

int main(int argc, char** argv) {
    string usage = "Usage:\n\t divEval -q <qrels_file> -r <run_file>\n";
    if (argc < 3){
        cout << usage;
        exit(0);
    }
    string qrelsFile = argv[1], runFile = argv[2];
    if (runFile == "" || qrelsFile == ""){
        cout << usage;
        exit(0);
    }

    // Get Qrels File
    map<string, Qrels> qrels = readDiversityQrelsFile(qrelsFile);

    // Read the run file
    map<string, vector<Document> > run = readRunFile(runFile);

    // Simulation of preference judgments would go here




    double srecallSum, andcgSum, erriaSum;
    map<string, Qrels>::iterator it;
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        string query = it->first;
        Qrels qrels = it->second;
        int rank = 10;
        arma::mat run_matrix = judge_diversity(run.find(query)->second, qrels, rank);
        cout << query;

        // Get Subtopic-Recall Score
        double srecallScore = s_recall(run_matrix, qrels, rank);
        cout << " " << srecallScore;
        srecallSum += srecallScore;

        // Get Alpha-nDCG Score
        double andcgScore = andcg(run_matrix, qrels, rank);
        cout << " " << andcgScore;
        andcgSum += andcgScore;

        // Get ERR-IA Score
        double erriaScore = erria(run_matrix, qrels, rank);
        cout << " " << erriaScore;
        erriaSum += erriaScore;

/*
        double prfScore = pref_measure( run[query], qrels, rank);
        cout << " " << prfScore;*/

    }
    return 0;
}

