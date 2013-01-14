/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */
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

    map<int, Qrels> qrels = readDiversityQrelsFile(qrelsFile);

    // Read the run file
    map<int, vector<Document> > run = readRunFile(runFile);

    // Simulation of preference judgments would go here




    double srecallSum, andcgSum, erriaSum;
    arma::vec prfScoreSum = arma::zeros(8);
    double numOfQ = qrels.size();
    map<int, Qrels>::iterator it;
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        int query = it->first;
        Qrels qrels = it->second;
        int rank = 10;
        arma::mat run_matrix = judge_diversity(run.find(query)->second, qrels, rank);
        cout << query;
/**/
        // Get Subtopic-Recall Score
        arma::vec srecallScore = s_recall(run_matrix, qrels, rank);
        cout << "," << srecallScore(rank - 1);
        srecallSum += srecallScore(rank - 1);

        // Get Alpha-nDCG Score
        arma::vec andcgScore = andcg(run_matrix, qrels, rank);
        cout << "," << andcgScore(rank - 1);
        andcgSum += andcgScore(rank - 1);

        // Get ERR-IA Score
        arma::vec erriaScore = erria(run_matrix, qrels, rank);
        cout << "," << erriaScore(rank -1);
        erriaSum += erriaScore(rank - 1);


        map<string, arma::vec> prfScore = pref_measure(run.find(query)->second, qrels, rank);
        map<string, arma::vec>::iterator prefScore_iter;
        int i =0;
        for(prefScore_iter = prfScore.begin();prefScore_iter != prfScore.end(); prefScore_iter++ ){
                cout << "," << prefScore_iter->second(rank - 1);
                prfScoreSum(i++) += prefScore_iter->second(rank - 1);
        }


        cout << endl;

    }
    cout << "all," << srecallSum/numOfQ << "," << andcgSum/numOfQ << "," << erriaSum/numOfQ;
    map<string, arma::vec>::iterator prefScore_iter;
    for(int i =0; i< 8; i ++ )
            cout << "," << prfScoreSum(i)/numOfQ;
    cout << endl;
    return 0;
}

