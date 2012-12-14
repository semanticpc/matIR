/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */

#include <iostream>
#include "DivQrels.hpp"
#include "DivMeasures.hpp"
using namespace std;

map<string, vector<string> > readRunFile(string runFileName){
    map<string, vector<string> > run;
    std::ifstream runFile(runFileName.c_str(), ios_base::in);
    std::string docid;
    std::string query;
    std::string q0;
    std::string rank;
    std::string score;
    std::string runid;

    string curQuery = "";
    bool doneParsing = false;
    runFile >> query >> q0 >> docid >> rank >> score >> runid;
    curQuery = query;
    vector<string> documents;
    documents.push_back(docid);
    while (true) {
        while(curQuery == query){
            runFile >> query >> q0 >> docid >> rank >> score >> runid;
            documents.push_back(docid);
            if( runFile.eof() ) {
                doneParsing = true;
                break;
            }

        }
        run.insert(make_pair(curQuery,documents));

        if(doneParsing) break;
        curQuery = query;
        vector<string> documents;
        documents.push_back(docid);

        }
    return run;
}

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

    // Get Qrels
    map<string, DivQrels*> qrels = DivQrels::readDiversityQrelsFile(qrelsFile);

    // Simulation of preference judgments would go here


    // Read the run file
    map<string, vector<string> > run = readRunFile(runFile);

    double srecallSum, andcgSum, erriaSum;
    map<string, DivQrels*>::iterator it;
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        string query = it->first;
        cout << query;
        DivQrels* qrels = it->second;
        arma::vec runMatrix =  qrels->getSubtopicMatrix(run[query]);
        int rank = 10;
        // Get Subtopic-Recall Score
        double srecallScore = s_recall(runMatrix, qrels, rank);
        cout << " " << srecallScore;
        srecallSum += srecallScore;

        // Get Alpha-nDCG Score
        double andcgScore = andcg(runMatrix, qrels, rank);
        cout << " " << andcgScore;
        andcgSum += andcgScore;

        // Get ERR-IA Score
        double erriaScore = erria(runMatrix, qrels, rank);
        cout << " " << erriaScore;
        erriaSum += erriaScore;


        double prfScore = pref_measure( run[query], qrels, rank);
        cout << " " << prfScore;

    }
    return 0;
}

