/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */
#include <iostream>
#include <fstream>
#include <armadillo>
#include "DivQrels.hpp"

using namespace std;

DivQrels::DivQrels(){
}

DivQrels::DivQrels(string q){
    query = q;
    numOfSubtopics = 0;
    numOfRelDocuments = 0;
}

void DivQrels::generateMatrix(){
    string docid;
    int subtopic;
    int indexCounter = -1;
    int currentIndex;
    subtopicMatrix.set_size(numOfRelDocuments, numOfSubtopics);
    for(std::vector<int>::size_type i = 0; i != _tmpVector.size(); i++) {
        docid = _tmpVector[i].first;
        subtopic = _tmpVector[i].second;
        currentIndex = docid_index_map.find(docid)->second;
        subtopicMatrix(currentIndex, subtopic-1) = 1;
    }

    _tmpVector.clear();

}

void DivQrels::addDocument(string docid, int subtopic, int rel){
    if (rel >= 1){
        _tmpVector.push_back(make_pair(docid, subtopic));
        if ( docid_index_map.find(docid) == docid_index_map.end()){
            docid_index_map.insert(make_pair(docid, numOfRelDocuments++));
        }

    }
    if(subtopic > numOfSubtopics)
        numOfSubtopics = subtopic;

}


arma::mat DivQrels::getSubtopicMatrix(vector<string> run){
        // Set a maximum to 10,000 so its easy to slice the matrix later
        arma::mat runMatrix = arma::zeros<arma::mat>(10, numOfSubtopics);
        int rank = 0;

        for(vector<string>::iterator iter = run.begin(); iter != run.end(); ++iter){
            map<string, int>::iterator docFound = docid_index_map.find((*iter));
            if(docFound != docid_index_map.end())
                runMatrix.row(rank) = subtopicMatrix.row(docFound->second);
            else
                runMatrix.row(rank) = arma::zeros<arma::vec>(numOfSubtopics);
            ++rank;
        }
        return runMatrix;
    }







