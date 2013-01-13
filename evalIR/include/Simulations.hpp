/*
 * File:   PrefMeasure.hpp
 * Author: semanticpc
 *
 * Created on December 13, 2012, 9:12 AM
 */

#ifndef SIMULATIONS_HPP
#define	SIMULATIONS_HPP

#include "Utils.hpp"


using namespace std;

static int get_preference(arma::rowvec vectorA, arma::rowvec vectorB, arma::rowvec seen){
    double alpha = 0.5;
    double scoreA = 0, scoreB = 0;
    for(int i=0; i < vectorA.n_cols; i++){
        if(vectorA(i) == 1)
            scoreA += pow((1- alpha), seen(i));
    }
    for(int i=0; i < vectorB.n_cols; i++){
        if(vectorB(i) == 1)
            scoreB +=  pow((1 - alpha), seen(i));
    }
    if (scoreA > scoreB)
        return 1;
    else if (scoreA < scoreB)
        return 0;
    else
        return 1;


}

static arma::vec simulate_level(Qrels* qrels, vector<int> rankedDocs=vector<int>()){


    // Initialize Seen Subtopic Counts
    arma::rowvec seen = arma::zeros<arma::rowvec>(qrels->numOfSubtopics);
    for(vector<int>::iterator i= rankedDocs.begin(); i != rankedDocs.end(); i++){
        //cout << endl;
        //cout << endl;
        //cout << "topDoc:" << (*i) << endl;
        if((*i) != -1)
            seen += qrels->subtopicMatrix.row((*i));
    }

    arma::vec lvl_score = arma::zeros(qrels->numOfRelDocuments);
    arma::vec appearance_count = arma::zeros(qrels->numOfRelDocuments);

    for(int i=0; i < qrels->numOfRelDocuments; i++){
        if(find(rankedDocs.begin(), rankedDocs.end(), i) != rankedDocs.end())
            continue;
        double score = 0.0;
        int appearance = 0;
        for(int j= i; j < qrels->numOfRelDocuments; j++){
            if( i == j || (find(rankedDocs.begin(), rankedDocs.end(), j) != rankedDocs.end()) )
                continue;
            if(get_preference(qrels->subtopicMatrix.row(i),
                            qrels->subtopicMatrix.row(j), seen) == 1)
                lvl_score(i) += 1;
            else
                lvl_score(j) += 1;

            appearance_count(i) += 1;
            appearance_count(j) += 1;
        }
    }
    //lvl_score.t().print("level count");
    //appearance_count.t().print("appear");
    lvl_score /= appearance_count;
    return lvl_score;
}


#endif	/* SIMULATIONS_HPP */

