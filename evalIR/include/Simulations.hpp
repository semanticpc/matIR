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

class PrefSimulation{
private:
    map<string, arma::vec>  _utilScores;
    Qrels& _qrels;

private:
    int get_preference(arma::rowvec vectorA, arma::rowvec vectorB, arma::rowvec seen);
    arma::vec simulate_level(vector<int> rankedDocs);

public:
    double get_UtilityScore(int docIndex);
    pair<int, double> get_BestUtilityDoc();
    double get_UtilityScore(int docIndex, int prevDocIndex);
    pair<int, double> get_BestUtilityDoc(int prevDocIndex);
    //double get_UtilityScore(int docIndex, vector<int> rankedDocs);

    PrefSimulation(Qrels& qrels);

};
#endif	/* SIMULATIONS_HPP */

