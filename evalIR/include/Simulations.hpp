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
    vector<Qrels> _qrels_vector;
    double _missing_rate;
    double _error_rate;
    int _error_pairs;
    int _total_pairs;

private:
    int get_preference(arma::rowvec vectorA, arma::rowvec vectorB, arma::rowvec seen);
    pair<arma::vec,arma::vec> get_simulation_scores(Qrels &qrels, vector<int> rankedDocs=vector<int>());
    arma::vec simulate_level(vector<int> rankedDocs);

public:
    double get_UtilityScore(int docIndex);
    pair<int, double> get_BestUtilityDoc();
    double get_UtilityScore(int docIndex, int prevDocIndex);
    pair<int, double> get_BestUtilityDoc(int prevDocIndex);
    void printCounts();
    int getTotalRelDocs();
    int getTotalPairs();

    PrefSimulation(Qrels& qrels, double error_rate=0, double missing_rate=0);
    PrefSimulation(Qrels& qrels, vector<Qrels> qrels_vector, double error_rate=0, double missing_rate=0);

};
#endif	/* SIMULATIONS_HPP */

