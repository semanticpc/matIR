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
struct Triplet{
    string code;
    int topDoc;
    int leftDoc;
    int rightDoc;
};


struct triplet_comparison{
    inline bool operator () ( const Triplet& lhs, const Triplet& rhs ) const{
        return lhs.code < rhs.code;

    }
};
class AMTSimulation{


private:
    map<string, arma::vec>  _utilScores;
    map<string, arma::vec>  _appearanceCounts;
    Qrels& _qrels;
    vector<Qrels> _qrels_vector;
    int _error_pairs;
    int _total_pairs;





private:
    int get_preference(arma::rowvec vectorA, arma::rowvec vectorB, arma::rowvec seen);
    void simulateScores();
    map<Triplet, int, triplet_comparison> sampleTriplets(vector<int> );

    void updateScores(int topDoc, int leftDoc, int rightDoc, bool preferLeft);
    bool checkPairs(map<Triplet, int, triplet_comparison> &seen_pairs,
                        int left, int right);

public:
    double get_UtilityScore(int docIndex);
    pair<int, double> get_BestUtilityDoc();
    double get_UtilityScore(int docIndex, int prevDocIndex);
    map<Document, int, docid_comparison> doc_index_map;
    vector<Document> allDocs;
    void printCounts();
    int getTotalRelDocs();
    int getTotalPairs();

    AMTSimulation(Qrels& qrels);
    AMTSimulation(Qrels& qrels, vector<Qrels> qrels_vector);

};
#endif	/* SIMULATIONS_HPP */

