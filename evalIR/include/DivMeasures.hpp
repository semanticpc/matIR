/*
 * File:   DivMeasures.hpp
 * Author: semanticpc
 *
 * Created on December 11, 2012, 1:10 PM
 */

#ifndef DIVMEASURES_HPP
#define	DIVMEASURES_HPP

#include <armadillo>
#include <set>
#include <math.h>
#include "DivQrels.hpp"
#include "Simulations.hpp"

using namespace std;

static double s_recall(arma::mat runMatrix, DivQrels* qrels, int rank){
    // Use a greedy approach to find the max subtopics at rank k
    arma::vec tmpVector;

    arma::rowvec seenSubtopics = arma::zeros<arma::rowvec>(qrels->numOfSubtopics);
    set<int> idealRankList;
    set<int>::iterator it;
    for(int i = 0; i < qrels->numOfRelDocuments; i++){
        if(i > rank)
            break;
        int maxNewSubtopics = 0;
        int pickedDoc = -1;

        for(int j = 0; j < qrels->numOfRelDocuments; j++){
            it = idealRankList.find(j);
            if(it != idealRankList.end() )
                continue;

            int numOfNewSubtopics = arma::sum((qrels->subtopicMatrix.row(j) + seenSubtopics) >
                                arma::zeros<arma::vec>(qrels->numOfSubtopics));

            if(numOfNewSubtopics > maxNewSubtopics){
                maxNewSubtopics = numOfNewSubtopics;
                pickedDoc = j;
            }
        }
        seenSubtopics += qrels->subtopicMatrix.row(pickedDoc);
        idealRankList.insert(pickedDoc);
    }
    int i = 0;
    arma::mat idealMatrix = arma::zeros<arma::mat>(rank, qrels->numOfSubtopics);
    for (it=idealRankList.begin(); it!=idealRankList.end(); it++)
        idealMatrix.row(i++) = qrels->subtopicMatrix.row(*it);

    double maxSubtopics = arma::sum(arma::sum(idealMatrix) >
                                        arma::zeros<arma::vec>(qrels->numOfSubtopics));

    // Consider only the top 'rank' number of documents
    runMatrix = runMatrix.rows(0,rank - 1);
    double retSubtopics = arma::sum(arma::sum(runMatrix) >
                                        arma::zeros<arma::vec>(qrels->numOfSubtopics));
    return retSubtopics/maxSubtopics;
}




// Alpha-nDCG Computation
static arma::rowvec compute_andcg_discount(arma::rowvec subtopics, double alpha){
    // Iterate through the elements in the vector and compute (1- alpha)^x
    arma::rowvec discount = arma::zeros<arma::rowvec>(subtopics.n_cols);
    int i =0;
    for(arma::rowvec::iterator it=subtopics.begin(); it!=subtopics.end(); ++it){
        discount(i++) = pow((1 - alpha),(*it));
    }
    return discount;
}


static arma::mat ideal_andcg_matrix(DivQrels* qrels, int rank, double alpha){
    arma::mat idealMatrix = arma::zeros<arma::mat>(rank, qrels->numOfSubtopics);
    arma::rowvec seenSubtopics = arma::zeros<arma::rowvec>(qrels->numOfSubtopics);
    set<int> idealRankList;
    set<int>::iterator it;
    for(int i = 0; i < qrels->numOfRelDocuments; i++){
        if(i > rank)
            break;
        double maxGain = 0;
        int pickedDoc = -1;
        for(int j = 0; j < qrels->numOfRelDocuments; j++){
            it = idealRankList.find(j);
            if(it != idealRankList.end() )
                continue;
            arma::rowvec J_d = qrels->subtopicMatrix.row(j);
            // Compute the discount
            arma::rowvec discount = compute_andcg_discount(seenSubtopics, alpha);

            double gain = arma::sum(J_d % discount) / (log(2 + i)/log(2));
            if(gain > maxGain){
                maxGain = gain;
                pickedDoc = j;
            }
        }
        seenSubtopics += qrels->subtopicMatrix.row(pickedDoc);
        idealRankList.insert(pickedDoc);
        idealMatrix.row(i) = qrels->subtopicMatrix.row(pickedDoc);
    }

    return idealMatrix;
}

static double compute_dcg(arma::mat matrix, int rank, double alpha){
    arma::rowvec seenSubtopics = arma::zeros<arma::rowvec>(matrix.n_cols);
    set<int>::iterator it;
    double dcg = 0;
    for(int i = 0; i < rank; i++){
        arma::rowvec J_d = matrix.row(i);
        arma::rowvec discount = compute_andcg_discount(seenSubtopics, alpha);
        double gain = arma::sum(J_d % discount) / (log(2 + i)/log(2));
        dcg += gain;
        seenSubtopics += matrix.row(i);
    }
    return dcg;
}
static double andcg(arma::mat run_matrix, DivQrels* qrels, int rank, double alpha=0.5){
    // Use a greedy approach to find the ideal rank list

    arma::mat ideal_matrix = ideal_andcg_matrix(qrels, rank, alpha);
    return compute_dcg(run_matrix, rank, alpha) / compute_dcg(ideal_matrix, rank, alpha) ;
}


// ERR Computation
static double computeERR(arma::vec grades, int rank, double alpha){

 double score = 0.0;
 double decay = 1.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {

    double grade = (pow(2.0, grades(i)) - 1) / pow(2.0, 1);
    score += ( grade * decay) / (i + 1);
    decay *= (1 - grade);
 }
 return score;
}


static double computeDCG(arma::vec grades, int rank, double alpha){

 double score = 0.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {
    if(grades(i) == 1)
        r = 1;
    else
        r = 0;
    score += r / log(i + 2);
 }
 return score;
}

static double computePrecision(arma::vec grades, int rank, double alpha){

 double score = 0.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {
    if(grades(i) == 1)
        r = 1;
    else
        r = 0;
    score += r;
 }
 return score/rank;
}


static double erria(arma::mat matrix, DivQrels* qrels, int rank, double alpha=0.5){
    int num_of_subtopics = arma::sum(arma::sum(qrels->subtopicMatrix) > arma::zeros(qrels->numOfSubtopics));
    //int num_of_subtopics = qrels->numOfSubtopics;
    double erria = 0;
    arma::rowvec subtopicGain = arma::ones(matrix.n_cols);
    for (int i = 0; i < rank; i++) {
        double score = 0.0;
        for(int st = 0; st < matrix.n_cols ; st++){
            if(matrix(i, st) == 1){
                score += subtopicGain(st);
                subtopicGain(st) *= (1 - alpha);
            }
        }
        erria += score / (i + 1);
    }
    double norm;
    double normGain = num_of_subtopics;
    for (int i = 0; i < rank; i++){
      norm += normGain /((double)(i + 1));
      normGain *= (1.0 - alpha);
    }
    norm -= num_of_subtopics;
    return erria / norm;
}




// Preference Based Utility Measures
static double f_function(arma::vec utilities, string type="ave"){
    if (type == "ave")
        return arma::mean(utilities);
    else if (type == "min")
        return arma::min(utilities);
}


static double p_function(int rank, string type="none"){
    if (type == "none")
        return 1;
    else if (type == "min")
        return 0;
}


static double get_docUtility_score(int docIndex, int prevDocIndex,
                        DivQrels* qrels, map<string, arma::vec>& cache){

    // Generate a code for the current rankedDocs
    string code; ostringstream convert;
    convert << prevDocIndex;
    code = convert.str();

    if(docIndex == -1)
        return 0;
    else if (docIndex != -1 && prevDocIndex == -1)
        return cache.find("")->second(docIndex);
    else if (docIndex != -1 && prevDocIndex != -1){
        if(cache.find(code) != cache.end())
            return cache.find(code)->second(docIndex);
        else{
            vector<int> rankedDocs;
            rankedDocs.push_back(prevDocIndex);
            arma::vec lvl_score = simulate_level(qrels, rankedDocs);
            cache.insert(make_pair(code, lvl_score));
            return lvl_score(docIndex);
        }
    }
}


static arma::vec get_doc_utilites(DivQrels* qrels, vector<string> rankList,
                                        int rank, map<string, arma::vec>& cache){

    arma::vec doc_utility = arma::zeros(rank);
    vector<int> rankList_docIndex;

    if(rankList.size() < rank){
        for(int i=rankList.size(); i < rank; i++)
            rankList.push_back("dummyDoc");
    }

    // Get the document utilities
    for (int i = 0; i < rank; i++) {
        int docIndex;
        // Obtain the index of the document
        map<string, int>::iterator doc = qrels->docid_index_map.find(rankList[i]);
        if(doc == qrels->docid_index_map.end())
            docIndex = -1;
        else
            docIndex = doc->second;

        // Calculate the document utilities
        if(i == 0)
            doc_utility(i) = get_docUtility_score(docIndex, -1, qrels, cache);
        else{
            arma::vec this_doc_utilities = arma::zeros(i);
            for(int seenDoc=0; seenDoc < i; seenDoc++)
                this_doc_utilities(seenDoc) = get_docUtility_score(docIndex,
                                        rankList_docIndex[seenDoc], qrels, cache);

            doc_utility(i) = f_function(this_doc_utilities);
        }
        rankList_docIndex.push_back(docIndex);
    }
    return doc_utility;
}



static arma::vec get_ideal_utilites(DivQrels* qrels, int rank, map<string, arma::vec>& cache){

    arma::vec ideal_utility = arma::zeros(rank);
    vector<int> rankList_docIndex;
    // Get the document utilities
    for (int i = 0; i < rank; i++) {
        // Get Ideal Document at rank 1
        if(i == 0){
            arma::uvec indices = arma::sort_index(cache.find("")->second, 1);
            ideal_utility(i) = cache.find("")->second(indices(0));
            rankList_docIndex.push_back(indices(0));
        } else{
            double max_utility = 0;
            int max_doc = -1;
            for(int docIndex=0; docIndex < qrels->numOfRelDocuments; docIndex++){
                if(find(rankList_docIndex.begin(), rankList_docIndex.end(), docIndex)!=rankList_docIndex.end())
                    continue;
                arma::vec this_doc_utilities = arma::zeros(i);
                for(int seenDoc=0; seenDoc < i; seenDoc++)
                    this_doc_utilities(seenDoc) = get_docUtility_score(docIndex,
                                        rankList_docIndex[seenDoc],qrels, cache);

                double score =  f_function(this_doc_utilities);
                if(score > max_utility){
                    max_utility = score;
                    max_doc = docIndex;
                }
            }
            ideal_utility(i) = max_utility;
            rankList_docIndex.push_back(max_doc);
        }

    }
    return ideal_utility;
}






static vector<double> compute_pref_score(arma::vec doc_utility, arma::vec ideal_utility,
                                                    int rank, string pType){
    // Compute Preference scores
    vector<double> pref_score;

    double cummulativeScore = 0.0,cummulativeIdealScore = 0.0;
    for (int k = 0; k < rank; k++){
        double totalUtility = 0;
        for (int i = 0; i < k; i++)
                totalUtility += doc_utility(i);
        totalUtility = totalUtility * p_function(k + 1);

        double totalIdealUtility = 0.0;
        for (int i = 0; i < k; i++)
            totalIdealUtility += ideal_utility(i);
        totalIdealUtility = totalIdealUtility * p_function(k + 1);


        cummulativeScore += totalUtility;
        cummulativeIdealScore += totalIdealUtility;
        pref_score.push_back(cummulativeScore / cummulativeIdealScore);
    }
    return pref_score;

}

static double pref_measure(vector<string> rankList, DivQrels* qrels, int rank){


    map<string, arma::vec>  cache;
    cache.insert(make_pair("",simulate_level(qrels)));
    arma::vec ideal_utility = get_ideal_utilites(qrels, rank, cache);
    arma::vec doc_utility = get_doc_utilites(qrels, rankList, rank, cache);

    // Average and None
    string ptype = "none";
    vector<double> ave_none = compute_pref_score(doc_utility, ideal_utility, rank, ptype);
    //vector<double> ave_none = compute_cummulative_pref_scores(doc_utility, ideal_utility, ptype);

    return ave_none[4];
}

#endif	/* DIVMEASURES_HPP */

