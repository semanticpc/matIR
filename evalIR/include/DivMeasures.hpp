/*
 * File:   DivMeasures.hpp
 * Author: semanticpc
 *
 * Created on December 11, 2012, 1:10 PM
 */

#ifndef DIVMEASURES_HPP
#define	DIVMEASURES_HPP


#include <set>
#include <math.h>
#include <armadillo>
#include "Simulations.hpp"

using namespace std;

static arma::vec s_recall(arma::mat &runMatrix, Qrels &qrels, int rank){
    // Use a greedy approach to find the max subtopics at rank k

    // Initialize working variables
    int numOfSubtopics = int(*qrels.subtopics.rbegin());
    int numOfRows = min(rank, qrels.numOfRelDocs);
    arma::vec maxSubtopics = arma::zeros(rank);
    arma::vec tmpVector;
    arma::rowvec seenSubtopics = arma::zeros<arma::rowvec>(numOfSubtopics);
    set<int> idealRankList;
    set<int>::iterator it;

    // Greedy approach searches for the best document at each rank
    for(int i = 0; i < numOfRows; i++){
        int maxNewSubtopics = 0;
        int pickedDoc = -1;
        // Iterate through the set of relevant documents to find
        //   the best document
        for(int j = 0; j < qrels.matrix.n_rows; j++){

            // Ignore the document already picked
            it = idealRankList.find(j);
            if(it != idealRankList.end() )
                continue;

            // Compute the number of new subtopics the document contains
            int numOfNewSubtopics = arma::sum((qrels.matrix.row(j) + seenSubtopics) >
                                arma::zeros<arma::vec>(numOfSubtopics));

            // Keep track of the best document/ max number of subtopics
            if(numOfNewSubtopics > maxNewSubtopics){
                maxNewSubtopics = numOfNewSubtopics;
                pickedDoc = j;
            }
        }

        // Add the best document to the ideal rank list
        //  and keep track of subtopics seen
        seenSubtopics += qrels.matrix.row(pickedDoc);
        idealRankList.insert(pickedDoc);
        maxSubtopics(i) = arma::sum(seenSubtopics >
                                        arma::zeros<arma::vec>(numOfSubtopics));
    }
    if(numOfRows < rank){
        for(int i=numOfRows; i < rank; i++)
            maxSubtopics(i) = maxSubtopics(numOfRows - 1);
    }

    // Consider only the top 'rank' number of documents
    arma::vec s_recall = arma::zeros(rank);
    for(int i=0; i < rank; i++){
        arma::mat matrix = runMatrix.rows(0,i);
        double retSubtopics = arma::sum(arma::sum(matrix) >
                                        arma::zeros<arma::vec>(numOfSubtopics));
        s_recall(i) = (retSubtopics/maxSubtopics(i));
    }
    return s_recall;
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


static arma::mat ideal_andcg_matrix(Qrels &qrels, int rank, double alpha){
    int numOfSubtopics = int(*qrels.subtopics.rbegin());
    arma::mat idealMatrix = arma::zeros<arma::mat>(rank, numOfSubtopics);
    arma::rowvec seenSubtopics = arma::zeros<arma::rowvec>(numOfSubtopics);
    set<int> idealRankList;
    set<int>::iterator it;
    for(int i = 0; i < min(rank, qrels.numOfRelDocs); i++){
        double maxGain = 0;
        int pickedDoc = -1;
        for(int j = 0; j < qrels.matrix.n_rows; j++){
            it = idealRankList.find(j);
            if(it != idealRankList.end() )
                continue;
            arma::rowvec J_d = qrels.matrix.row(j);
            // Compute the discount
            arma::rowvec discount = compute_andcg_discount(seenSubtopics, alpha);

            double gain = arma::sum(J_d % discount) / (log(2 + i)/log(2));
            if(gain > maxGain){
                maxGain = gain;
                pickedDoc = j;
            }
        }
        seenSubtopics += qrels.matrix.row(pickedDoc);
        idealRankList.insert(pickedDoc);
        idealMatrix.row(i) = qrels.matrix.row(pickedDoc);
    }


    return idealMatrix;
}

static arma::vec compute_dcg(arma::mat matrix, int rank, double alpha){
    arma::vec dcg_vector = arma::zeros(rank);
    arma::rowvec seenSubtopics = arma::zeros<arma::rowvec>(matrix.n_cols);
    set<int>::iterator it;
    double dcg = 0;
    for(int i = 0; i < rank; i++){
        arma::rowvec J_d = matrix.row(i);
        arma::rowvec discount = compute_andcg_discount(seenSubtopics, alpha);
        double gain = arma::sum(J_d % discount) / (log(2 + i)/log(2));
        dcg += gain;
        seenSubtopics += matrix.row(i);
        dcg_vector(i) = dcg;
    }
    return dcg_vector;
}
static arma::vec andcg(arma::mat &run_matrix, Qrels &qrels, int rank, double alpha=0.5){
    // Use a greedy approach to find the ideal rank list
    arma::mat ideal_matrix = ideal_andcg_matrix(qrels, rank, alpha);

    return (compute_dcg(run_matrix, rank, alpha) / compute_dcg(ideal_matrix, rank, alpha)) ;
}


// ERR Computation
static arma::vec computeERR(arma::vec grades, int rank, double alpha){
 arma::vec err_vector = arma::zeros(rank);
 double score = 0.0;
 double decay = 1.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {

    double grade = (pow(2.0, grades(i)) - 1) / pow(2.0, 1);
    score += ( grade * decay) / (i + 1);
    decay *= (1 - grade);
    err_vector(score) = score;

 }
 return err_vector;
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


static arma::vec erria(arma::mat &matrix, Qrels &qrels, int rank, double alpha=0.5){
    arma::vec erria_vector = arma::zeros(rank);
    int num_of_subtopics = arma::sum(arma::sum(qrels.matrix) > arma::zeros(qrels.matrix.n_cols));
    double erria = 0;
    arma::rowvec subtopicGain = arma::ones(matrix.n_cols);
    for (int i = 0; i < rank; i++) {
        double score = 0.0;
        for(int st = 0; st < matrix.n_cols ; st++){
            if(matrix(i, st) > 0){
                score += subtopicGain(st);// * qrels.subtopicImportance(st));
                subtopicGain(st) *= (1 - alpha);
            }
        }
        erria += (score / (i + 1));
        erria_vector(i) = erria;
    }
    double norm = 0;
    double normGain = num_of_subtopics;
    arma::vec norm_vector = arma::zeros(rank);
    for (int i = 0; i < rank; i++){
      norm += (normGain /((double)(i + 1)));
      norm_vector(i) = norm;
      normGain *= (1.0 - alpha);
    }
    return (erria_vector / norm_vector);
}




// Preference Based Utility Measures
static double f_function(arma::vec utilities, string type){
    if (type == "ave")
        return arma::mean(utilities);
    else if (type == "min")
        return arma::min(utilities);
}


static double p_function(int rank, string type){
    if (type == "none")
        return 1;
    else if (type == "RBP"){
        double theta = 0.5; // RBP
        return pow((1.0 - theta),rank) * theta; // NRBP
    } else if (type == "RR")
        return 1.0/ (rank * (rank + 1.0)); // Rank
    else if (type == "DCG")
        return (1.0/ log(rank + 1.0)) - (1.0/ log(rank + 2)); // LogProb
}





static arma::vec get_ideal_utilites(Qrels &qrels, PrefSimulation &utility_scores,
                                                        int rank, string fType){

    arma::vec ideal_utility = arma::zeros(rank);
    vector<int> rankList_docIndex;
    // Get the document utilities
    for (int i = 0; i < rank; i++) {
        // Get Ideal Document at rank 1
        if(i == 0){
            pair<int, double> bestDoc = utility_scores.get_BestUtilityDoc();
            ideal_utility(i) = bestDoc.second;
            rankList_docIndex.push_back(bestDoc.first);
        } else{
            double max_utility = 0;
            int max_doc = -1;
            for(int docIndex=0; docIndex < qrels.matrix.n_rows; docIndex++){
                if(find(rankList_docIndex.begin(), rankList_docIndex.end(), docIndex)!=rankList_docIndex.end())
                    continue;

                arma::vec this_doc_utilities = arma::zeros(i);
                for(int seenDoc=0; seenDoc < i; seenDoc++)
                    this_doc_utilities(seenDoc) = utility_scores.get_UtilityScore(docIndex, rankList_docIndex[seenDoc]);

                double score =  f_function(this_doc_utilities, fType);
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




static arma::vec get_doc_utilites(Qrels &qrels, vector<Document> run,
                                PrefSimulation& utility_scores, int rank, string fType){

    arma::vec utility = arma::zeros(rank);
    vector<int> rankList_docIndex;
    map<Document, int>::iterator find_doc;
    int docIndex = -1;
    // Get the document utilities
    for (int i = 0; i < min(rank, int(run.size())); i++) {
        find_doc = qrels.relDocs.find(run.at(i));
        if(find_doc != qrels.relDocs.end())
            docIndex = find_doc->second;
        else
            docIndex = -1;
        // Get Ideal Document at rank 1
        if(i == 0){
            utility(i) = utility_scores.get_UtilityScore(docIndex);
            rankList_docIndex.push_back(docIndex);
        } else{
            arma::vec this_doc_utilities = arma::zeros(i);
            for(int seenDoc=0; seenDoc < i; seenDoc++)
                this_doc_utilities(seenDoc) = utility_scores.get_UtilityScore(docIndex, rankList_docIndex[seenDoc]);

            utility(i) = f_function(this_doc_utilities, fType);
            rankList_docIndex.push_back(docIndex);
        }

    }
    return utility;
}

 static arma::vec compute_pref_score(arma::vec doc_utility, arma::vec ideal_utility,
                                                    int rank, string pType){
    // Compute Preference scores
    arma::vec pref_score = arma::zeros(rank);

    double cummulativeScore = 0.0,cummulativeIdealScore = 0.0;
    for (int k = 0; k < rank; k++){
        double totalUtility = 0;
        for (int i = 0; i < k; i++)
                totalUtility += doc_utility(i);
        totalUtility = totalUtility * p_function(k + 1, pType);

        double totalIdealUtility = 0.0;
        for (int i = 0; i < k; i++)
            totalIdealUtility += ideal_utility(i);
        totalIdealUtility = totalIdealUtility * p_function(k + 1,pType);


        cummulativeScore += totalUtility;
        cummulativeIdealScore += totalIdealUtility;
        pref_score(k) = (cummulativeScore / cummulativeIdealScore);
    }
    return pref_score;

}
// Some Testing Functions
static void simulationTesting(PrefSimulation& utility_scores, Qrels& qrels){
    cout << endl;
    cout << endl;
    for(int i=0; i< qrels.matrix.n_rows; i++){
        cout << i << " " << utility_scores.get_UtilityScore(i) << endl;
    }
    for(int k=0; k< qrels.matrix.n_rows; k++){
    cout << endl;
    cout << endl;
    cout << k << endl;
    for(int i=0; i< qrels.matrix.n_rows; i++){
        cout << i << " " << utility_scores.get_UtilityScore(i, k) << endl;
    }}

    cout << endl;
    cout << endl;
    map<Document, int>::iterator it;
    for(it=qrels.relDocs.begin(); it != qrels.relDocs.end(); it++)
        cout << it->first.docid << " " << it->second << endl;



}
static map<string, arma::vec> pref_measure(vector<Document>& run, int rank, vector<Qrels> qrel_vector, Qrels& qrels){
    PrefSimulation utility_scores(qrels, qrel_vector);
    //Qrels qrels = qrel_vector.front();
    map<string, arma::vec> res;
    //simulationTesting(utility_scores, qrels);


    arma::vec ideal_utility_ave = get_ideal_utilites(qrels, utility_scores, rank, "ave");
    arma::vec doc_utility_ave = get_doc_utilites(qrels, run, utility_scores, rank, "ave");

    arma::vec ave_none = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "none");
    arma::vec ave_RBP = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "RBP");
    arma::vec ave_RR = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "RR");
    arma::vec ave_DCG = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "DCG");

    res.insert(make_pair("ave_none", ave_none));
    res.insert(make_pair("ave_RBP", ave_RBP));
    res.insert(make_pair("ave_RR", ave_RR));
    res.insert(make_pair("ave_DCG", ave_DCG));


    arma::vec ideal_utility_min = get_ideal_utilites(qrels, utility_scores, rank, "min");
    arma::vec doc_utility_min = get_doc_utilites(qrels, run, utility_scores, rank, "min");

    arma::vec min_none = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "none");
    arma::vec min_RBP = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "RBP");
    arma::vec min_RR = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "RR");
    arma::vec min_DCG = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "DCG");

    res.insert(make_pair("min_none", min_none));
    res.insert(make_pair("min_RBP", min_RBP));
    res.insert(make_pair("min_RR", min_RR));
    res.insert(make_pair("min_DCG", min_DCG));

    return res;
}
static map<string, arma::vec> pref_measure(vector<Document>& run, Qrels& qrels, int rank){
    PrefSimulation utility_scores(qrels);
    map<string, arma::vec> res;
    //simulationTesting(utility_scores, qrels);


    arma::vec ideal_utility_ave = get_ideal_utilites(qrels, utility_scores, rank, "ave");

    arma::vec doc_utility_ave = get_doc_utilites(qrels, run, utility_scores, rank, "ave");

    arma::vec ave_none = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "none");
    arma::vec ave_RBP = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "RBP");
    arma::vec ave_RR = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "RR");
    arma::vec ave_DCG = compute_pref_score(doc_utility_ave, ideal_utility_ave, rank, "DCG");

    res.insert(make_pair("ave_none", ave_none));
    res.insert(make_pair("ave_RBP", ave_RBP));
    res.insert(make_pair("ave_RR", ave_RR));
    res.insert(make_pair("ave_DCG", ave_DCG));


    arma::vec ideal_utility_min = get_ideal_utilites(qrels, utility_scores, rank, "min");
    arma::vec doc_utility_min = get_doc_utilites(qrels, run, utility_scores, rank, "min");

    arma::vec min_none = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "none");
    arma::vec min_RBP = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "RBP");
    arma::vec min_RR = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "RR");
    arma::vec min_DCG = compute_pref_score(doc_utility_min, ideal_utility_min, rank, "DCG");

    res.insert(make_pair("min_none", min_none));
    res.insert(make_pair("min_RBP", min_RBP));
    res.insert(make_pair("min_RR", min_RR));
    res.insert(make_pair("min_DCG", min_DCG));

    return res;
}




#endif	/* DIVMEASURES_HPP */