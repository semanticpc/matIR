//============================================================================
// Name        : matIR.cpp
// Author      : Praveen Chandar
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include "matIR/ResultStats.hpp"
#include "matIR/Postretrieval.hpp"
#include "indri/greedy_vector"
#include <queue>
#include <math.h>

void matIR::postretrieval::query_clarity( matIR::LanguageModel& lm_model,
        indri::api::QueryEnvironment& env,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    double clarityScore = 0;
    lm_model.generateRelevanceModel();
    indri::utility::greedy_vector< std::pair< string, double > > terms = lm_model.getScoredTerms();
    int count = 0;
    double sum=0, ln_Pr=0;
    for( size_t j=0; j < terms.size(); j++ ) {
        std::string t = terms[j].first;
        count++;
        // query-clarity = SUM_w{P(w|Q)*log(P(w|Q)/P(w))}
        // P(w)=cf(w)/|C|
        // the relevance model uses stemmed terms, so use stemCount
        double pw = ((double)env.stemCount(t)/(double)env.termCount());
        // P(w|Q) is a prob computed by any model, e.g. relevance models
        double pwq = terms[j].second;
        sum += pwq;
        clarityScore += (pwq)*log(pwq/pw);
    }
    clarityScore = (clarityScore/(sum ? sum : 1.0)/log(2.0));

    feature_scores.push_back(std::make_pair("queryClarity", clarityScore));
    return;
}

void matIR::postretrieval::NQC( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){
    double mu = arma::mean(stats.docScores);
    int k = stats.docScores.size();
    double score_D = arma::sum(stats.getQueryStats().getQuerycTFs()
                                                / stats.collectionLength);

    double nqc = arma::sum((stats.docScores - mu) / score_D);
    nqc /= k;

    feature_scores.push_back(std::make_pair("NQC", nqc));
}

void matIR::postretrieval::retScore_related( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    feature_scores.push_back(std::make_pair("maxRetrievalScore", arma::max(stats.docScores)));
    feature_scores.push_back(std::make_pair("meanRetrievalScore", arma::mean(stats.docScores)));
    feature_scores.push_back(std::make_pair("varRetrievalScore", arma::var(stats.docScores)));
}



void matIR::postretrieval::weighted_info_gain( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    arma::vec query_ctf = stats.getQueryStats().getQuerycTFs();
    double lambda = (1 / sqrt(stats.getQueryStats().getQueryTermIds().size()));
    int k = stats.docLength.size();

    matIR::scoring::TermScoreFunction* termScorer =
                        matIR::scoring::TermScoreFunctionFactory::get( "dirichlet",
                        query_ctf, stats.collectionLength, 0 , 0 );

    // Obtain the tf matrix of retrieved document for only the query terms
    arma::mat prob_t_d = stats.tfMat.cols(stats.getQueryStats().getQueryTermIds());

    // Probability(term | document)
    arma::vec termFrequency_Md;
    int length_Md;
    for (int i=0; i < prob_t_d.n_rows ; i++) {
        // Need to transpose to be consistent with other vectors
        termFrequency_Md = arma::trans(prob_t_d.row(i));
        length_Md = (int) stats.docLength[i];
        prob_t_d.row(i) = arma::exp(termScorer->scoreOccurrence(termFrequency_Md, length_Md));
    }


    // Probability(term | collection)
    arma::vec prob_t_collection = stats.getQueryStats().getQuerycTFs() / stats.collectionLength;
    prob_t_d.each_row() /= prob_t_collection;
    double score = arma::sum(arma::sum((arma::log(prob_t_d) * lambda)));
    score /= k;
    feature_scores.push_back(std::make_pair("weightedInfoGain", score));

}

void matIR::postretrieval::query_feedback( matIR::LanguageModel& lm_model,
        indri::api::QueryEnvironment& env,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    double clarityScore = 0;
    lm_model.setSmoothing("kl");
    lm_model.generateLanguageModel();
    indri::utility::greedy_vector< std::pair< string, double > > terms = lm_model.getScoredTerms();
    int count = 0;
    double sum=0, ln_Pr=0;
    for( size_t j=0; j < terms.size(); j++ ) {
        std::string t = terms[j].first;
        cout << t << endl;
    }


    //feature_scores.push_back(std::make_pair("queryClarity", clarityScore));
    return;
}




