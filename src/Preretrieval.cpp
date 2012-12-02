//============================================================================
// Name        : matIR.cpp
// Author      : Praveen Chandar
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================
#include "indri/QueryEnvironment.hpp"
#include "matIR/ResultStats.hpp"
#include "matIR/Preretrieval.hpp"
#include "indri/greedy_vector"

#include <queue>

void matIR::preretrieval::idfRelated( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    // Compute IDF from qstats
    arma::vec df = stats.getQueryStats().getQueryDFs();
    arma::vec idf = arma::log((stats.documentCount + 1) / (df + 0.5));
    feature_scores.push_back(std::make_pair("idfSum", arma::sum(idf)));
    feature_scores.push_back(std::make_pair("idfStdDev", arma::stddev(idf)));
    feature_scores.push_back(std::make_pair("idfAve", arma::mean(idf)));
    feature_scores.push_back(std::make_pair("idfVar", arma::var(idf)));
    feature_scores.push_back(std::make_pair("idfMax", arma::max(idf)));

    feature_scores.push_back(std::make_pair("idfMaxMinRatio", arma::max(idf)));
    return;
}


void matIR::preretrieval::ctfRelated( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    // Compute IDF from qstats
    arma::vec collTermFreq = ( /
                                                stats.collectionLength);


    double score = arma::log(stats.collectionLength /
                                stats.getQueryStats().getQuerycTFs());

    feature_scores.push_back(std::make_pair("ictfAve", arma::mean(score)));
    feature_scores.push_back(std::make_pair("ictfMax", arma::max(score)));
    feature_scores.push_back(std::make_pair("ictfVar", arma::var(score)));
    return;
}


void matIR::preretrieval::simplified_query_clarity( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    // Get Probability( term | Query) --> Pml
    arma::vec pQueryLang = (stats.getQueryStats().getTermFreqInQuery() /
                                                stats.getQueryStats().queryLength);

    arma::vec pCollLang = (stats.getQueryStats().getQuerycTFs() /
                                                stats.collectionLength);


    double score = arma::sum(pQueryLang * arma::log( pQueryLang / pCollLang  ));
    feature_scores.push_back(std::make_pair("simplifiedQueryClarity", score));
    return;
}

void matIR::preretrieval::simple_features( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){


    std::map<string, int> queryText = stats.getQueryStats().getQueryText();

    // Average Query Length
    arma::vec tokenLength(queryText.size());
    int i = 0;
    map<std::string, int>::const_iterator iter;
    for (iter=queryText.begin(); iter != queryText.end(); ++iter)
        tokenLength(i++) = (iter->first).length();

    tokenLength = tokenLength * stats.getQueryStats().getTermFreqInQuery();
    feature_scores.push_back(std::make_pair("averageTokenLength",
                                                        arma::mean(tokenLength)));

    feature_scores.push_back(std::make_pair("tokenCount",
                                                stats.getQueryStats().queryLength));

    feature_scores.push_back(std::make_pair("uniqueTokenCount", queryText.size()));


    return;
}

