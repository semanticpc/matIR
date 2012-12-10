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
#include <math.h>
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

    feature_scores.push_back(std::make_pair("idfMaxMinRatio", arma::min(idf)/arma::max(idf)));
    return;
}


void matIR::preretrieval::ctfRelated( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    // Compute IDF from qstats
    arma::vec ictf = arma::log(stats.collectionLength /
                                stats.getQueryStats().getQuerycTFs());

    feature_scores.push_back(std::make_pair("ictfAve", arma::mean(ictf)));
    feature_scores.push_back(std::make_pair("ictfMax", arma::max(ictf)));
    feature_scores.push_back(std::make_pair("ictfVar", arma::var(ictf)));
    return;
}


void matIR::preretrieval::query_scope( ResultStats& stats ,indri::api::QueryEnvironment& env,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    // Query Scope

    // Make an Boolean OR query and obtain the number of documents containing
    //   the query expression
    std::map<string, int> queryTokens = stats.getQueryStats().getQueryTokens();

    // Average Query Length
    std:string queryExpression = " #combine(";
    int i = 0;
    map<std::string, int>::const_iterator iter;
    for (iter=queryTokens.begin(); iter != queryTokens.end(); ++iter){
        queryExpression = queryExpression + " #1( " + iter->first + ") " ;
        //break;
    }
    queryExpression += " )";


    double nq = env.expressionCount(queryExpression);
    double queryScore = -1 * log(nq / stats.documentCount);

    feature_scores.push_back(std::make_pair("queryScope", queryScore));
    return;
}


void matIR::preretrieval::simplified_query_clarity( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){

    // Get Probability( term | Query) --> Pml
    arma::vec pQueryLang = (stats.getQueryStats().getTermFreqInQuery() /
                                        stats.getQueryStats().queryString.size());

    arma::vec pCollLang = (stats.getQueryStats().getQuerycTFs() /
                                                stats.collectionLength);


    double score = arma::sum(pQueryLang % arma::log( pQueryLang / pCollLang  ));
    feature_scores.push_back(std::make_pair("simplifiedQueryClarity", score));
    return;
}

void matIR::preretrieval::simple_features( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){


    std::map<string, int> queryTokens = stats.getQueryStats().getQueryTokens();
    // Average Query Length
    arma::vec tokenLength(queryTokens.size());
    int i = 0;
    map<std::string, int>::const_iterator iter;
    for (iter=queryTokens.begin(); iter != queryTokens.end(); ++iter)
        tokenLength(i++) = (iter->first).length();

    tokenLength = tokenLength * stats.getQueryStats().getTermFreqInQuery();
    feature_scores.push_back(std::make_pair("averageTokenLength",
                                                        arma::mean(tokenLength)));

    feature_scores.push_back(std::make_pair("tokenCount",
                                                stats.getQueryStats().queryString.size()));

    feature_scores.push_back(std::make_pair("uniqueTokenCount", queryTokens.size()));


    return;
}



void matIR::preretrieval::collection_query_similarity( ResultStats& stats ,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){


    arma::vec idf = arma::log((stats.documentCount + 1) /
                                        (stats.getQueryStats().getQueryDFs() + 0.5));


    // % is scalar multiplication of two vectors in Aramidillo Library
    arma::vec SCQ = (arma::log(stats.getQueryStats().getQuerycTFs()) + 1) % idf;

    feature_scores.push_back(std::make_pair("SCQAvg", arma::mean(SCQ)));
    feature_scores.push_back(std::make_pair("SCQSum", arma::sum(SCQ)));
    feature_scores.push_back(std::make_pair("SCQMax", arma::max(SCQ)));
    return;
}


void matIR::preretrieval::pmi( ResultStats& stats , indri::api::QueryEnvironment& env,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores){



    // Make an Boolean OR query and obtain the number of documents containing
    //   the query expression
    std::vector<string> queryString = stats.getQueryStats().queryString;
    arma::vec query_ctf = stats.getQueryStats().getQuerycTFs();

    if(queryString.size() == 1){
        feature_scores.push_back(std::make_pair("avgPMI", 0));
        return;
    }


    // Average Query Length
    arma::vec pmi_score = arma::zeros(queryString.size() - 1);

    for(int i =0; (i + 1) < queryString.size(); i++){
        std:string queryExpression = " #1( " + queryString[i] + " " +
                                        queryString[i + 1] + " )";

        double p_t1_t2_D = env.expressionCount(queryExpression) / stats.collectionLength;
        double p_t1_D = query_ctf[i] / stats.collectionLength;
        double p_t2_D = query_ctf[i+1] / stats.collectionLength;
        if(p_t1_D == 0 || p_t2_D == 0 || p_t1_t2_D  == 0 )
            pmi_score(i) = 0;
        else
            pmi_score(i) = log(p_t1_t2_D / (p_t1_D * p_t2_D));
    }

    feature_scores.push_back(std::make_pair("avgPMI", arma::mean(pmi_score)));

    return;
}