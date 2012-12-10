//============================================================================
// Name        : matIR.cpp
// Author      : Praveen Chandar
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================
#include "indri/Parameters.hpp"
#include "indri/QueryEnvironment.hpp"
#include "indri/QueryParserFactory.hpp"
#include "matIR/ResultStats.hpp"
#include "matIR/LanguageModel.hpp"
#include "matIR/Preretrieval.hpp"
#include "matIR/Postretrieval.hpp"
#include <queue>


struct query_t {

    struct greater {

        bool operator() (query_t* one, query_t * two) {
            return one->index > two->index;
        }
    };

    query_t(int _index, std::string _number, const std::string& _text, const std::string &queryType, std::vector<std::string> workSet, std::vector<std::string> FBDocs) :
    index(_index),
    number(_number),
    text(_text), qType(queryType), workingSet(workSet), relFBDocs(FBDocs) {
    }

    query_t(int _index, std::string _number, const std::string & _text) :
    index(_index),
    number(_number),
    text(_text) {
    }

    std::string number;
    int index;
    std::string text;
    std::string qType;
    // working set to restrict retrieval
    std::vector<std::string> workingSet;
    // Rel fb docs
    std::vector<std::string> relFBDocs;
};

static bool copy_parameters_to_string_vector(std::vector<std::string>& vec, indri::api::Parameters p, const std::string& parameterName) {
    if (!p.exists(parameterName))
        return false;

    indri::api::Parameters slice = p[parameterName];

    for (size_t i = 0; i < slice.size(); i++) {
        vec.push_back(slice[i]);
    }

    return true;
}

void push_queue(std::queue< query_t* >& q, indri::api::Parameters& queries,
        int queryOffset) {

    for (size_t i = 0; i < queries.size(); i++) {
        std::string queryNumber;
        std::string queryText;
        std::string queryType = "indri";
        if (queries[i].exists("type"))
            queryType = (std::string) queries[i]["type"];
        if (queries[i].exists("text"))
            queryText = (std::string) queries[i]["text"];
        if (queries[i].exists("number")) {
            queryNumber = (std::string) queries[i]["number"];
        } else {
            int thisQuery = queryOffset + int(i);
            std::stringstream s;
            s << thisQuery;
            queryNumber = s.str();
        }
        if (queryText.size() == 0)
            queryText = (std::string) queries[i];

        // working set and RELFB docs go here.
        // working set to restrict retrieval
        std::vector<std::string> workingSet;
        // Rel fb docs
        std::vector<std::string> relFBDocs;
        copy_parameters_to_string_vector(workingSet, queries[i], "workingSetDocno");
        copy_parameters_to_string_vector(relFBDocs, queries[i], "feedbackDocno");

        q.push(new query_t(i, queryNumber, queryText, queryType, workingSet, relFBDocs));

    }
}

static void open_indexes(indri::api::QueryEnvironment& environment, indri::api::Parameters& param) {
    if (param.exists("index")) {
        indri::api::Parameters indexes = param["index"];

        for (unsigned int i = 0; i < indexes.size(); i++) {
            environment.addIndex(std::string(indexes[i]));
        }
    }

    if (param.exists("server")) {
        indri::api::Parameters servers = param["server"];

        for (unsigned int i = 0; i < servers.size(); i++) {
            environment.addServer(std::string(servers[i]));
        }
    }

    std::vector<std::string> smoothingRules;
    if (copy_parameters_to_string_vector(smoothingRules, param, "rule"))
        environment.setScoringRules(smoothingRules);
}


void preretrievalFeatures(matIR::ResultStats& rs_stats, indri::api::QueryEnvironment& env,
       indri::utility::greedy_vector< std::pair< string, double > >& features_scores){


    // Some simple features such as queryLength, average Token Length, etc
    matIR::preretrieval::simple_features(rs_stats, features_scores);

    // IDF Related Features
    matIR::preretrieval::idfRelated(rs_stats, features_scores);

    // cTF Related Features
    matIR::preretrieval::ctfRelated(rs_stats, features_scores);

    // Simplified Query Clarity
    matIR::preretrieval::simplified_query_clarity(rs_stats, features_scores);

    // Collection Query Similarity
    matIR::preretrieval::collection_query_similarity(rs_stats, features_scores);

    // Point-wise Mutual Information
    matIR::preretrieval::pmi(rs_stats, env, features_scores);

    // Query Scope
    //matIR::preretrieval::query_scope(rs_stats, env, features_scores);



}


void postretrievalFeatures(matIR::ResultStats& rs_stats,
        indri::api::QueryEnvironment& env,
        indri::api::Parameters& param,
        indri::utility::greedy_vector< std::pair< string, double > >& features_scores){
        int termLimit = 10;
        std::string rmSmoothing = "";
        matIR::LanguageModel lm(rmSmoothing, termLimit, rs_stats);
        //lm.generateRelevanceModel();

        matIR::postretrieval::query_clarity(lm, env, features_scores);
        //matIR::postretrieval::query_feedback(lm, env, features_scores);
        //matIR::postretrieval::NQC(rs_stats, features_scores);
        //matIR::postretrieval::retScore_related(rs_stats, features_scores);
        //matIR::postretrieval::weighted_info_gain(rs_stats, features_scores);

}

void generateFeatures(std::queue< query_t* >& queries, indri::api::Parameters& param,
        indri::api::QueryEnvironment& env) {

    int documents = (int) param[ "documents" ];
    int termLimit = 10;//(int) param[ "termLimit" ];

    std::string rmSmoothing = ""; // eventually, we should offer relevance model smoothing

    bool header = true;
    while (queries.size() > 0) {
        matIR::ResultStats stats(env, documents);
        query_t* q = queries.front();
        // Initialize Statistics for the current query
        if(q->workingSet.size() > 0){
            std::vector<lemur::api::DOCID_T> docids;
            docids = env.documentIDsFromMetadata("docno", q->workingSet);
            stats.init(q->text, docids);
        }else{
            stats.init(q->text);
        }

        //indri::utility::greedy_vector< std::pair< string, double > > model = lm.getScoredTerms();
        indri::utility::greedy_vector< std::pair< string, double > > features_scores;

        // Generate Pre-retrieval features
        preretrievalFeatures(stats, env,  features_scores);

        //postretrievalFeatures(stats, env, param, features_scores);

        // Generate Post-retrieval features
        // Generate Document features

        // Print query related features

        if(header){
            header = false;
            cout << "topic";
            for(int i=0;i < features_scores.size(); i++)
                cout << "," <<  features_scores[i].first;
            cout << endl;
        }
        cout << q->number;
        for(int i=0;i < features_scores.size(); i++)
            cout <<  "," <<  features_scores[i].second;

        queries.pop();
    }

}

static void usage(indri::api::Parameters param) {
    if (!param.exists("query") || !(param.exists("index") ||
            param.exists("server")) || !param.exists("documents")) {
        std::cerr << "featureVector usage: " << std::endl
                << "   featureVector -query=myquery -index=myindex -documents=10 -maxGrams=2" << std::endl
                << "     myquery: a valid Indri query (be sure to use quotes around it if there are spaces in it)" << std::endl
                << "     myindex: a valid Indri index" << std::endl
                << "     documents: the number of documents to use to build the language model" << std::endl
                << "     termLimit: the number of terms in the language model to  output" << std::endl
                << "     maxGrams (optional): maximum length (in words) of phrases to be added to the model, default is 1 (unigram)" << std::endl;
        exit(-1);
    }
}

int main(int argc, char * argv[]) {
    try {
        indri::api::Parameters& param = indri::api::Parameters::instance();
        param.loadCommandLine(argc, argv);
        usage(param);
        if (param.get("version", 0)) {
            std::cout << INDRI_DISTRIBUTION << std::endl;
        }

        if (!param.exists("query"))
            LEMUR_THROW(LEMUR_MISSING_PARAMETER_ERROR, "Must specify at least one query.");

        if (!param.exists("index") && !param.exists("server"))
            LEMUR_THROW(LEMUR_MISSING_PARAMETER_ERROR, "Must specify a server or index to query against.");

        if (param.exists("baseline") && param.exists("rule"))
            LEMUR_THROW(LEMUR_BAD_PARAMETER_ERROR, "Smoothing rules may not be specified when running a baseline.");

        int threadCount = param.get("threads", 1);
        std::queue< query_t* > queries;


        // push all queries onto a queue
        indri::api::Parameters parameterQueries = param[ "query" ];
        int queryOffset = param.get("queryOffset", 0);
        push_queue(queries, parameterQueries, queryOffset);
        int queryCount = (int) queries.size();

        indri::api::QueryEnvironment environment;
        open_indexes(environment, param);

        generateFeatures(queries, param, environment);
    } catch (lemur::api::Exception& e) {
        LEMUR_ABORT(e);
    } catch (...) {
        std::cout << "Caught unhandled exception" << std::endl;
        return -1;
    }
    return 0;
}// </editor-fold>
