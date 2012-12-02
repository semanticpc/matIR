
#include "matIR/QueryStats.hpp"
#include "indri/QueryParserFactory.hpp"
#include "indri/RawScorerNodeExtractor.hpp"

matIR::QueryStats::QueryStats(){

}


void matIR::QueryStats::init(const std::string& query, indri::api::QueryEnvironment& environment)
{

    // Extract only the terms from the query and add to the vector
    indri::api::QueryParserWrapper *parser = indri::api::QueryParserFactory::get(query, "indri");
    indri::lang::ScoredExtentNode* rootNode = parser->query();
    indri::lang::RawScorerNodeExtractor extractor;
    rootNode->walk(extractor);
    std::vector<indri::lang::RawScorerNode*>& scorerNodes = extractor.getScorerNodes();

    queryLength = 0;
    for (int i = 0; i < scorerNodes.size(); i++){
        std::string qterm = scorerNodes[i]->queryText();
        queryLength += 1;
        if( _queryText.find(qterm) == _queryText.end() )
            _queryText.insert(make_pair( qterm, 1));
        else
            _queryText[qterm] += 1;
    }

    // Initialize vectors


    _query_collectionFrequency.set_size(_queryText.size());
    _query_documentFrequency.set_size(_queryText.size());



    // Now obtain the statistics
    int i = 0;
    map<std::string, int>::const_iterator iter;
    for (iter=_queryText.begin(); iter != _queryText.end(); ++iter) {
        std::string stem = environment.stemTerm(iter->first);
        _query_collectionFrequency(i) = (double) environment.stemCount(stem);
        _query_documentFrequency(i) = (double) environment.documentStemCount(stem);
        ++i;

    }

}

const arma::uvec matIR::QueryStats::getQueryTermIds() {
    return _query_matrix_indices;
}
const std::map<string, int> matIR::QueryStats::getQueryText() {
    return _queryText;
}
const arma::vec matIR::QueryStats::getQueryDFs() {
    return _query_documentFrequency;
}

const arma::vec matIR::QueryStats::getQuerycTFs() {
    return _query_collectionFrequency;
}// </editor-fold>

void matIR::QueryStats::setTermIDs(arma::uvec indices) {
    _query_matrix_indices = indices;
}


arma::vec matIR::QueryStats::getTermFreqInQuery() {
    arma::vec termFreqInQuery(_queryText.size());
    int i = 0;
    map<std::string, int>::const_iterator iter;
    for (iter=_queryText.begin(); iter != _queryText.end(); ++iter)
        termFreqInQuery(i++) = iter->second;
    return termFreqInQuery;
}




