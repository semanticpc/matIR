
/*==========================================================================
 * Copyright (c) 2004 University of Massachusetts.  All Rights Reserved.
 *
 * Use of the Lemur Toolkit for Language Modeling and Information Retrieval
 * is subject to the terms of the software license set forth in the LICENSE
 * file included with this software, and also available at
 * http://www.lemurproject.org/license.html
 *
 *==========================================================================
 */

    //
    // relevancemodel
    //
    // 23 June 2005 -- tds
    //
#include "indri/QueryEnvironment.hpp"
#include "matIR/ResultStats.hpp"
#include "matIR/QueryStats.hpp"
#include <math.h>


    //
    // RelevanceModel
    //

matIR::ResultStats::ResultStats(
                                indri::api::QueryEnvironment& environment,
                                int documents )
:
_environment(environment),
_documents(documents)
{

	documentCount = _environment.documentCount();
	collectionLength = _environment.termCount();


}

    //
    // ~RelevanceModel
    //

matIR::ResultStats::~ResultStats() {
	documentCount = -1;
	collectionLength = -1;
	tfMat.reset();
	dfVector.reset();
	ctfVector.reset();
	docLength.reset();
	docScores.reset();
}

    //
    // getQueryResults
    //

const std::vector<indri::api::ScoredExtentResult>& matIR::ResultStats::getQueryResults() const {
    return _results;
}

    //
    // _extractDocuments
    //

void matIR::ResultStats::_extractDocuments() {
    for( size_t i=0; i<_results.size(); i++ ) {
        _documentIDs.push_back( _results[i].document );
    }
}

    //
    // _countGrams
    //
    // Builds a hash table of grams, and counts the times that each
    // gram occurs in each query result.
    //

bool isValid(const string & word)
{
    size_t length = word.size();
    const char * chArray = word.c_str();
    size_t pos = 0;

    while (pos < length)
        {
        if(isalnum((unsigned char)*(chArray+pos)) == 0)
            {
            return false;
            }
        pos ++;
        }
    return true;
}


void matIR::ResultStats::_countGrams() {
    docLength.reset();
    docScores.reset();
    docLength.set_size(_results.size());
    docScores.set_size(_results.size());
        // for each query result
    for( size_t i=0; i<_results.size(); i++ ) {
            // run through the text, extracting n-grams
        indri::api::ScoredExtentResult& result = _results[i];
        indri::api::DocumentVector* v = _vectors[i];
        std::vector<int>& positions = v->positions();
        std::vector<std::string>& stems = v->stems();
        if (result.end == 0) result.end = positions.size();

            // for each word position in the text
        for( int j = result.begin; j < result.end; j++ ) {
                //int maxGram = std::min( _maxGrams, result.end - j );

                // extract every possible n-gram that starts at this position
                // up to _maxGrams in length

            GramCounts* newCounts = new GramCounts;
            bool containsOOV = false;

                // build the gram

            if( positions[ j ] == 0 || (! isValid(stems[ positions[ j ] ])) ) {
                containsOOV = true;
                continue;
            }

            newCounts->gram.term =  stems[ positions[ j ] ] ;


            if( containsOOV ) {
                    // if this contanied OOV, all larger n-grams
                    // starting at this point also will
                delete newCounts;
                break;
            }

            GramCounts** gramCounts = 0;
            gramCounts = _gramTable.find( &newCounts->gram );

            if( gramCounts == 0 ) {
                _gramTable.insert( &newCounts->gram, newCounts );
                gramCounts = &newCounts;
            } else {
                delete newCounts;
            }

            if( (*gramCounts)->counts.size() && (*gramCounts)->counts.back().first == i ) {
                    // we already have some counts going for this query result, so just add this one
                (*gramCounts)->counts.back().second++;
            } else {
                    // no counts yet in this document, so add an entry
                (*gramCounts)->counts.push_back( std::make_pair( i, 1 ) );
            }

        }
        // Store the values
        docScores(i) = _results[i].score;
        docLength(i) = _environment.documentLength(_results[i].document);

    }
}


void matIR::ResultStats::_buildStats() {
        // Initialize the Hash Map for the vocabulary
        //   <key, value> -> <termString, matrixIndex>
        //vocab.clear();
	HGram::iterator iter;

	tfMat = arma::zeros<arma::mat>(_results.size(), _gramTable.size());


        // Initialize the
	dfVector.set_size(_gramTable.size());
	ctfVector.set_size(_gramTable.size());

	int tmpTermID = -1;
	for( iter = _gramTable.begin(); iter != _gramTable.end(); iter++ ) {
            double gramCount = 0;
            ++tmpTermID;
            Gram* gram = *iter->first;
            GramCounts* gramCounts = *iter->second;
            gram->internal_termID = tmpTermID;

            ctfVector(tmpTermID) = _environment.stemCount(gram->term);
            dfVector(tmpTermID) =  _environment.documentStemCount(gram->term);
            size_t c;
            size_t r;
            for( r = 0, c = 0; r < _results.size() && c < gramCounts->counts.size(); r++ ) {
                if( gramCounts->counts[c].first == r ) {
                // we have counts for this result
                    tfMat(r, tmpTermID) = gramCounts->counts[c].second;
                    c++;
                }
            }
        }

    return;

}


arma::mat matIR::ResultStats::getTFIDFMatrix() {
	arma::mat tfidfMat = tfMat;
        arma::vec idf = arma::log((documentCount - dfVector + 0.5) / (dfVector + 0.5));
	tfidfMat.each_row() %= idf;
	return tfidfMat;
}

const arma::mat& matIR::ResultStats::getTFMatrix() const {
	return tfMat;
}
    // In:  log(x1) log(x2) ... log(xN)
    // Out: x1/sum, x2/sum, ... xN/sum
    //
    // Extra care is taken to make sure we don't overflow
    // machine precision when taking exp (log x)
    // This is done by adding a constant K which cancels out
    // Right now K is set to maximally preserve the highest value
    // but could be altered to a min or average, or whatever...

static void _logtoposterior(std::vector<indri::api::ScoredExtentResult> &res) {
    if (res.size() == 0) return;
    std::vector<indri::api::ScoredExtentResult>::iterator iter;
    iter = res.begin();
    double K = (*iter).score;
        // first is max
    double sum=0;

    for (iter = res.begin(); iter != res.end(); iter++) {
        sum += (*iter).score=exp(K+(*iter).score);
    }
    for (iter = res.begin(); iter != res.end(); iter++) {
        (*iter).score/=sum;
    }
}

matIR::QueryStats& matIR::ResultStats::getQueryStats(){
    return _queryStatistics;
}

void matIR::ResultStats::sortGrams(){

        // copy grams into a _grams vector
    HGram::iterator iter;

    terms.clear();

    for( iter = _gramTable.begin(); iter != _gramTable.end(); iter++ ) {
        terms.push_back( *(iter->first) );
    }

    std::sort( terms.begin(), terms.end(), Gram::weight_greater() );
}

void matIR::ResultStats::_setQueryStatistics(const std::string& query){
        // Update query statistics
        _queryStatistics.init(query, _environment);
        std::map<string, int> queryTokens = _queryStatistics.getQueryTokens();
        // Temporary hack
        arma::uvec indices(queryTokens.size());
        GramCounts* newCounts = new GramCounts;

        int i = 0;
        map<string, int>::const_iterator iter;
        for (iter=queryTokens.begin(); iter != queryTokens.end(); ++iter) {
            newCounts->gram.term = iter->first;
            GramCounts** gramCounts = _gramTable.find( &newCounts->gram );

            if( gramCounts == 0 ) {
                newCounts->gram.internal_termID = tfMat.size() + 1;
                _gramTable.insert( &newCounts->gram, newCounts );
                gramCounts = &newCounts;
                // Add Dummy Column to the tfMatrix
                indices(i) = -1;

            }else{
            indices(i) = (*gramCounts)->gram.internal_termID;
            }
            ++i;
        }
        _queryStatistics.setTermIDs(indices);
}
    //
    // generate
    //

void matIR::ResultStats::init( const std::string& query ) {
    try {

        // run the query, get the document vectors
        _results = _environment.runQuery( query, _documents );
        _logtoposterior(_results);
        _extractDocuments();

        _vectors = _environment.documentVectors( _documentIDs );
        _countGrams();
        _buildStats();
        _setQueryStatistics(query);
        for (unsigned int i = 0; i < _vectors.size(); i++)
            delete _vectors[i];


    } catch( lemur::api::Exception& e ) {
        LEMUR_RETHROW( e, "Couldn't generate relevance model for '" + query + "' because: " );
    }
}

    //
    // generate
    //

void matIR::ResultStats::init( const std::string& query, const std::vector<lemur::api::DOCID_T>& documentSet  ) {
    try {
        _results = _environment.runQuery( query, documentSet, documentSet.size());
        _logtoposterior(_results);
        _extractDocuments();
        _vectors = _environment.documentVectors( _documentIDs );
        _countGrams();
        _buildStats();

        _setQueryStatistics(query);
        for (unsigned int i = 0; i < _vectors.size(); i++)
            delete _vectors[i];


    } catch( lemur::api::Exception& e ) {
        LEMUR_RETHROW( e, "Couldn't generate relevance model for '" + query + "' because: " );
    }
}

const std::vector<lemur::api::DOCID_T>& matIR::ResultStats::getDocumentIDs() const{
    return _documentIDs;
}

