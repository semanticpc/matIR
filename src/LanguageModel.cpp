/*
 * LanguageModel.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: semanticpc
 */


#include "matIR/LanguageModel.hpp"

#include <math.h>


//
// RelevanceModel
//

matIR::LanguageModel::LanguageModel(const std::string& smoothing,
                                    int termLimit,
                                    matIR::ResultStats& stats )
:
_smoothing(smoothing),
_termLimit(termLimit),
_stats(stats)
{
	_termScorer = matIR::scoring::TermScoreFunctionFactory::get( _smoothing,
                                                                _stats.ctfVector, _stats.collectionLength, 0 , 0 );
}

//
// ~RelevanceModel
//

matIR::LanguageModel::~LanguageModel() {
}

void matIR::LanguageModel::generateRelevanceModel(){

    int length_Md;
    arma::vec termFrequency_Md;
    arma::colvec lmscore = arma::zeros(_stats.tfMat.n_cols);

    for (int i=0; i < _stats.tfMat.n_rows ; i++) {
        // Need to transpose to be consistent with other vectors
        termFrequency_Md = arma::trans(_stats.tfMat.row(i));
        length_Md = (int)_stats.docLength[i];
        lmscore += (_stats.docScores[i] * arma::exp(_termScorer->scoreOccurrence(termFrequency_Md, length_Md)));
    }
    arma::uvec indices = arma::sort_index(lmscore, 1);
    _stats.sortGrams();
    _scoredTerms.clear();

    for(int i=0;i< std::min(_termLimit, (int)_stats.terms.size()); i++)
        _scoredTerms.push_back( std::make_pair(_stats.terms[indices(i)]->term,
                                                        (double)lmscore(indices(i))) );

}


void matIR::LanguageModel::generateLanguageModel(){


    // Need to transpose to be consistent with other vectors
	arma::vec termFrequency_Md = arma::sum(_stats.tfMat, 0).t();
        int length_Md;
        arma::colvec lmscore = arma::zeros(_stats.tfMat.n_cols);
        for (int i=0; i < _stats.tfMat.n_rows ; i++) {
        // Need to transpose to be consistent with other vectors
            termFrequency_Md = arma::trans(_stats.tfMat.row(i));
            length_Md = (int)_stats.docLength[i];
            lmscore += (arma::exp(_termScorer->scoreOccurrence(termFrequency_Md, length_Md)));
        }

	arma::uvec indices = arma::sort_index(lmscore,1);
	_stats.sortGrams();
	_scoredTerms.clear();
	for(int i=0;i< std::min(_termLimit, (int)_stats.terms.size()); i++)
		_scoredTerms.push_back( std::make_pair(_stats.terms[indices(i)]->term, (double)lmscore(indices(i))) );
}

const indri::utility::greedy_vector< std::pair< string, double > > matIR::LanguageModel::getScoredTerms() {
	return _scoredTerms;
}


void matIR::LanguageModel::setSmoothing(const std::string& smoothing){
    _smoothing.assign(smoothing);
    _termScorer = matIR::scoring::TermScoreFunctionFactory::get( _smoothing,
                                                                _stats.ctfVector, _stats.collectionLength, 0 , 0 );
}