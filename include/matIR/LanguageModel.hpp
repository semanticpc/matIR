/*
 * LanguageModel.hpp
 *
 *  Created on: Nov 30, 2012
 *      Author: semanticpc
 */

#ifndef LANGUAGEMODEL_HPP_
#define LANGUAGEMODEL_HPP_
#include <string>
#include <vector>
#include <armadillo>

#include "indri/greedy_vector"
#include "indri/QueryEnvironment.hpp"

#include "matIR/TermScoreFunctionFactory.hpp"
#include "matIR/ResultStats.hpp"

namespace matIR {
    class LanguageModel {
    public:



    private:
    	int _termLimit;
    	std::string _smoothing;
        matIR::ResultStats& _stats;
        matIR::scoring::TermScoreFunction* _termScorer;
    	indri::utility::greedy_vector< std::pair< string, double > > _scoredTerms;

    public:
        LanguageModel(  const std::string& smoothing,
                      int termLimit,
                      matIR::ResultStats& stats);
        ~LanguageModel();

        void generateRelevanceModel();
        void generateLanguageModel();
        void setSmoothing(const std::string& smoothing);
        const indri::utility::greedy_vector< std::pair< string, double > > getScoredTerms();

    };
}


#endif /* LANGUAGEMODEL_HPP_ */
