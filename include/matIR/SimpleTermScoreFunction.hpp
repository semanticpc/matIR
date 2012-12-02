//
//  SimpleTermScoreFunction.cpp
//  matIR
//
//  Created by Praveen Chandar on 11/30/12.
//  Copyright (c) 2012 Praveen Chandar. All rights reserved.
//



//
// SimpleTermScoreFunction
//
// 26 January 2004 - tds
//

#ifndef MATIR_SIMPLETERMSCOREFUNCTION_HPP
#define MATIR_SIMPLETERMSCOREFUNCTION_HPP

namespace matIR
{
    namespace scoring    {
        
        class SimpleTermScoreFunction : public TermScoreFunction {
        private:

            
        public:
            SimpleTermScoreFunction( ) {}
            
            arma::vec scoreOccurrence( arma::vec& termFrequency_Md, int length_Md ) {
                return (termFrequency_Md ) / ( double(length_Md));
            }
            
            // needs to be checks do not use
            arma::vec scoreOccurrence( arma::vec& occurrences, int contextSize, arma::vec& documentOccurrences, int documentLength ) {
                return scoreOccurrence(occurrences, contextSize);
            }
        };
    }
}

#endif // MATIR_SIMPLETERMSCOREFUNCTION_HPP
