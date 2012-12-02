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
// DirichletTermScoreFunction
//
// 26 January 2004 - tds
//

#ifndef MATIR_DIRICHLETTERMSCOREFUNCTION_HPP
#define MATIR_DIRICHLETTERMSCOREFUNCTION_HPP

#include <math.h>
namespace matIR
{
    /*! Query processing, smoothing, and scoring classes. */
    namespace scoring
    {
    
    class DirichletTermScoreFunction : public TermScoreFunction {
    private:
        double _mu;
        double _docmu;
        arma::vec _prob_term_Collection;
        arma::vec _mu_TIMES_prob_term_Collection;
        
    public:
        DirichletTermScoreFunction( double mu, arma::vec p_t_C, double docmu=-1.0 ) {
            _prob_term_Collection = p_t_C;
            _mu = mu;
            _mu_TIMES_prob_term_Collection = _mu * _prob_term_Collection;
            _docmu = docmu;
        }
        
        
        
        
        arma::vec scoreOccurrence( arma::vec& termFrequency_Md, int length_Md ) {
            return log(( termFrequency_Md + _mu_TIMES_prob_term_Collection ) / ( double(length_Md) + _mu ));
            
        }
        
        
        
        // needs to be checks do not use
        arma::vec scoreOccurrence( arma::vec& occurrences, int contextSize, arma::vec& documentOccurrences, int documentLength ) {
            //two level Dir Smoothing!
            //        tf_E + documentMu*P(t|D)
            //P(t|E)= ------------------------
            //         extentlen + documentMu
            //                 mu*P(t|C) + tf_D
            //where P(t|D)= ---------------------
            //                  doclen + mu
            // if the _docmu parameter is the default, do collection level
            // smoothing only.
            if (_docmu < 0)
                return scoreOccurrence(occurrences, contextSize);
            else {
                
                return log(((occurrences+_docmu)*(_mu_TIMES_prob_term_Collection+documentOccurrences)/(double(documentLength)+_mu))/(double(contextSize)+_docmu));
            }
        }
    };
    }
}

#endif // INDRI_DIRICHLETTERMSCOREFUNCTION_HPP
