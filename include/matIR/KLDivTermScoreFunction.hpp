/*
 * File:   KLDivTermScoreFunction.hpp
 * Author: semanticpc
 *
 * Created on December 3, 2012, 1:02 AM
 */

#ifndef KLDIVTERMSCOREFUNCTION_HPP
#define	KLDIVTERMSCOREFUNCTION_HPP


#include <math.h>
namespace matIR
{
    /*! Query processing, smoothing, and scoring classes. */
    namespace scoring
    {

    class DirichletTermScoreFunction : public TermScoreFunction {
    private:
        arma::vec _prob_term_Collection;

    public:
        DirichletTermScoreFunction( arma::vec p_t_C ) {
            _prob_term_Collection = p_t_C;
        }




        arma::vec scoreOccurrence( arma::vec& termFrequency_Md, int length_Md ) {
            arma::vec p_t_d = (termFrequency_Md / length_Md);
            return (p_t_d * log( p_t_d / _prob_term_Collection ));

        }



        // needs to be checks do not use
        arma::vec scoreOccurrence( arma::vec& occurrences, int contextSize, arma::vec& documentOccurrences, int documentLength ) {
        }

    };
    }
}


#endif	/* KLDIVTERMSCOREFUNCTION_HPP */

