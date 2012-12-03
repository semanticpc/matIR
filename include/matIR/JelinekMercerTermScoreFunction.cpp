/*
 * File:   JelinekMercerTermScoreFunction.cpp
 * Author: semanticpc
 *
 * Created on December 3, 2012, 1:58 AM
 */

#ifndef MATIR_JELINEKMERCERTERMSCOREFUNCTION_CPP
#define	MATIR_JELINEKMERCERTERMSCOREFUNCTION_CPP


#include <math.h>
namespace matIR
{
    /*! Query processing, smoothing, and scoring classes. */
    namespace scoring
    {

    class JelinekMercerTermScoreFunction : public TermScoreFunction {
    private:
      double _lambda;
      double _backgroundLambda;
      double _oneLevelCollectionComponent;
      double _contextLambda;
      double _collectionLambda;
      double _documentLambda;
      double _foregroundLambda;

      arma::vec _prob_term_Collection;
      arma::vec _collectionComponent;

    public:
        JelinekMercerTermScoreFunction( arma::vec p_t_C, double collectionLambda, double documentLambda = 0.0  ) {

        _prob_term_Collection = p_t_C;
        _contextLambda = (1 - collectionLambda - documentLambda);

        _collectionLambda = collectionLambda;
        _documentLambda = documentLambda;
        _foregroundLambda = (1 - _collectionLambda);


        _collectionComponent = _collectionLambda * _prob_term_Collection;
        }



        //
        //             [                      occurrences                                             ]
        // score = log [ foregroundLambda * ---------------  + collectionLambda * collectionFrequency ]
        //             [                      contextSize                                             ]
        //



        arma::vec scoreOccurrence( arma::vec& termFrequency_Md, int length_Md ) {
            arma::vec zero_vector = arma::zeros(termFrequency_Md.size());
            arma::vec contextFrequency = length_Md ? (termFrequency_Md / double(length_Md)) : zero_vector;
            return arma::log( (_foregroundLambda * contextFrequency) + _collectionComponent );

        }



        // needs to be checks do not use
        arma::vec scoreOccurrence( arma::vec& occurrences, int contextSize, arma::vec& documentOccurrences, int documentLength ) {
            arma::vec zero_vector = arma::zeros(occurrences.size());
        arma::vec contextFrequency = contextSize ? occurrences / double(contextSize) : zero_vector;
        arma::vec documentFrequency = documentLength ? documentOccurrences / double(documentLength) : occurrences;
        return log( _contextLambda * contextFrequency + _documentLambda * documentFrequency + _collectionComponent );
        }
    };
    }
}

#endif	/* JELINEKMERCERTERMSCOREFUNCTION_CPP */

