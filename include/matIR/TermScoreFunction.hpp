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
// TermScoreFunction
//
// 23 January 2004 -- tds
//

#ifndef MATIR_TERMSCOREFUNCTION_HPP
#define MATIR_TERMSCOREFUNCTION_HPP
#include <armadillo>
namespace matIR
{
  namespace scoring
  {

    /*! Abstract base class for all term scoring and smoothing functions.
      See <a href="IndriParameters.html#rule">the <tt>rule</tt> parameter
      format</a> for a description of the rule parameter format. @see
      TermScoreFunctionFactory for a description of how to add a new scoring
      function.
    */
    class TermScoreFunction {
    public:
    	virtual ~TermScoreFunction(){};
      virtual arma::vec scoreOccurrence( arma::vec& occurrences, int contextLength ) = 0;
      virtual arma::vec scoreOccurrence( arma::vec& occurrences,
      int contextLength, arma::vec& documentOccurrences, int documentLength ) = 0;
    };
  }
}

#endif // INDRI_TERMSCOREFUNCTION_HPP

