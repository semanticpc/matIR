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
// TermScoreFunctionFactory
//
// 19 August 2004 -- tds
//

#ifndef MATIR_TERMSCOREFUNCTIONFACTORY_HPP
#define MATIR_TERMSCOREFUNCTIONFACTORY_HPP

#include "matIR/TermScoreFunction.hpp"
#include <string>
#include <armadillo>
namespace matIR
{
  namespace scoring
  {
    class TermScoreFunctionFactory {
    public:
      static TermScoreFunction* get( const std::string& spec, arma::vec& p_term_given_Collection,
    		  double collectionLength_Md, int documentOccurrences = 0, int documentCount = 0 );
    };
  }
}

#endif // INDRI_TERMSCOREFUNCTIONFACTORY_HPP

