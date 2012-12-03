
//
//  QueryStats.hpp
//  matIR
//
//  Created by Praveen Chandar on 11/30/12.
//  Copyright (c) 2012 Praveen Chandar. All rights reserved.
//
/*
 * QueryStats.hpp
 *
 *  Created on: Nov 29, 2012
 *      Author: semanticpc
 */



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



#ifndef MATIR_QUERYSTATS_HPP_
#define MATIR_QUERYSTATS_HPP_

#include <armadillo>

#include "indri/QueryEnvironment.hpp"


namespace matIR {

    class QueryStats {
    public:
        std::vector<string> queryString;
    private:
        std::map<string, int> _queryTokens;

        // Internal term ID of the query
        arma::uvec _query_matrix_indices;

        arma::vec _query_collectionFrequency;
        arma::vec _query_documentFrequency;

    public:

        QueryStats();
        //~QueryStats();


        void init(const std::string& , indri::api::QueryEnvironment&);
        const arma::vec getQueryDFs();
        const arma::vec getQuerycTFs();
        const arma::uvec getQueryTermIds();
        const std::map<std::string, int> getQueryTokens();
        arma::vec getTermFreqInQuery();
        void setTermIDs(arma::uvec);

    };

}

#endif // INDRI_DIRICHLETTERMSCOREFUNCTION_HPP// </editor-fold>// </editor-fold>
