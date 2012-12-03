/*
 * File:   Preretrieval.hpp
 * Author: semanticpc
 *
 * Created on December 1, 2012, 11:43 PM
 */

#ifndef PRERETRIEVAL_HPP
#define	PRERETRIEVAL_HPP

#include <vector>
#include "matIR/QueryStats.hpp"
#include "indri/greedy_vector"


namespace matIR {
    namespace preretrieval {

        // Functions for Pre-retrieval Features goes here

        void idfRelated( ResultStats& stats,
            indri::utility::greedy_vector< std::pair< string, double > >& features);

        void ctfRelated( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

        void simplified_query_clarity( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

        void simple_features( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

        void query_scope( ResultStats& stats , indri::api::QueryEnvironment& env,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);


        void pmi( ResultStats& stats , indri::api::QueryEnvironment& env,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);


        void collection_query_similarity( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

    }
}

#endif	/* PRERETRIEVAL_HPP */// </editor-fold>
