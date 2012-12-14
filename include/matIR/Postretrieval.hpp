/*
 * File:   Postretrieval.hpp
 * Author: semanticpc
 *
 * Created on December 1, 2012, 11:43 PM
 */

#ifndef POSTRETRIEVAL_HPP
#define	POSTRETRIEVAL_HPP
#include <vector>

#include "matIR/LanguageModel.hpp"
#include "indri/QueryEnvironment.hpp"
#include "matIR/QueryStats.hpp"
#include "indri/greedy_vector"


namespace matIR {
    namespace postretrieval {

        // Functions for Pre-retrieval Features goes here

        void query_clarity( std::string query,
        indri::api::QueryEnvironment& env,
        indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

        void weighted_info_gain( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

        void query_feedback( matIR::LanguageModel& lm_model,
                indri::api::QueryEnvironment& env,
                indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

        void NQC( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

        void retScore_related( ResultStats& stats ,
            indri::utility::greedy_vector< std::pair< string, double > >& feature_scores);

    }
}


#endif	/* POSTRETRIEVAL_HPP */

