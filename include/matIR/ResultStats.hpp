/*
 * QueryStats.hpp
 *
 *  Created on: Nov 29, 2012
 *      Author: semanticpc
 */

#ifndef MATIR_RESULTSTATS_HPP_
#define MATIR_RESULTSTATS_HPP_
#include <string>
#include <vector>
#include <armadillo>

#include "indri/HashTable.hpp"
#include "indri/greedy_vector"
#include "indri/QueryEnvironment.hpp"

#include "matIR/QueryStats.hpp"


namespace matIR {
class ResultStats {
public:
	struct Gram {
		//std::vector<std::string> terms;
		std::string term;
		int internal_termID;

		struct hash {
			int operator() ( const Gram* one ) const {
				indri::utility::GenericHash<const char*> h;
				int accumulator = 0;

				//for( size_t i=0; i<one->terms.size(); i++ ) {
					accumulator *= 7;
					accumulator += h( one->term.c_str() );
					//}

					return accumulator;
			}
		};


		struct weight_greater {
			bool operator() ( const Gram* o, const Gram* t ) const {
				return t->internal_termID > o->internal_termID;
			}
		};

		struct string_comparator {
			int operator() ( const Gram* o, const Gram* t ) const {
				const Gram& one = *o;
				const Gram& two = *t;

				/*if( one.terms.size() != two.terms.size() ) {
    	                    if( one.terms.size() < two.terms.size() ) {
    	                      return 1;
    	                    } else {
    	                      return -1;
    	                    }
    	                  }*/

				//for( size_t i=0; i<one.terms.size(); i++ ) {
				const std::string& oneString = one.term;//s[i];
				const std::string& twoString = two.term;//s[i];

				if( oneString != twoString ) {
					if( oneString < twoString )
						return -1;
					else
						return 1;
				}
				//}

				return 0;
			}
		};
	};
	arma::mat tfMat;
	arma::vec dfVector;
	arma::vec ctfVector;

	arma::vec docScores;
	arma::vec docLength;


	// Index Statistics
	int documentCount;
	int collectionLength;


	std::vector<Gram*> terms;



        struct GramCounts {
             Gram gram;
             indri::utility::greedy_vector< std::pair< int, int > > counts;
         };

	typedef indri::utility::HashTable< Gram*, GramCounts*, Gram::hash, Gram::string_comparator > HGram;
	HGram _gramTable;

private:




        matIR::QueryStats _queryStatistics;
	indri::api::QueryEnvironment& _environment;
	int _documents;

	// Datastructure to store the results retrieved for the query
	std::vector<indri::api::ScoredExtentResult> _results;
	std::vector<lemur::api::DOCID_T> _documentIDs;
	std::vector<indri::api::DocumentVector*> _vectors;

	void _countGrams();
	void _buildStats();
	void _extractDocuments();
        void _setQueryStatistics(const std::string& query);

public:
	ResultStats( indri::api::QueryEnvironment& environment,
			int documents );
	~ResultStats();


	void init( const std::string& query );
	// generate from an existing result set
	void init( const std::string &query , const std::vector<lemur::api::DOCID_T>& docids );
	const std::vector<indri::api::ScoredExtentResult>& getQueryResults() const;
        matIR::QueryStats& getQueryStats();
	arma::mat getTFIDFMatrix();
	void sortGrams();
        const std::vector<lemur::api::DOCID_T>& getDocumentIDs() const;
	const arma::mat& getTFMatrix() const;
};
}


#endif /* MATIR_RESULTSTATS_HPP_ */
