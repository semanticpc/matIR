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

#include "matIR/TermScoreFunctionFactory.hpp"
#include "matIR/DirichletTermScoreFunction.hpp"
#include "matIR/JelinekMercerTermScoreFunction.cpp"
#include "matIR/SimpleTermScoreFunction.hpp"
#include "indri/Parameters.hpp"


static void termscorefunctionfactory_parse( indri::api::Parameters& converted, const std::string& spec );


matIR::scoring::TermScoreFunction* matIR::scoring::TermScoreFunctionFactory::get(
                                                                                 const std::string& stringSpec, arma::vec& collectionFrequency,
                                                                                 double contextSize, int documentOccurrences, int documentCount ) {

    indri::api::Parameters spec;
    termscorefunctionfactory_parse( spec, stringSpec );
    std::string method = spec.get( "method", "" );

    arma::vec p_t_C = collectionFrequency /contextSize;

    if( method == "dirichlet" || method == "d" || method == "dir" ) {
        double mu = spec.get( "mu", 2500 );
        double docmu=spec.get("documentMu",-1.0); // default is no doc-level smoothing
        return new matIR::scoring::DirichletTermScoreFunction( mu, p_t_C, docmu );
    } else if( method == "linear" || method == "jm" || method == "jelinek-mercer" ) {
        // jelinek-mercer -- can take parameters collectionLambda (or just lambda) and documentLambda
        double documentLambda = spec.get( "documentLambda", 0.0 );
        double collectionLambda;

        if( spec.exists( "collectionLambda" ) )
            collectionLambda = spec.get( "collectionLambda", 0.5 );
        else
            collectionLambda = spec.get( "lambda", 0.5 );

        return new matIR::scoring::JelinekMercerTermScoreFunction( p_t_C, collectionLambda, documentLambda );
    }else {
        return new matIR::scoring::SimpleTermScoreFunction();
    }


}

static void termscorefunctionfactory_parse( indri::api::Parameters& converted, const std::string& spec ) {
    int nextComma = 0;
    int nextColon = 0;
    int  location = 0;

    for( location = 0; location < spec.length(); ) {
        nextComma = spec.find( ',', location );
        nextColon = spec.find( ':', location );

        std::string key = spec.substr( location, nextColon-location );
        std::string value = spec.substr( nextColon+1, nextComma-nextColon-1 );

        converted.set( key, value );

        if( nextComma > 0 )
            location = nextComma+1;
        else
            location = spec.size();
    }
}

