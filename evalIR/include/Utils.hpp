/*
 * File:   utils.hpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:50 PM
 */

#ifndef Qrels_HPP
#define	Qrels_HPP
#include <armadillo>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

using namespace std;

// Data Structure

struct Document{
    string docid;
    string query;

    // Judgments Related
    set<int> subtopics;
    int grade;

    // Run Related
    double score;
    int rank;
    string runid;

    bool operator < (const Document& b) const{
        if (score == b.score){
            return docid < b.docid;
        }else
            return score > b.score;
    }

};

struct docid_comparison{
    //string docid;
    //find_doc(string docid) : docid(docid) {}
    bool operator () ( const Document& a, const Document& b ) const{
        return a.docid < b.docid;
    }
};

struct Qrels{
    string query;
    map<Document, int, docid_comparison> relDocs;
    set<Document> nonRelDocs;
    set<int> subtopics;
    arma::mat matrix;
};






map<string, vector<Document> > readRunFile(string runFileName){
    map<string, vector<Document> > run;

    std::ifstream runFile(runFileName.c_str(), ios_base::in);

    std::string q0;
    string query;

    string curQuery = "";
    bool doneParsing = false;

    // Parse First Line
    vector<Document> documents;
    Document d;
    runFile >> query >> q0 >> d.docid >> d.rank >> d.score >> d.runid;
    query = d.query;
    curQuery = d.query;

    documents.push_back(d);


    while (true) {
        while(curQuery == query){
            Document d;
            runFile >> d.query >> q0 >> d.docid >> d.rank >> d.score >> d.runid;
            query = d.query;
            documents.push_back(d);
            if( runFile.eof() ) {
                doneParsing = true;
                break;
            }

        }
        std::sort(documents.begin(), documents.end());
        run.insert(make_pair(curQuery,documents));

        if(doneParsing) break;
        curQuery = query;
        vector<Document> documents;
        documents.push_back(d);
        }
    return run;
}


map<string, Qrels> readDiversityQrelsFile(string qrelsFileName){
    std::ifstream qrelsFile(qrelsFileName.c_str(), ios_base::in);


    // Initialization

    map<string, Qrels> qrels;
    bool doneParsing = false;
    std::string curQuery = "";
    std::string docid;
    std::string query;
    int subtopic;
    int rel;


    // Initialize first query
    Qrels q;
    map<string, Document> relDocuments;


    // Read the first line
    qrelsFile >> query >> subtopic >> docid >> rel;
    Document d;
    d.query = query;
    d.docid = docid;
    d.subtopics.insert(subtopic);
    d.grade = rel;

    curQuery = query;

    // Add the document to the relevant or non-relevant list
    if(rel <= 0)
        q.nonRelDocs.insert(d);
    else
        relDocuments.insert(make_pair(docid, d));

    // Loop through add all other documents
    while (true) {

        // Loop through documents in the same query
        while(curQuery == query){

            // Read document
            qrelsFile >> query >> subtopic >> docid >> rel;
            Document d;
            d.query = query;
            d.docid = docid;
            d.grade = rel;


            // If the document was seen before simple add the subtopic
            map<string,Document>::iterator it = relDocuments.find(docid);
            if (it != relDocuments.end()){
                it->second.subtopics.insert(subtopic);
                if(q.subtopics.find(subtopic) == q.subtopics.end())
                    q.subtopics.insert(subtopic);
            }else{
                d.subtopics.insert(subtopic);
                if(q.subtopics.find(subtopic) == q.subtopics.end())
                    q.subtopics.insert(subtopic);
                if(rel <= 0)
                    q.nonRelDocs.insert(d);
                else
                    relDocuments.insert(make_pair(docid, d));
            }


            // Check for end of file
            if( qrelsFile.eof() ) {
                doneParsing = true;
                break;
            }

        }
        // Iterate through the map and add it to the relevant set
        q.matrix = arma::zeros(int(relDocuments.size()), int(*q.subtopics.rbegin()));

        map<string,Document>::iterator doc_it;
        int doc_index = 0;
        for(doc_it = relDocuments.begin(); doc_it != relDocuments.end(); ++doc_it ) {
            set<int>::iterator st_it;
            q.relDocs.insert(make_pair(doc_it->second,doc_index));
            for(st_it= doc_it->second.subtopics.begin();
                    st_it != doc_it->second.subtopics.end(); st_it++ )
                q.matrix(doc_index, (*st_it)-1)  = 1;

            doc_index++;
        }
        qrels.insert(make_pair(curQuery, q));


        if(doneParsing) break;

        // If Not End of file initialize another query
        Qrels q;
        map<string, Document> relDocuments;
        curQuery = query;

        // Read the next document
        qrelsFile >> query >> subtopic >> docid >> rel;
        Document d;
        d.query = query;
        d.docid = docid;
        d.subtopics.insert(subtopic);
        d.grade = rel;

        // Add the document to the relevant or non-relevant list
        if(q.subtopics.find(subtopic) == q.subtopics.end())
            q.subtopics.insert(subtopic);
        if(rel <= 0)
            q.nonRelDocs.insert(d);
        else
            relDocuments.insert(make_pair(docid, d));

    }
        return qrels;
}


arma::mat judge_diversity(vector<Document> run, Qrels qrels, int rank){
    arma::mat run_matrix = arma::zeros(std::min(rank, int(run.size())),int(*qrels.subtopics.rbegin()) );
    vector<Document>::iterator doc;
    int doc_index = 0;
    for(doc = run.begin(); doc != run.end(); ++doc ) {

        //map<Document, int>::iterator index_iter = std::find_if(qrels.relDocs.begin(),
          //      qrels.relDocs.end(), find_doc(doc->docid));
        map<Document, int>::iterator index_iter = qrels.relDocs.find(*doc);
        if(index_iter != qrels.relDocs.end()){

            run_matrix.row(doc_index) = qrels.matrix.row(int(index_iter->second));
        }
        doc_index++;
    }
    return run_matrix;
}

#endif	/* UTILS_HPP */

