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
#include <regex.h>
#include <string>

using namespace std;

// Data Structure

struct Document{
    string docid;
    int query;

    // Judgments Related
    set<int> subtopics;
    int grade;

    // Run Related
    double score;
    int rank;
    string runid;

    bool operator < (const Document& b) const{
        if (rank == b.rank){
            return docid < b.docid;
        }else
            return rank < b.rank;
    }

};

struct docid_comparison{
    //string docid;
    //find_doc(string docid) : docid(docid) {}
    bool operator () ( const Document& a, const Document& b ) const{
        return a.docid < b.docid;
    }
};

struct strdocid_comparison{
    //string docid;
    //find_doc(string docid) : docid(docid) {}
    bool operator () ( const string& a, const string& b ) const{
        return a < b;
    }
};

struct Qrels{
    int query;
    int numOfRelDocs;
    map<Document, int, docid_comparison> relDocs;
    set<Document> nonRelDocs;
    set<int> subtopics;
    arma::mat matrix;
    arma::vec subtopicImportance;

};

struct Profiles{
    int query;
    map<int, double> subtopic_importance;
};


// Some Utility Functions
static bool checkForDouble(std::string const& s) {
  std::istringstream ss(s);
  double d;
  return (ss >> d) && (ss >> std::ws).eof();
}


// Parsing Functions
static map<int, vector<Document> > readRunFile(string runFileName){
    map<int, vector<Document> > run;

    std::ifstream runFile(runFileName.c_str(), ios_base::in);

    std::string q0;
    int query;

    int curQuery;
    bool doneParsing = false;

    // Parse First Line
    vector<Document> documents;
    Document d;
    runFile >> d.query >> q0 >> d.docid >> d.rank >> d.score >> d.runid;
    query = d.query;
    curQuery = d.query;
    documents.push_back(d);

    while(curQuery == query){
        Document d;
        runFile >> d.query >> q0 >> d.docid >> d.rank >> d.score >> d.runid;
        query = d.query;
        if(curQuery != query) {
            std::sort(documents.begin(), documents.end());
            run.insert(make_pair(curQuery,documents));
            curQuery = query;
            documents.clear();
        }

        documents.push_back(d);
        if( runFile.eof()) {
            doneParsing = true;
            break;
        }
    }

    std::sort(documents.begin(), documents.end());
    run.insert(make_pair(curQuery,documents));
    return run;
}


static map<int, Qrels> readDiversityQrelsFile(string qrelsFileName){
    std::ifstream qrelsFile(qrelsFileName.c_str(), ios_base::in);


    // Initialization

    map<int, Qrels> qrels;
    bool doneParsing = false;
    int curQuery;
    std::string docid;
    int query;
    int subtopic;
    int rel;


    // Initialize first query
    Qrels q;
    q.relDocs.clear();
    map<string, Document, strdocid_comparison> relDocuments;


    // Read the first line
    qrelsFile >> query >> subtopic >> docid >> rel;
    Document d;
    d.query = query;
    d.docid = docid;
    d.subtopics.insert(subtopic);
    d.grade = rel;

    curQuery = query;
    q.query = curQuery;
    q.subtopics.insert(subtopic);

    // Add the document to the relevant or non-relevant list
    if(rel <= 0)
        q.nonRelDocs.insert(d);
    else
        relDocuments.insert(make_pair(docid, d));

    // Loop through add all other documents
    while(curQuery == query){

        // Read document
        qrelsFile >> query >> subtopic >> docid >> rel;
        Document d;
        d.query = query;
        d.docid = docid;
        d.grade = rel;

        if(curQuery != query) {
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
            q.numOfRelDocs = int(q.relDocs.size());
            q.subtopicImportance = arma::ones(q.subtopics.size());
            qrels.insert(make_pair(curQuery, q));


            curQuery = query;
            q.matrix.reset();
            q.subtopics.clear();
            q.nonRelDocs.clear();
            q.relDocs.clear();
            q.query = curQuery;
            q.numOfRelDocs = 0;
            q.subtopicImportance.clear();
            relDocuments.clear();


        }

        // If the document was seen before simple add the subtopic
        map<string,Document>::iterator it = relDocuments.find(docid);
        if (it != relDocuments.end()){
            if(rel > 0 && subtopic != 0){
                it->second.subtopics.insert(subtopic);
                if(q.subtopics.find(subtopic) == q.subtopics.end())
                        q.subtopics.insert(subtopic);
            }
        }else{


            if(rel > 0){
                d.subtopics.insert(subtopic);
                relDocuments.insert(make_pair(docid, d));
                if(q.subtopics.find(subtopic) == q.subtopics.end())
                    q.subtopics.insert(subtopic);
            }
            else
                q.nonRelDocs.insert(d);

        }

        // Check for end of file
         if( qrelsFile.eof() ) {
             doneParsing = true;
             break;
         }
    }

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
    q.numOfRelDocs = int(q.relDocs.size());
    q.subtopicImportance = arma::ones(q.subtopics.size());
    qrels.insert(make_pair(curQuery, q));



    return qrels;
}


static arma::mat judge_diversity(vector<Document> &run, Qrels &qrels, int rank){

    arma::mat run_matrix = arma::zeros(rank ,int(*qrels.subtopics.rbegin()) );
    vector<Document>::iterator doc;
    int doc_index = 0;
    map<Document, int>::iterator doc_iter = qrels.relDocs.begin();

    for(doc = run.begin(); doc != run.end(); ++doc ) {
        if(doc_index >= rank) break;
        map<Document, int>::iterator index_iter = qrels.relDocs.find(*doc);
        if(index_iter != qrels.relDocs.end()){
            run_matrix.row(doc_index) = qrels.matrix.row(int(index_iter->second));
        }
        doc_index++;
    }
    return run_matrix;
}

static Qrels update_SubtopicImportance(Qrels& qrels, map<int, double> &_subtopicImportance){
    Qrels new_qrels = qrels;
    new_qrels.subtopicImportance = arma::zeros(qrels.subtopics.size());
    arma::vec multiplier = arma::zeros(qrels.subtopics.size());
    int i =0 ;
    set<int>::iterator st_it;
    map<int, double>::iterator find_st;
    for(st_it= qrels.subtopics.begin(); st_it != qrels.subtopics.end(); st_it++ ){
        find_st = _subtopicImportance.find(*st_it);
        if(find_st != _subtopicImportance.end()){
            new_qrels.subtopicImportance(i) = find_st->second;
            multiplier(i) = 1;
            i++;
        }
    }
    new_qrels.matrix.each_row() %= multiplier;
    new_qrels.numOfRelDocs = arma::sum(arma::sum(new_qrels.matrix, 1) >
                                        arma::zeros(new_qrels.matrix.n_rows));
    return new_qrels;
}


static Qrels update_SubtopicImportance_binary(Qrels& qrels, map<int, double> &_subtopicImportance){
    Qrels new_qrels = qrels;
    new_qrels.subtopicImportance = arma::zeros(qrels.matrix.n_cols);


    set<int>::iterator find_st;
    map<int, double>::iterator st_it;
    for(st_it= _subtopicImportance.begin(); st_it != _subtopicImportance.end(); st_it++ ){
        find_st = qrels.subtopics.find(st_it->first);
        if(find_st != qrels.subtopics.end()){
            if(st_it->second > 0)
                new_qrels.subtopicImportance(st_it->first - 1) = 1;
        }

    }

    new_qrels.matrix.each_row() %= new_qrels.subtopicImportance;
    new_qrels.numOfRelDocs = arma::sum(arma::sum(new_qrels.matrix, 1) >
                                        arma::zeros(new_qrels.matrix.n_rows));
    return new_qrels;
}



// User Profile Related Functions
static map<int, map<int, Profiles* > > read_userProfiles(string filename){
    std::ifstream  data(filename.c_str(), ios_base::in);
    map<int, map<int, Profiles* > > userProfiles;
    map<int, map<int, Profiles* > >::iterator iter;
    map<int, Profiles* >::iterator userProfiles_iter;
    std::string line;
    while(std::getline(data,line)){
        std::stringstream  lineStream(line);
        std::string        cell;
        int index = 0;
        int query, subtopic;
        while(std::getline(lineStream,cell,','))
        {
            if(checkForDouble(cell)){
                if(index == 0){
                    std::stringstream  lineStream(cell);
                    std::string value;
                    std::getline(lineStream,value,'.');

                    query = atoi(value.c_str());
                    std::getline(lineStream,value,'.');
                    subtopic = atoi(value.c_str());
                }
                else{
                    userProfiles_iter = userProfiles[index].find(query);
                    if(userProfiles_iter == userProfiles[index].end()){
                        Profiles *profile = new Profiles();
                        profile->query = query;
                        profile->subtopic_importance.insert(make_pair(subtopic, atof(cell.c_str())));
                        userProfiles[index].insert(make_pair(query, profile));

                    }else
                        userProfiles_iter->second->subtopic_importance.insert(make_pair(subtopic, atof(cell.c_str())));
                }
                index++;
            }
        }
    }
    return userProfiles;
}

#endif	/* UTILS_HPP */

