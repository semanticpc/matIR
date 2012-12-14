/*
 * File:   utils.hpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:50 PM
 */

#ifndef DivQrels_HPP
#define	DivQrels_HPP


#include <armadillo>
#include <map>

using namespace std;

class DivQrels{
public:
    string query;
    arma::mat subtopicMatrix;
    map<string, int> docid_index_map;
    int numOfSubtopics;
    int numOfRelDocuments;

private:
    vector< pair<string,int> > _tmpVector;


public:
    DivQrels();
    DivQrels(string q);
    void generateMatrix();
    void addDocument(string docid, int subtopic, int rel);
    arma::mat getSubtopicMatrix(vector<string> run);


    static map<string, DivQrels*> readDiversityQrelsFile(string qrelsFileName){
        map<string, DivQrels*> qrels;
        std::ifstream qrelsFile(qrelsFileName.c_str(), ios_base::in);
        std::string docid;
        std::string query;
        int subtopic;
        int rel;


        string curQuery = "";
        bool doneParsing = false;


        qrelsFile >> query >> subtopic >> docid >> rel;
        curQuery = query;

        DivQrels* q = new DivQrels(query);
        q->addDocument( docid, subtopic, rel);
        while (true) {

            while(curQuery == query){
                qrelsFile >> query >> subtopic >> docid >> rel;
                q->addDocument( docid, subtopic, rel);
                if( qrelsFile.eof() ) {
                    doneParsing = true;
                    break;
                }

            }

            q->generateMatrix();
            qrels.insert(make_pair(curQuery, q));


            if(doneParsing) break;

            qrelsFile >> query >> subtopic >> docid >> rel;
            curQuery = query;

            DivQrels* q = new DivQrels(query);
            q->addDocument( docid, subtopic, rel);

            }
        return qrels;
    }


};


#endif	/* DivQrels_HPP */

