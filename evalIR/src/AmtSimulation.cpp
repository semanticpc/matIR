/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */

#include "AmtSimulation.hpp"
#include <algorithm>    // std::random_shuffle


AMTSimulation::AMTSimulation(Qrels& qrels)
    :
    _qrels(qrels),
    _qrels_vector(vector<Qrels>()){
    _error_pairs = 0;
    _total_pairs = 0;
    srand ( time(NULL) );
    simulateScores();

}

AMTSimulation::AMTSimulation(Qrels& qrels, vector<Qrels> qrels_vector)
    :
    _qrels(qrels),
    _qrels_vector(qrels_vector){
    _error_pairs = 0;
    _total_pairs = 0;
    srand ( time(NULL) );
    simulateScores();

}

int AMTSimulation::get_preference(arma::rowvec vectorA, arma::rowvec vectorB, arma::rowvec seen){
    double alpha = 0.5;
    double scoreA = 0, scoreB = 0;
    for(int i=0; i < vectorA.n_cols; i++){
        if(vectorA(i) == 1)
            scoreA += pow((1- alpha), seen(i));
    }
    for(int i=0; i < vectorB.n_cols; i++){
        if(vectorB(i) == 1)
            scoreB +=  pow((1 - alpha), seen(i));
    }
    if (scoreA > scoreB)
        return 1;
    else if (scoreA < scoreB)
        return 0;
    else{
        // Resolve ties randomly
        int number = rand() % 10;
        if( number < 5)
            return 1;
        else
            return 0;
    }



}
map<Triplet, int, triplet_comparison> AMTSimulation::sampleTriplets(vector<int> rand_numbers){
    map<Triplet, int, triplet_comparison> triplets;


    int maxIndex = allDocs.size();
    int _random_number = -1;

    for(int i=0; i < 1000; i++){
        std::random_shuffle ( rand_numbers.begin(), rand_numbers.end() );
        Triplet t;
        t.leftDoc = rand_numbers.at(0);
        t.rightDoc = rand_numbers.at(1);;
        t.topDoc = rand_numbers.at(2);;


        ostringstream convert;
        convert << t.topDoc << ":" << t.leftDoc << ":" << t.rightDoc;
        t.code = convert.str();
        if(triplets.find(t) == triplets.end()){
            ostringstream convert;
            convert << t.topDoc << ":" << t.rightDoc << ":" << t.leftDoc;
            t.code = convert.str();
            if(triplets.find(t) == triplets.end())
                triplets.insert(make_pair(t,1));
        }
    }
    return triplets;

}



void AMTSimulation::updateScores(int topDoc, int leftDoc, int rightDoc, bool preferLeft){
    string code; ostringstream convert;
    if(topDoc != -1){
        convert << topDoc;
        code = convert.str();
    }else
        code = "";



    if(preferLeft){
        _utilScores[code](leftDoc) += 1;
        _appearanceCounts[code](leftDoc) += 1;
    }
    else{
        _utilScores[code](rightDoc) += 1;
        _appearanceCounts[code](rightDoc) += 1;
    }


}

bool AMTSimulation::checkPairs(map<Triplet, int, triplet_comparison> &seen_pairs,
                        int left, int right){
    Triplet t;
    ostringstream convert;
    convert << right << ":" << left;
    t.code = convert.str();

    if(seen_pairs.find(t) == seen_pairs.end()){
        ostringstream convert;
        convert << left << ":" << right;
        t.code = convert.str();
        if(seen_pairs.find(t) == seen_pairs.end()){

            seen_pairs.insert(make_pair(t,1));
            return false;
        }

    }
    return true;
}



void AMTSimulation::simulateScores(){


    // Gather all documents into a set
    int index = 0;
    vector<int> rand_numbers;
    map<Document, int>::iterator it = _qrels.relDocs.begin();
    for(;it != _qrels.relDocs.end(); it++){
        allDocs.push_back((*it).first);
        rand_numbers.push_back(index);
        doc_index_map.insert(make_pair((*it).first, index));

        // Initialize
        string code; ostringstream convert;
        convert << index;
        code = convert.str();
        _utilScores[code] = arma::zeros(_qrels.relDocs.size());
        _appearanceCounts[code] = arma::zeros(_qrels.relDocs.size());

        index++;
    }
    set<Document>::iterator itnonrel = _qrels.nonRelDocs.begin();
    for(;itnonrel != _qrels.nonRelDocs.end(); itnonrel++){
        //allDocs.push_back((*itnonrel));
        rand_numbers.push_back(index);
        //doc_index_map.insert(make_pair((*itnonrel), index));

        // Initialize
        string code; ostringstream convert;
        convert << index;
        code = convert.str();
        _utilScores[code] = arma::zeros(_qrels.relDocs.size() );
        _appearanceCounts[code] = arma::zeros(_qrels.relDocs.size());

        index++;
    }




    _utilScores[""] = arma::zeros(_qrels.relDocs.size());
    _appearanceCounts[""] = arma::zeros(_qrels.relDocs.size());

    if( _qrels.relDocs.size() < 3){
        int numOfSubtopics = int(*_qrels.subtopics.rbegin());
        arma::rowvec seen = arma::zeros<arma::rowvec>(numOfSubtopics);
        if(get_preference(_qrels.matrix.row(rand_numbers.at(0)), _qrels.matrix.row(rand_numbers.at(1)), seen)){
            _utilScores[""](0) += 1;
        }else
            _utilScores[""](1) += 1;
        _appearanceCounts[""](0) += 1;
        _appearanceCounts[""](1) += 1;
        _utilScores[""] /= _appearanceCounts[""];
        return;
    }

    // Simulate the preferences for only these 1000 triplets randomly

    vector<Qrels>::iterator qrels_iterator = _qrels_vector.begin();
    for(; qrels_iterator != _qrels_vector.end(); qrels_iterator++){
        Qrels currentQrels = *qrels_iterator;

        // Generate triplets and sample 1000 triplets randomly
        map<Triplet, int, triplet_comparison> triplets = sampleTriplets(rand_numbers);
        map<Triplet, int, triplet_comparison>::iterator triplet_iterator = triplets.begin();
        map<Triplet, int, triplet_comparison> seen_pairs;

        set<int> relevantDocIndices;
        for(;triplet_iterator != triplets.end(); triplet_iterator++){

            Triplet triplet = (*triplet_iterator).first;

            int numOfSubtopics = int(*currentQrels.subtopics.rbegin());
            arma::rowvec seen = arma::zeros<arma::rowvec>(numOfSubtopics);



            if(arma::sum(currentQrels.matrix.row(triplet.topDoc)) <= 0 || triplet.topDoc >=  currentQrels.matrix.n_rows )
                triplet.topDoc = -1;
            if(arma::sum(currentQrels.matrix.row(triplet.leftDoc)) <= 0 || triplet.leftDoc >=  currentQrels.matrix.n_rows)
                    triplet.leftDoc = -1;
            if(arma::sum(currentQrels.matrix.row(triplet.rightDoc)) <= 0 || triplet.rightDoc >=  currentQrels.matrix.n_rows)
                    triplet.rightDoc = -1;


            // Update the scores



            // Keep Track of level 1 pairs for this user profile
            if(triplet.topDoc == -1){
                if(triplet.leftDoc != -1 && triplet.rightDoc != -1){
                    if(checkPairs(seen_pairs, triplet.leftDoc, triplet.rightDoc))
                            continue;
                }else
                    continue;
            }else if(triplet.topDoc != -1){
                seen += currentQrels.matrix.row(triplet.topDoc);
            }



            if(triplet.leftDoc != -1 && triplet.rightDoc != -1){
                if(get_preference(currentQrels.matrix.row(triplet.leftDoc), currentQrels.matrix.row(triplet.rightDoc), seen))
                    updateScores(triplet.topDoc, triplet.leftDoc, triplet.rightDoc, true);
                else
                    updateScores(triplet.topDoc, triplet.leftDoc, triplet.rightDoc, false);
            }
            if(relevantDocIndices.find(triplet.topDoc) == relevantDocIndices.end() && triplet.topDoc != -1)
                relevantDocIndices.insert(triplet.topDoc);
            if(relevantDocIndices.find(triplet.leftDoc) == relevantDocIndices.end() && triplet.leftDoc != -1)
                relevantDocIndices.insert(triplet.leftDoc);
            if(relevantDocIndices.find(triplet.rightDoc) == relevantDocIndices.end() && triplet.rightDoc != -1)
                relevantDocIndices.insert(triplet.rightDoc);

        }
        map<string, arma::vec>::iterator codes = _utilScores.begin();
        for(;codes!= _utilScores.end(); codes++){
            set<int>::iterator relevantDocIterator = relevantDocIndices.begin();

            for(; relevantDocIterator != relevantDocIndices.end(); relevantDocIterator++){
                if((*codes).first != ""){
                    if((*relevantDocIterator) != atoi((*codes).first.c_str()))
                        _utilScores[(*codes).first](*relevantDocIterator) += 1;
                }else
                    _utilScores[(*codes).first](*relevantDocIterator) += 1;
            }
            _appearanceCounts[(*codes).first] += 2;
        }

    }


    map<string, arma::vec>::iterator codes = _utilScores.begin();
    for(;codes!= _utilScores.end(); codes++){

        _utilScores[(*codes).first] /= (_appearanceCounts[(*codes).first] + 2);

    }
}



double AMTSimulation::get_UtilityScore(int docIndex){
    if(docIndex == -1 || docIndex > _utilScores.find("")->second.n_rows)
        return 0;
    else
        return _utilScores.find("")->second(docIndex);
}


pair<int, double> AMTSimulation::get_BestUtilityDoc(){
    arma::uword  index;
    if(_utilScores.find("") == _utilScores.end())
        return make_pair(-1, 0);
    else{
        _utilScores.find("")->second.max(index);
        return make_pair(index, arma::max(_utilScores.find("")->second));
    }
}


double AMTSimulation::get_UtilityScore(int docIndex, int prevDocIndex){
    // Generate a code for the current rankedDocs
    string code; ostringstream convert;
    convert << prevDocIndex;
    code = convert.str();
    if(docIndex == -1)
        return 0;
    else if (docIndex != -1 && prevDocIndex == -1){
        if(_utilScores.find("") != _utilScores.end())
            return _utilScores.find("")->second(docIndex);
        else
            return 0;
    }
    else if (docIndex != -1 && prevDocIndex != -1){
        if(_utilScores.find(code) != _utilScores.end())
            return _utilScores.find(code)->second(docIndex);
        else
            return _utilScores.find("")->second(docIndex);
    }
}

void AMTSimulation::printCounts(){
    cout << "\n Total Pairs : " << _total_pairs << " \nError Pairs : " << _error_pairs << endl;
}

int AMTSimulation::getTotalPairs(){
    return _total_pairs;
}

int AMTSimulation::getTotalRelDocs(){
    return _qrels.matrix.n_rows;
}