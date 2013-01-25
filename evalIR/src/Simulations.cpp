/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */

#include "Simulations.hpp"
#include <algorithm>    // std::random_shuffle
PrefSimulation::PrefSimulation(Qrels& qrels, double error_rate, double missing_rate)
    :
    _error_rate(error_rate),
    _missing_rate(missing_rate),
    _qrels(qrels),
    _qrels_vector(vector<Qrels>()){
    _error_pairs = 0;
    _total_pairs = 0;
    srand ( time(NULL) );
    _utilScores.insert(make_pair("",simulate_level(vector<int>())));

}

PrefSimulation::PrefSimulation(Qrels& qrels, vector<Qrels> qrels_vector, double error_rate, double missing_rate)
    :
    _error_rate(error_rate),
    _missing_rate(missing_rate),
    _qrels(qrels),
    _qrels_vector(qrels_vector){
    _error_pairs = 0;
    _total_pairs = 0;
    srand ( time(NULL) );
    _utilScores.insert(make_pair("",simulate_level(vector<int>())));
}

int PrefSimulation::get_preference(arma::rowvec vectorA, arma::rowvec vectorB, arma::rowvec seen){
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
arma::vec PrefSimulation::simulate_level(vector<int> rankedDocs=vector<int>()){
    pair<arma::vec,arma::vec> utils;
    arma::vec lvl_score = arma::zeros(_qrels.numOfRelDocs);
    arma::vec appearance_count = arma::zeros(_qrels.numOfRelDocs);;
    if(_qrels_vector.size() > 0){
        vector<Qrels>::iterator it = _qrels_vector.begin();
        for(; it != _qrels_vector.end(); it++){
            utils =get_simulation_scores(*it, rankedDocs);
            for(int i=0; i < (*it).matrix.n_rows; i++){
                if(arma::sum((*it).matrix.row(i)) > 0){
                    utils.first(i) += 1;
                    utils.second(i) += 2;
                }
            }
            lvl_score += utils.first;
            appearance_count += utils.second;
        }
        lvl_score /= appearance_count;
    }else{
        utils =get_simulation_scores(_qrels, rankedDocs);
        lvl_score = utils.first;
        appearance_count = utils.second;
        // Add one only to relevant documents
        for(int i=0; i < _qrels.matrix.n_rows; i++){
            if(appearance_count(i) > 0){
                lvl_score(i) += 1;
            }
        }
        lvl_score /= (appearance_count+2);

    }

    // All documents in the rank list must get a score of zero
    for(int i=0;i<rankedDocs.size();i++)
        lvl_score(rankedDocs.at(i)) = 0.0;

    return lvl_score;
}
pair<arma::vec,arma::vec> PrefSimulation::get_simulation_scores(Qrels &qrels, vector<int> rankedDocs){

    int numOfSubtopics = int(*qrels.subtopics.rbegin());
    // Initialize Seen Subtopic Counts
    arma::rowvec seen = arma::zeros<arma::rowvec>(numOfSubtopics);
    for(vector<int>::iterator i= rankedDocs.begin(); i != rankedDocs.end(); i++){
        if((*i) != -1)
            seen += qrels.matrix.row((*i));
    }

    arma::vec lvl_score = arma::zeros(qrels.matrix.n_rows);
    arma::vec appearance_count = arma::zeros(qrels.matrix.n_rows);

    for(int i= 0; i < qrels.matrix.n_rows; i++){
        // Ignore documents already present in the rank list
        if(find(rankedDocs.begin(), rankedDocs.end(), i) != rankedDocs.end())
            continue;


        double score = 0.0;
        int appearance = 0;

        if(arma::sum(qrels.matrix.row(i)) <= 0)
            continue;
        // Generate pairs
        for(int j= i; j < qrels.matrix.n_rows; j++){

            // Ignore documents already present in the rank list
            if( i == j || (find(rankedDocs.begin(), rankedDocs.end(), j) != rankedDocs.end()) )
                continue;

            if(arma::sum(qrels.matrix.row(j)) <= 0)
                continue;

            // Obtain the preference for each pair
            int missing_random_number = ((rand() * 1.0 / RAND_MAX) * 100 + 1);
            if(_missing_rate < missing_random_number){
                int eror_random_number = ((rand() * 1.0 / RAND_MAX) * 100) + 1;
                if(eror_random_number >= _error_rate){
                    if(get_preference(qrels.matrix.row(i), qrels.matrix.row(j), seen))
                        lvl_score(i) += 1;
                    else
                        lvl_score(j) += 1;
                }else {
                    // We are flipping the preference here
                    _error_pairs++;
                    if(get_preference(qrels.matrix.row(i), qrels.matrix.row(j), seen))
                        lvl_score(j) += 1;
                    else
                        lvl_score(i) += 1;

                }
                _total_pairs++;
                appearance_count(i) += 1;
                appearance_count(j) += 1;
            }
        }
    }

    return make_pair(lvl_score,appearance_count) ;
}


double PrefSimulation::get_UtilityScore(int docIndex){
    if(docIndex == -1 || docIndex > _utilScores.find("")->second.n_rows)
        return 0;
    else
        return _utilScores.find("")->second(docIndex);
}


pair<int, double> PrefSimulation::get_BestUtilityDoc(){
    arma::uword  index;
    _utilScores.find("")->second.max(index);
    return make_pair(index, arma::max(_utilScores.find("")->second));
}


double PrefSimulation::get_UtilityScore(int docIndex, int prevDocIndex){
    // Generate a code for the current rankedDocs
    string code; ostringstream convert;
    convert << prevDocIndex;
    code = convert.str();

    if(docIndex == -1)
        return 0;
    else if (docIndex != -1 && prevDocIndex == -1)
        return _utilScores.find("")->second(docIndex);
    else if (docIndex != -1 && prevDocIndex != -1){
        if(_utilScores.find(code) != _utilScores.end())
            return _utilScores.find(code)->second(docIndex);
        else{
            vector<int> rankedDocs;
            rankedDocs.push_back(prevDocIndex);
            arma::vec lvl_score = simulate_level(rankedDocs);
            _utilScores.insert(make_pair(code, lvl_score));
            return lvl_score(docIndex);
        }
    }
}


pair<int, double> PrefSimulation::get_BestUtilityDoc(int prevDocIndex){
    // Generate a code for the current rankedDocs
    string code; ostringstream convert;
    convert << prevDocIndex;
    code = convert.str();

    if(prevDocIndex == -1)
        return get_BestUtilityDoc();
    else if (_utilScores.find(code) == _utilScores.end()){
        arma::uword  index;
        _utilScores.find(code)->second.max(index);
        return make_pair(index, arma::max(_utilScores.find(code)->second));
    }else{
            vector<int> rankedDocs;
            rankedDocs.push_back(prevDocIndex);
            arma::vec lvl_score = simulate_level(rankedDocs);
            _utilScores.insert(make_pair(code, lvl_score));
            arma::uword  index;
            lvl_score.max(index);
            return make_pair(index, arma::max(lvl_score));
    }
}

void PrefSimulation::printCounts(){
    cout << "\n Total Pairs : " << _total_pairs << " \nError Pairs : " << _error_pairs << endl;
}

int PrefSimulation::getTotalPairs(){
    return _total_pairs;
}

int PrefSimulation::getTotalRelDocs(){
    return _qrels.matrix.n_rows;
}