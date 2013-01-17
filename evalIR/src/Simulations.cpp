/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */

#include "Simulations.hpp"
PrefSimulation::PrefSimulation(Qrels& qrels, int error_rate, int missing_rate)
    :
    _error_rate(error_rate),
    _missing_rate(missing_rate),
    _qrels(qrels),
    _qrels_vector(vector<Qrels>()){
    _utilScores.insert(make_pair("",simulate_level(vector<int>())));

}

PrefSimulation::PrefSimulation(Qrels& qrels, vector<Qrels> qrels_vector, int error_rate, int missing_rate)
    :
    _error_rate(error_rate),
    _missing_rate(missing_rate),
    _qrels(qrels),
    _qrels_vector(qrels_vector){
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

        //int number = rand() % 10;
        int number = 5;
        if( number < 5)
            return 1;
        else
            return 1;
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
            lvl_score += utils.first;
            appearance_count += utils.second;
        }
    }else{
        utils =get_simulation_scores(_qrels, rankedDocs);
        lvl_score = utils.first;
        appearance_count = utils.second;
    }

    // Plus One Smoothing to avoid zero utility  scores
    lvl_score += 1;
    lvl_score /= (appearance_count + 1);

    // All documents in the rank list must get a score of zero
    for(int i=0;i<rankedDocs.size();i++)
        lvl_score(rankedDocs.at(i)) = 0.0;


    return lvl_score;
}
pair<arma::vec,arma::vec> PrefSimulation::get_simulation_scores(Qrels &qrels, vector<int> rankedDocs){
    srand ( time(NULL) );
    int numOfSubtopics = int(*qrels.subtopics.rbegin());
    // Initialize Seen Subtopic Counts
    arma::rowvec seen = arma::zeros<arma::rowvec>(numOfSubtopics);
    for(vector<int>::iterator i= rankedDocs.begin(); i != rankedDocs.end(); i++){
        if((*i) != -1)
            seen += qrels.matrix.row((*i));
    }

    arma::vec lvl_score = arma::zeros(qrels.matrix.n_rows);
    arma::vec appearance_count = arma::zeros(qrels.matrix.n_rows);

    for(int i=0; i < qrels.matrix.n_rows; i++){
        // Ignore documents already present in the rank list
        if(find(rankedDocs.begin(), rankedDocs.end(), i) != rankedDocs.end())
            continue;


        double score = 0.0;
        int appearance = 0;

        ;
        // Generate pairs
        for(int j= i; j < qrels.matrix.n_rows; j++){
            // Ignore documents already present in the rank list
            if( i == j || (find(rankedDocs.begin(), rankedDocs.end(), j) != rankedDocs.end()) )
                continue;

            if(arma::sum(_qrels.matrix.row(i)) <= 0)
                continue;

            // Obtain the preference for each pair
            //srand ( time(NULL) );
            //int missing_random_number = (rand() * 1.0 / RAND_MAX) * 100;
            int missing_random_number = 10;
            if(missing_random_number > _missing_rate){
            // If there are more than one user profile consider them as well
                //srand ( time(NULL) );
                //int eror_random_number = (rand() * 1.0 / RAND_MAX) * 100;
                int eror_random_number= 10;
                if(eror_random_number > _error_rate){
                    if(get_preference(_qrels.matrix.row(i), _qrels.matrix.row(j), seen) == 1)
                        lvl_score(i) += 1;
                    else
                        lvl_score(j) += 1;
                }else{
                    // We are flipping the preference here
                    if(get_preference(_qrels.matrix.row(i), _qrels.matrix.row(j), seen) == 1)
                        lvl_score(j) += 1;
                    else
                        lvl_score(i) += 1;

                }
            }


            appearance_count(i) += 1;
            appearance_count(j) += 1;
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