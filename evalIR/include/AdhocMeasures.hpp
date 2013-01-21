/*
 * File:   AdhocMeasures.hpp
 * Author: semanticpc
 *
 * Created on January 16, 2013, 12:51 PM
 */

#ifndef ADHOCMEASURES_HPP
#define	ADHOCMEASURES_HPP


// ERR Computation
static double err(arma::vec grades, int rank){
 //arma::vec err;
 double score = 0.0;
 double decay = 1.0;
 double r = 0.0;
 int g;
 for (int i = 0; i < rank; i++) {
     if(grades(i) > 0)
        g = 3;
     else
        g = 0;
    double grade = (pow(2.0, g) - 1) / pow(2.0, 4);
    score += ( grade * decay) / (i + 1);
    decay *= (1 - grade);
    //err_vector(score) = score;

 }
 return score;
}





static double DCG(arma::vec run, int rank){

 double score = 0.0;
 double r = 0.0;
 for (int i = 0; i <= rank; i++) {
    if(run(i) > 0)
        r = (pow(2.0, 3) - 1);// / pow(2.0, 4);
    else
        r = 0;
    score +=  r/ (log (i + 2)/ log(2));

 }
 return score;
}

static double nDCG(arma::vec run, arma::vec qrels, int rank){
    return DCG(run, rank) / DCG(arma::sort(qrels, 1), rank);
}

static double precision(arma::vec grades, int rank){

 double score = 0.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {
    if(grades(i) > 0)
        r = 1;
    else
        r = 0;
    score += r;
 }
 return score/rank;
}


#endif	/* ADHOCMEASURES_HPP */

