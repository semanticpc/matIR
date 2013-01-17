/*
 * File:   AdhocMeasures.hpp
 * Author: semanticpc
 *
 * Created on January 16, 2013, 12:51 PM
 */

#ifndef ADHOCMEASURES_HPP
#define	ADHOCMEASURES_HPP


// ERR Computation
static arma::vec computeERR(arma::vec grades, int rank, double alpha){
 arma::vec err_vector = arma::zeros(rank);
 double score = 0.0;
 double decay = 1.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {

    double grade = (pow(2.0, grades(i)) - 1) / pow(2.0, 1);
    score += ( grade * decay) / (i + 1);
    decay *= (1 - grade);
    err_vector(score) = score;

 }
 return err_vector;
}


static double computeDCG(arma::vec grades, int rank, double alpha){

 double score = 0.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {
    if(grades(i) == 1)
        r = 1;
    else
        r = 0;
    score += r / log(i + 2);
 }
 return score;
}

static double computePrecision(arma::vec grades, int rank, double alpha){

 double score = 0.0;
 double r = 0.0;
 for (int i = 0; i < rank; i++) {
    if(grades(i) == 1)
        r = 1;
    else
        r = 0;
    score += r;
 }
 return score/rank;
}


#endif	/* ADHOCMEASURES_HPP */

