/*
 * File:   main.cpp
 * Author: semanticpc
 *
 * Created on December 4, 2012, 6:45 PM
 */
#include "DivMeasures.hpp"
using namespace std;

static void printResults(map<int, vector<Document> > run, map<int, Qrels> qrels){


    double numOfQ = qrels.size();
    map<int, Qrels>::iterator it;

    double srecallSum5, srecallSum10, srecallSum20;
    cout << "topic,srecall@5,srecall@10,srecall@20";

    double andcgSum5, andcgSum10, andcgSum20;
    cout << "andcg@5,andcg@10,andcg@20";
    double erriaSum5, erriaSum10, erriaSum20;
    cout << "erria@5,erria@10,erria@20";


    arma::vec prfScoreSum5 = arma::zeros(8);
    arma::vec prfScoreSum10 = arma::zeros(8);
    arma::vec prfScoreSum20 = arma::zeros(8);
    cout << "prf_ave_none@5,prf_ave_none@10,prf_ave_none@20";
    cout << "prf_ave_RBP@5,prf_ave_RBP@10,prf_ave_RBP@20";
    cout << "prf_ave_RR@5,prf_ave_RR@10,prf_ave_RR@20";
    cout << "prf_ave_DCG@5,prf_ave_DCG@10,prf_ave_DCG@20";

    cout << "prf_min_none@5,prf_min_none@10,prf_min_none@20";
    cout << "prf_min_RBP@5,prf_min_RBP@10,prf_min_RBP@20";
    cout << "prf_min_RR@5,prf_min_RR@10,prf_min_RR@20";
    cout << "prf_min_DCG@5,prf_min_DCG@10,prf_min_DCG@20";

    cout << endl;
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        int query = it->first;
        Qrels qrels = it->second;
        int rank = 20;
        arma::mat run_matrix = judge_diversity(run.find(query)->second, qrels, rank);
        cout << query;
/**/
        // Get Subtopic-Recall Score
        arma::vec srecallScore = s_recall(run_matrix, qrels, rank);
        cout << "," << srecallScore(4);
        cout << "," << srecallScore(9);
        cout << "," << srecallScore(19);
        srecallSum5 += srecallScore(4);
        srecallSum10 += srecallScore(9);
        srecallSum20 += srecallScore(19);

        // Get Alpha-nDCG Score
        arma::vec andcgScore = andcg(run_matrix, qrels, rank);
        cout << "," << andcgScore(4);
        cout << "," << andcgScore(9);
        cout << "," << andcgScore(19);
        andcgSum5 += andcgScore(4);
        andcgSum10 += andcgScore(9);
        andcgSum20 += andcgScore(19);

        // Get ERR-IA Score
        arma::vec erriaScore = erria(run_matrix, qrels, rank);
        cout << "," << erriaScore(4);
        cout << "," << erriaScore(9);
        cout << "," << erriaScore(19);
        erriaSum5 += erriaScore(4);
        erriaSum10 += erriaScore(9);
        erriaSum20 += erriaScore(19);

        map<string, arma::vec> prfScore = pref_measure(run.find(query)->second, qrels, rank);
        map<string, arma::vec>::iterator prefScore_iter;
        int i =0;
        for(prefScore_iter = prfScore.begin();prefScore_iter != prfScore.end(); prefScore_iter++ ){
                cout << "," << prefScore_iter->second(4);
                cout << "," << prefScore_iter->second(9);
                cout << "," << prefScore_iter->second(19);
                prfScoreSum5(i) += prefScore_iter->second(4);
                prfScoreSum10(i) += prefScore_iter->second(9);
                prfScoreSum20(i) += prefScore_iter->second(19);
                i++;
        }


        cout << endl;

    }
    cout << "all," << srecallSum5/numOfQ << "," << srecallSum10/numOfQ << "," << srecallSum20/numOfQ
            << "," << andcgSum5/numOfQ <<  "," << andcgSum10/numOfQ << "," << andcgSum20/numOfQ
            << "," << erriaSum5/numOfQ << "," << erriaSum10/numOfQ   << "," << erriaSum20/numOfQ;
    map<string, arma::vec>::iterator prefScore_iter;
    for(int i =0; i< 8; i ++ ){
        cout << "," << prfScoreSum5(i)/numOfQ;
        cout << "," << prfScoreSum10(i)/numOfQ;
        cout << "," << prfScoreSum20(i)/numOfQ;
    }
    cout << endl;
}




static void printResults(map<int, vector<Document> > &run, map<int, Qrels> &qrels, map<int, map<int, Profiles* > > &profiles){
    double srecallSum, andcgSum, erriaSum;
    arma::vec prfScoreSum = arma::zeros(8);
    double numOfQ = qrels.size();
    map<int, Qrels>::iterator it;


     double srecallSum5, srecallSum10, srecallSum20;
    cout << "topic,srecall@20";


    cout << ",andcg@20";

    cout << ",erria@20";

    cout << ",prf_ave_none@20";
    cout << ",prf_ave_RBP@20";
    cout << ",prf_ave_RR@20";
    cout << ",prf_ave_DCG@20";

    cout << ",prf_min_none@20";
    cout << ",prf_min_RBP@20";
    cout << ",prf_min_RR@20";
    cout << ",prf_min_DCG@20";
    cout << endl;
    // Iterate through each query
    for ( it=qrels.begin() ; it != qrels.end(); it++ ){
        int query = it->first;
        Qrels qrels = it->second;
        int rank = 20;
        //if(query != 52)
          //  continue;
        cout << query;


        // Iterate user profiles
        map<int, map<int, Profiles* > >::iterator iter;
        double srecallSum_query = 0, andcgSum_query = 0, erriaSum_query = 0;
        int numOfProfiles = 0;
        vector<Qrels> new_qrels_vector;
        for(iter=profiles.begin();iter!=profiles.end();iter++){

            Qrels new_qrels;
            map<int, Profiles* >::iterator find_query = iter->second.find(query);
            if(find_query != iter->second.end())
                new_qrels = update_SubtopicImportance(qrels, find_query->second->subtopic_importance);
            else
                continue;

            if(arma::sum(new_qrels.subtopicImportance) <= 0)
                continue;
            new_qrels_vector.push_back(new_qrels);
            arma::mat run_matrix = judge_diversity(run.find(query)->second, new_qrels, rank);


            // Get Subtopic-Recall Score
            arma::vec srecallScore = s_recall(run_matrix, new_qrels, rank);
            srecallSum_query += srecallScore(rank - 1);

            // Get Alpha-nDCG Score
            arma::vec andcgScore = andcg(run_matrix, new_qrels, rank);
            andcgSum_query += andcgScore(rank - 1);

            // Get ERR-IA Score
            arma::vec erriaScore = erria(run_matrix, new_qrels, rank);
            erriaSum_query += erriaScore(rank - 1);

            numOfProfiles++;

        }

        srecallSum += srecallSum_query/numOfProfiles;
        andcgSum += andcgSum_query/numOfProfiles;
        erriaSum += erriaSum_query/numOfProfiles;
        cout << "," << srecallSum_query/numOfProfiles
             << "," << andcgSum_query/numOfProfiles
             << "," << erriaSum_query/numOfProfiles;

        map<string, arma::vec> prfScore = pref_measure(run.find(query)->second, rank, new_qrels_vector, qrels);
        map<string, arma::vec>::iterator prefScore_iter;
        int i =0;
        for(prefScore_iter = prfScore.begin();prefScore_iter != prfScore.end(); prefScore_iter++ ){
                cout << "," << prefScore_iter->second(rank - 1);
                prfScoreSum(i++) += prefScore_iter->second(rank - 1);
        }

        cout << endl;
        //break;

    }
    cout << "all," << srecallSum/numOfQ << "," << andcgSum/numOfQ << "," << erriaSum/numOfQ;
    map<string, arma::vec>::iterator prefScore_iter;
    for(int i =0; i< 8; i ++ )
            cout << "," << prfScoreSum(i)/numOfQ;
    cout << endl;
}

int main(int argc, char** argv) {
    string usage = "Usage:\n\t divEval <qrels_file> <run_file>\n";
    if (argc < 3){
        cout << usage;
        exit(0);
    }
    string qrelsFile = argv[1], runFile = argv[2];
    string profilesFile = "";
    double error, missing;
    int nextIndex = 3;

    for(int nextIndex= 3; nextIndex < argc; nextIndex++){
        if(argv[nextIndex] == "-e"){
            error = atof(argv[nextIndex + 1]);
            nextIndex += 1;
        }
        else if(argv[nextIndex] == "-m" ){
            missing = atof(argv[nextIndex + 1]);
            nextIndex += 2;
        }
        else{
            profilesFile = argv[nextIndex];
            nextIndex += 1;
        }
    }

   if (runFile == "" || qrelsFile == ""){
        cout << usage;
        exit(0);
   }




    // Get Qrels File
    map<int, Qrels> qrels = readDiversityQrelsFile(qrelsFile);

    // Read the run file
    map<int, vector<Document> > run = readRunFile(runFile);

    if (profilesFile == "")
        printResults(run, qrels);
    else{
        map<int, map<int, Profiles* > > profiles = read_userProfiles(profilesFile);
        printResults(run, qrels, profiles);
    }


    return 0;
}

