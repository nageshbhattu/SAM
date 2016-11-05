/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SentimentAnalysis;

import aspectSegmenter.Analyzer;
import lara.LRR;

/**
 *
 * @author nageshbhattu
 */
public class SALARATest {
    public static void main(String[] args){
        /*Analyzer analyzer = new Analyzer("Data/Seeds/hotel_bootstrapping.dat",
                                        "Data/Seeds/stopwords.dat", "Data/Model/NLP/en-sent.zip",
                                        "Data/Model/NLP/en-token.zip", "Data/Model/NLP/en-pos-maxent.bin");
        analyzer.LoadDirectory("Data/Reviews/", ".dat");
        analyzer.BootStrapping("Data/Seeds/hotel_bootstrapping_sel_new.dat");
    
        analyzer.Save2Vectors("Data/Vectors/vector_CHI_4000.dat");*/
        LRR model = new LRR(500, 1e-2, 5000, 1e-2, 2.0);
        model.EM_est("Data/Vectors/vector_CHI_4000.dat", 10, 1e-4);
        model.SavePrediction("Data/Results/prediction.dat");
        model.SaveModel("Data/Model/model_hotel.dat");
    }
}
