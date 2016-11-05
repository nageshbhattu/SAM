/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lara;

/**
 *
 * @author agnibesh
 */
public class WordCode{
    int wordIndex;
    double[] topicCode;
    public WordCode(int wIndex,double[] tc){
        topicCode = tc;
        wordIndex = wIndex;
    }
    public WordCode(int wIndex,int numTopics){
        wordIndex = wIndex;
        topicCode = new double[numTopics];
    }
    public double[] getTopicCode(){
        return topicCode;
    }
    public int getWordIndex(){
        return wordIndex;
    }
}