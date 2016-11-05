/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lara;

import java.util.ArrayList;

/**
 *
 * @author agnibesh
 */
public class DocWordCode {
    ArrayList<WordCode> wordCodes;
    int[] indices;
    int[] counts;
    public DocWordCode(int numWordCodes){
        wordCodes = new ArrayList<WordCode>();
    }
    public void addWordCode(int wordInd,double[] wCode){
        wordCodes.add(new WordCode(wordInd,wCode));
    }
    public ArrayList<WordCode> getWordCodes(){
        return wordCodes;
    }
    public WordCode getWordCode(int ind){
        return wordCodes.get(ind);
    }
}
