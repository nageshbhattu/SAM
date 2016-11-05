/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aspectSegmenter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Vector;

/**
 *
 * @author agnibesh
 */
public class Reviews2Vectors {
    public static ArrayList<Document> dlist = new ArrayList<>();
    public static double[][] aspectRatingArray;
    public static int vocabularySize;
    public static int numAspects=Analyzer.ASPECT_CONTENT_CUT;
    public static int numUsers;
    public static int numItems;
   /* public Reviews2Vectors(int na){
        numAspects = na;
    }*/
    public  void createReviewVectors(Vector<Hotel> hotelList, Hashtable<String,Integer> vocabulary){
        int  wordID, outputSize=0, reviewSize=0; 
        String reviewID;
        int itemID = 0;
        for(Hotel hotel:hotelList){
            for(Review r:hotel.m_reviews){//aggregate all the reviews
                //collect the vectors
                HashMap<Integer,Integer> wordCounts = new HashMap<Integer,Integer>();
                for(Review.Sentence stn:r.m_stns){

                    for(Review.Token t:stn.m_tokens){//select the in-vocabulary word
                        if (vocabulary.containsKey(t.m_lemma)){
                            wordID = vocabulary.get(t.m_lemma);
                            wordCounts.put(wordID, wordCounts.containsKey(wordID)? wordCounts.get(wordID)+1:1);
                        }
                    }
                }
                int[] wordIds = wordCounts.keySet().stream().mapToInt(i -> i).toArray();
                int[] counts = wordCounts.values().stream().mapToInt(i -> i).toArray();
                Document doc = new Document(wordIds,counts,intToDouble(r.m_ratings),r.m_userID, itemID);
                dlist.add(doc);
              // System.out.println(doc);
              
              System.out.println("abc  "+r.toString());
            }
            itemID++;
        }
        vocabularySize = vocabulary.size();
        
    }
    public static double[] intToDouble(int[] ratings){
        double[] doubleRatings = new double[ratings.length];
        for(int i  =0;i<ratings.length;i++){
            doubleRatings[i] = ratings[i];
        }
        return doubleRatings;
    }
    public static void saveReviewVectors(String filename) throws IOException{
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
        bw.write(numAspects+ "\n");
        for(int di = 0;di<dlist.size();di++){
            Document doc = dlist.get(di);
            bw.write(doc.toString());
        }
    }
    public static void readReviewVectors(String filename) throws FileNotFoundException, IOException{
        BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
        String line;
        line = br.readLine();
        numAspects = Integer.parseInt(line);
        int maxUserID = 0; int maxItemID = 0;
        int maxWordID = 0;
        while((line=br.readLine())!=null){
            
            String[] pairs = line.split(" ");
            int[] wordIds = new int[pairs.length-(numAspects+3)];
            int[] wCounts = new int[pairs.length-(numAspects+3)];
            double[] ratings = new double[numAspects+1];
            
            for(int ri = 0;ri<=numAspects; ri++){
                pairs[ri] = pairs[ri].replace(",", "");
                pairs[ri] = pairs[ri].replace("]", "");
                pairs[ri] = pairs[ri].replace("[", "");
                ratings[ri]= Double.parseDouble(pairs[ri]);
            }
            int userID = Integer.parseInt(pairs[numAspects+1]);
            int itemID = Integer.parseInt(pairs[numAspects+2]);
            
            if(userID>maxUserID)
                maxUserID = userID;
            
            if(itemID>maxItemID)
                maxItemID = itemID;
            
            for(int pi = ratings.length + 2;pi<pairs.length;pi++){
                int wi = pi-numAspects-3;
                if(pairs[pi].contains(":")){
                    String[] wcs = pairs[pi].split(":");
                    wordIds[wi] = Integer.parseInt(wcs[0]);
                    wCounts[wi] = Integer.parseInt(wcs[1]);
                    if(maxWordID<wordIds[wi])
                        maxWordID = wordIds[wi];
                }
            }
            Document d = new Document(wordIds,wCounts,ratings,userID,itemID);
            
            dlist.add(d);
            
        }
        numUsers = maxUserID+1;
        numItems = maxItemID+1;
        vocabularySize = maxWordID+1;
    }
    public double[][] getAspectArray(){
        aspectRatingArray = new double[dlist.size()][numAspects];
        for(int di =0;di<dlist.size();di++){
        Document doc = dlist.get(di);
            for(int ai=0;ai<numAspects;ai++){
            double val = doc.ratings[ai];
            aspectRatingArray[di][ai]=val;
            }
        }
        return aspectRatingArray;
        // System.out.println(Arrays.toString(aspectRatingArray));       
    }
    public static void main(String[] args) throws IOException{
        Analyzer analyzer = new Analyzer("Data/Seeds/hotel_bootstrapping.dat", "Data/Seeds/stopwords.dat", 
				"Data/Model/NLP/en-sent.zip", "Data/Model/NLP/en-token.zip", "Data/Model/NLP/en-pos-maxent.bin");
		//analyzer.LoadVocabulary("Data/Seeds/hotel_vocabulary_CHI.dat");
        analyzer.LoadDirectory("Data/Reviews/", ".dat");
        Reviews2Vectors r2v = new Reviews2Vectors();
        
        r2v.createReviewVectors(analyzer.m_hotelList,analyzer.m_vocabulary);
        saveReviewVectors("Data/hotel_reviews.dat");
       // getAspectArray();
    }
}
