package sacm;


import aspectSegmenter.Analyzer;
import aspectSegmenter.Document;
import aspectSegmenter.Reviews2Vectors;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import lara.DocWordCode;
import lara.WordCode;
import optimizer.LBFGS;

/**
 * Created by AgniBesh on 20/6/16.
*/
public class SACModel  {
    
    public static void main(String[] args) throws IOException{
    SACModel s = new SACModel();
    s.init();
    s.initWordCodes();
    s.sacm();
    s.test();
    
    }
    
    
    
    
    int numAspects = Analyzer.ASPECT_CONTENT_CUT;
    Comparision cpr;
    Reviews2Vectors rv = new Reviews2Vectors();
    double[][] reviewAspectRatings;
    double[] reviewAspectRatingGradient;
    double[] UserCodeGradient;
    double[] AlphaCodeGradient;
    double[] ItemCodeGradient;
    double[] reviewAspectRatingsLR;
    double[] wordCodeGradient;
    double[][] beta;
    double[][] userCodes;
    double[][] itemCodes;
    double[] overallRatings;
    ArrayList<DocWordCode> docWordCodes ;
    ArrayList<Document> docList;
    double var = 0.0;
    double sum;
    int numDocs;
    double ratingVariance = 0.1;
    double alpha = 1; 
    int vocabularySize;
    double gamma  = 0.5;
    double rho = 5e-4;
    double lambda = 0.5;
    private double m_dPoisOffset = 0.0001;
    
    
    
    void test(){
      double[][] x = reviewAspectRatings;
      double[][] y = rv.getAspectArray();
        cpr.mse(x, y);
       for(int di=0;di<numDocs;di++){
        System.out.println(Arrays.toString(x[di]));
        System.out.println(Arrays.toString(y[di]));
        }
    }
   
    void initWordCodes(){
       for(int di = 0;di<docList.size();di++){  
            Document d = docList.get(di);  
            int[] wIndices = d.indices;
            DocWordCode dwc = new DocWordCode(wIndices.length);        
            for(int wi=0;wi<wIndices.length;wi++){
                double[] tc = new double[numAspects]; 
                dwc.addWordCode(wIndices[wi], tc);
                    for(int ai = 0;ai<numAspects;ai++)
                    tc[ai] = 1/vocabularySize;
            }
            docWordCodes.add(dwc);
          
       }
    }
    
    private int numUsers;
    private int numItems;
    public void init() throws IOException{
        Reviews2Vectors.readReviewVectors("Data/hotel_reviews.dat");
        cpr = new Comparision();
        numDocs = Reviews2Vectors.dlist.size();
        numUsers = Reviews2Vectors.numUsers;
        numItems = Reviews2Vectors.numItems;
        vocabularySize = Reviews2Vectors.vocabularySize;
        docList = Reviews2Vectors.dlist;
        docWordCodes = new ArrayList<DocWordCode>();
        AlphaCodeGradient = new double[numDocs];
        reviewAspectRatings = new double[numDocs][numAspects];
        reviewAspectRatingGradient = new double[numAspects];
        wordCodeGradient = new double[numAspects];
        UserCodeGradient = new double[numAspects];
        ItemCodeGradient = new double[numAspects];
        userCodes = new double[numDocs][numAspects];
        itemCodes = new double[numDocs][numAspects];
        overallRatings = new double[numDocs];
        beta = new double[vocabularySize][numAspects];
        betaGradient = new double[numAspects * vocabularySize];
    }
    public void sacm() {
        for(int di = 0;di<numDocs;di++){
            for(int ai= 0;ai<numAspects;ai++){
                reviewAspectRatings[di][ai] = 3.0;
                userCodes[di][ai] = 1;
                itemCodes[di][ai]  = 1;
            }
        }
        double[] thetad = new double[numAspects];
        double[] etad = new double[numAspects];
        for(int iter = 0;iter<1000;iter++)
        for(int di=0;di<numDocs;di++){
            if(di%10000==0)
            System.out.println("DocIndex"+ di);
            Document doc = docList.get(di);
            double sum = 0.0;
            overallRatings[di] = doc.getOverallRating();
            double meanAspectRating = 0.0;
            for(int ai = 0;ai<numAspects;ai++){
                thetad[ai] = userCodes[di][ai]*itemCodes[di][ai];
                etad[ai]  = Math.exp(thetad[ai]);
                sum+=etad[ai];
                
            }
            for(int ai = 0;ai<numAspects;ai++){
                etad[ai]/=sum;
                meanAspectRating +=etad[ai]*reviewAspectRatings[di][ai];
            }
            for(int ai = 0;ai<numAspects;ai++){
                //gradient with -ve sign (descent direction)
                reviewAspectRatingGradient[ai] = -(1/ratingVariance )* (overallRatings[di]-meanAspectRating)*etad[ai]+
                                                (1/(alpha*alpha*userCodes[di][ai]*userCodes[di][ai])) * 
                                                (reviewAspectRatings[di][ai]-itemCodes[di][ai]);
                if(di==0){
                System.out.println("etad: "+etad[ai]);
                System.out.println("usercode: "+userCodes[di][ai]);
                System.out.println("itemcode: "+itemCodes[di][ai]);
                }
            }
            double alphaArmijo = getAlphaAspectRatingsArmijo(overallRatings[di],reviewAspectRatings[di],
                                        reviewAspectRatingGradient,etad,itemCodes[di],userCodes[di]);
            for(int ai = 0;ai<numAspects;ai++){
                reviewAspectRatings[di][ai] += reviewAspectRatingGradient[ai]*alphaArmijo;
//                System.out.println("Armijo: "+alphaArmijo);
  //              System.out.println("Overall: "+overallRatings[di]);
               
            }
            

            ArrayList<WordCode> currentDocWordCodes = docWordCodes.get(di).getWordCodes();
            for(int vi=0;vi<currentDocWordCodes.size();vi++){
                WordCode wc = currentDocWordCodes.get(vi);
                double[] tc = wc.getTopicCode();
                int wIndex = wc.getWordIndex();
                double dotProduct = 0.0;
                for(int ai = 0;ai<numAspects;ai++){
                    dotProduct += beta[wIndex][ai]*tc[ai];
                }
                // Dot Product contains s_dn^t . beta_n
                
                for(int ai = 0;ai<numAspects;ai++){
                    
                    // gradient with -ve sign(descent direction)
                   
                    
                    wordCodeGradient[ai] = -(gamma *2*(tc[ai]-thetad[ai])+ rho + 
                                            (1/dotProduct * beta[wIndex][ai]) + beta[wIndex][ai]);
                    wordCodeGradient[ai] = Math.max(-tc[ai], wordCodeGradient[ai]-rho);
                    
                }
                double alphaWordCodeArmijo = getAlphaWordCodeArmijo(tc,thetad,beta[wIndex],wordCodeGradient, 0);
                for(int ai = 0;ai<numAspects;ai++){
                    tc[ai] += wordCodeGradient[ai]*alphaWordCodeArmijo;
                }
            }
            
        }
        
        //for each user
        // 3: Update the T Matrix using gradient Descent
      /*  for(int u = 0;u<numUsers;u++){
          
            for(int ai = 0;ai<numAspects;ai++){
           var = (1/(alpha*alpha*userCodes[u][ai]*userCodes[u][ai]*userCodes[u][ai])) * 
                    ((reviewAspectRatings[u][ai]-itemCodes[u][ai])*
                    (reviewAspectRatings[u][ai]-itemCodes[u][ai]));
           
           
            }
            for(int ai = 0;ai<numAspects;ai++){
            UserCodeGradient[ai] = -(lambda+var+(1/(alpha*userCodes[u][ai])));
            UserCodeGradient[ai] = Math.max(UserCodeGradient[ai],UserCodeGradient[ai]-rho);
            
            }
            System.out.println(Arrays.toString(reviewAspectRatings[u]));
            double alphaUserCodeArmijo = getAlphaUserCodeArmijo(userCodes[u],reviewAspectRatings[u],itemCodes[u],UserCodeGradient);
            for(int ai = 0;ai<numAspects;ai++){
                userCodes[u][ai] += UserCodeGradient[ai]*alphaUserCodeArmijo;
            }
        }
        */
        
        
        for(int di=0;di<numDocs;di++){
           /* Document doc = docList.get(di);
            int u = doc.userID;
            int h = doc.itemID;*/
            for(int ai = 0;ai<numAspects;ai++){
            var = (1/(alpha*alpha*userCodes[di][ai]*userCodes[di][ai]*userCodes[di][ai])) * 
                    ((reviewAspectRatings[di][ai]-itemCodes[di][ai])*
                    (reviewAspectRatings[di][ai]-itemCodes[di][ai]));
            //System.out.println(userCodes[di][ai]+"  "+var+"  "+reviewAspectRatings[di][ai]+"  "+itemCodes[di][ai]);
            }
            for(int ai = 0;ai<numAspects;ai++){ 
            UserCodeGradient[ai] = -(lambda+var+(1/(alpha*userCodes[di][ai])));
            UserCodeGradient[ai] = Math.max(UserCodeGradient[ai],UserCodeGradient[ai]-rho);
            
            }
           
            double alphaUserCodeArmijo = getAlphaUserCodeArmijo(userCodes[di],reviewAspectRatings[di],itemCodes[di],UserCodeGradient);
            for(int ai = 0;ai<numAspects;ai++){
                userCodes[di][ai] += UserCodeGradient[ai]*alphaUserCodeArmijo;
            }
        }
        
        //4: for each item
        
        for(int h=0;h<numItems;h++ ){
            for(int ai = 0;ai<numAspects;ai++){
                double hotel = -(1/(alpha*alpha*userCodes[h][ai]*userCodes[h][ai]) * 
                        (reviewAspectRatings[h][ai]-itemCodes[h][ai]));
                ItemCodeGradient[ai] = hotel;
                
            }
            double alphaItemCodeArmijo = getAlphaUserCodeArmijo(reviewAspectRatings[h],itemCodes[h],userCodes[h],ItemCodeGradient);
            for(int ai = 0;ai<numAspects;ai++){
                itemCodes[h][ai] += UserCodeGradient[ai]*alphaItemCodeArmijo;
            }
        
        }
        dict_learn();
        
        // alpha part
        
        double alphaGradient;
        double variable = 0.0;
        
        for (int di=0; di<numDocs; di++){
            for(int ai=0;ai<numAspects;ai++){
                double var1 =  (reviewAspectRatings[di][ai] - itemCodes[di][ai])/userCodes[di][ai];
                double var2 = var1*var1;
                variable = var2/(1/Math.abs(numDocs)*Math.abs(numAspects));
                variable += variable;
            }
            alphaGradient = Math.sqrt(variable);
            alphaGradient+=alphaGradient;
            AlphaCodeGradient[di]=alphaGradient;
        }
          
         
    }
    double armijoBeta = 0.1;
    double stepSize = 0.8;
    double sigma = 0.8;
    
    
    public double getAlphaAspectRatingsArmijo(double overallRating,double[] aspectRatings, double[] aspectRatingGradient,
                                                double[] etad, double[] itemCode,double[] userCode){
        double alphaAspectRatings = 1.0;
        double newObjectiveValue;
        double[] newAspectRatings =new double[numAspects];
        double oldObjectiveValue = getAspectRatingObjectiveValue(overallRating,aspectRatings, etad,itemCode,userCode);
        double gradientInnerProduct = 0.0;
        for(int ai = 0;ai<numAspects;ai++)
            gradientInnerProduct += aspectRatingGradient[ai]*aspectRatingGradient[ai];
        for(int i =1;i<100;i++){
            for(int ai = 0;ai<numAspects;ai++)
                newAspectRatings[ai] = aspectRatings[ai] + Math.pow(armijoBeta,i)*stepSize*aspectRatingGradient[ai];
            newObjectiveValue = getAspectRatingObjectiveValue(overallRating,newAspectRatings, etad,itemCode,userCode);
            if(oldObjectiveValue-newObjectiveValue>=sigma*Math.pow(armijoBeta,i)*stepSize*gradientInnerProduct)
            {
                return Math.pow(armijoBeta,i)*stepSize;
            }
        }
        return alphaAspectRatings;
    }
    double getAspectRatingObjectiveValue(double overAllRating, double[] aspectRatings, double[] etad, 
                                            double[] itemCode,double[] userCode){
        double objective = 0.0;
        double ratingPrediction = 0.0;
        double userItemRatingFactor = 0.0;
        for(int ai = 0;ai<numAspects;ai++){
            ratingPrediction += aspectRatings[ai]*etad[ai];
            userItemRatingFactor += 1/(userCode[ai]*userCode[ai]) * ( aspectRatings[ai]-itemCode[ai]) * ( aspectRatings[ai]-itemCode[ai]);
        }
        objective = 1/(2*ratingVariance*ratingVariance) * (overAllRating-ratingPrediction)*(overAllRating-ratingPrediction) +
                1/(2*alpha*alpha)*userItemRatingFactor; 
        return objective;
    }
    
    
    public double getAlphaWordCodeArmijo(double[] wordCodes, double[] thetad,double[] betan, 
            double[] wordCountGradient, double wordCount){
        double alphaWordCode = 1.0;
        double newObjectiveValue;
        double[] newWordCodes =new double[numAspects];
        double oldObjectiveValue = getWordCodeObjectiveValue(wordCodes,thetad,betan,wordCount);
        double gradientInnerProduct = 0.0;
        for(int ai = 0;ai<numAspects;ai++)
            gradientInnerProduct += wordCountGradient[ai]*wordCountGradient[ai];
        for(int i =1;i<100;i++){
            for(int ai = 0;ai<numAspects;ai++)
                newWordCodes[ai] = wordCodes[ai] + Math.pow(armijoBeta,i)*stepSize*wordCodeGradient[ai];
            newObjectiveValue = getWordCodeObjectiveValue(newWordCodes,thetad,betan,wordCount);
            if(oldObjectiveValue-newObjectiveValue>=sigma*Math.pow(armijoBeta,i)*stepSize*gradientInnerProduct)
            {
                return Math.pow(armijoBeta,i)*stepSize;
            }
        }
        return alphaWordCode;
    }
    
    public double getWordCodeObjectiveValue(double[] wordCodes, double[] thetad, double[] betan,double wordCount){
        double objectiveValue = 0.0;
        double dotProduct = 0.0;
        for(int ai = 0;ai < numAspects;ai++){
            objectiveValue += gamma* (wordCodes[ai]-thetad[ai])*(wordCodes[ai]-thetad[ai]) + rho * Math.abs(wordCodes[ai]) ;
            dotProduct += wordCodes[ai]*betan[ai]; 
        }
        objectiveValue -=wordCount*Math.log(dotProduct) + dotProduct;
        return objectiveValue;
     }
    
    
    public double getAlphaUserCodeArmijo(double[] userCodes,double[] aspectRatings, double[] itemCodes,
            double[] UserCodeGradiant){
        double alphaUserCode = 1.0;
        double newObjectiveValue;
        double[] newUserCodes =new double[numAspects];
        double oldObjectiveValue = getUserCodeObjectiveValue(userCodes,aspectRatings,itemCodes);
        double gradientInnerProduct = 0.0;
        for(int ai = 0;ai<numAspects;ai++)
            gradientInnerProduct += UserCodeGradiant[ai]*UserCodeGradiant[ai];
        for(int i =1;i<100;i++){
            for(int ai = 0;ai<numAspects;ai++)
                newUserCodes[ai] = userCodes[ai] + Math.pow(armijoBeta,i)*stepSize*UserCodeGradient[ai];
            newObjectiveValue = getUserCodeObjectiveValue(newUserCodes,aspectRatings,itemCodes);
            if(oldObjectiveValue-newObjectiveValue>=sigma*Math.pow(armijoBeta,i)*stepSize*gradientInnerProduct)
            {
                return Math.pow(armijoBeta,i)*stepSize;
            }
        }
        return alphaUserCode;
    }
    
    double getUserCodeObjectiveValue(double[] userCode,double[] aspectRatings, double[] itemCodes){
        double objective = 0.0;
        double ratingPrediction = 0.0;
        double userItemRatingFactor = 0.0;
        double log_alpha_t = 0.0;
        
        for(int ai = 0;ai<numAspects;ai++){
            ratingPrediction += lambda*(Math.abs(userCode[ai]));
            userItemRatingFactor += 1/(userCode[ai]*userCode[ai]) * ( aspectRatings[ai]-itemCodes[ai]) * ( aspectRatings[ai]-itemCodes[ai]);
            log_alpha_t += Math.log(alpha*userCode[ai]);
        }
        objective = ratingPrediction + log_alpha_t +
                1/(2*alpha*alpha)*userItemRatingFactor; 
        return objective;
    }
    
    public double getAlphaItemCodeArmijo(double[] userCodes,double[] aspectRatings, double[] itemCodes,
            double[] ItemCodeGradiant){
        double alphaItemCode = 1.0;
        double newObjectiveValue;
        double[] newItemCodes =new double[numAspects];
        double oldObjectiveValue = getUserCodeObjectiveValue(userCodes,aspectRatings,itemCodes);
        double gradientInnerProduct = 0.0;
        for(int ai = 0;ai<numAspects;ai++)
            gradientInnerProduct += ItemCodeGradiant[ai]*ItemCodeGradiant[ai];
        for(int i =1;i<100;i++){
            for(int ai = 0;ai<numAspects;ai++)
                newItemCodes[ai] = itemCodes[ai] + Math.pow(armijoBeta,i)*stepSize*UserCodeGradient[ai];
            newObjectiveValue = getUserCodeObjectiveValue(newItemCodes,aspectRatings,itemCodes);
            if(oldObjectiveValue-newObjectiveValue>=sigma*Math.pow(armijoBeta,i)*stepSize*gradientInnerProduct)
            {
                return Math.pow(armijoBeta,i)*stepSize;
            }
        }
        return alphaItemCode;
    }
    
    double getItemCodeObjectiveValue(double[] userCode,double[] aspectRatings, double[] itemCodes){
        double objective = 0.0;
        double userItemRatingFactor = 0.0;
        
        for(int ai = 0;ai<numAspects;ai++){
            
            userItemRatingFactor += 1/(userCode[ai]*userCode[ai]) * ( aspectRatings[ai]-itemCodes[ai]) * ( aspectRatings[ai]-itemCodes[ai]);
            
        }
        objective = 1/(2*alpha*alpha)*userItemRatingFactor; 
        return objective;
    }
    double[] betaGradient; 
    long VAR_MAX_ITER = 100000L;
    double MIN_BETA = 1E-30;
    public boolean dict_learn() 
    {
        
	// run projected gradient descent.
	int opt_size = numAspects * vocabularySize;
	int opt_iter = 0, ndim = opt_size + 10;
	double pref = 1.0;
	double f=0, eps = 1.0e-5, xtol = 1.0e-16;
	int m = 5;
        int [] iprint = new int[]{-1,0};
        int [] iflag = new int[]{0};
	boolean diagco = false;
        double[] mu_ = new double[vocabularySize];
        double[] optimVar = new double[vocabularySize*numAspects];
        double[] diagParam = new double[opt_size+1];
        do  {
            opt_iter ++;
            pref = f;
            
            computeFGradientBeta();  
            //get_param(x_, m_nK, m_nNumTerms);
            //f = fdf_beta(pC, theta, s, g_);
            setOptimParam(optimVar);
            
            if (Math.abs(pref/f - 1.0) < 1e-4) break;
            try	{
                    //m_pLBFGS->lbfgs( opt_size, m, x_, f, g_, diagco, null, iprint, eps, xtol, iflag );
                    LBFGS.lbfgs(opt_size, m, optimVar, f, betaGradient, diagco, diagParam, iprint, eps, xtol, iflag);
                    //lbfgs ( n , m ,  x , f ,  g , bdiag,  diag , int[] iprint , double eps , double xtol , int[] iflag )
            } catch(LBFGS.ExceptionWithIflag e) {
                    System.out.printf("exception in l-bfgs\n");
                    break;
            }
            setBeta(optimVar);
            //set_param(x_, m_nK, m_nNumTerms);
            //printf("\t\t projection \n");
            for ( int ai=0; ai<numAspects; ai++ ) {
                for ( int wi=0; wi<vocabularySize; wi++ ) {
                        mu_[wi] = beta[wi][ai];
                }

                project_beta2( mu_,  1.0, MIN_BETA );

                for ( int wi=0; wi<vocabularySize; wi++ ) {
                    beta[wi][ai] = mu_[wi];
                }
            }
	} while ( iflag[0] != 0 && opt_iter < VAR_MAX_ITER );

	return true;
    }

    double computeFGradientBeta()//Corpus *pC, double **theta, double ***s, double *g)
    {	
        Arrays.fill(betaGradient,0.0);
        int nWrd, wcount, gIx = 0, d, n, k;
        Document doc;
        double fVal = 0, dVal = 0;
        //double *pS = NULL, *bPtr = NULL;
        //double **pSS = NULL;
        for ( d=0; d<numDocs; d++ ) {
            doc = docList.get(d);
            DocWordCode dwc = docWordCodes.get(d);
            //pDoc = (pC->docs[d]);
            for ( n=0; n<doc.counts.length; n++ ) {
                nWrd = doc.indices[n];   // word index
                wcount = doc.counts[n];  // word count
                
                double[] wcode = dwc.getWordCode(n).getTopicCode();
                dVal = m_dPoisOffset;
                for ( int ai=0; ai<numAspects; ai++ ) {
                        if (  wcode[ai]> 0 ) 
                                dVal += wcode[ai] * beta[nWrd][ai];
                }

                // update function value.
                fVal += dVal - wcount * Math.log(dVal);

                // update gradients.
                dVal = 1 - wcount/dVal;
                gIx = nWrd * numAspects;
                for ( k=0; k<numAspects; k++ ) {
                        if ( wcode[k] > 0 )
                                betaGradient[gIx] += dVal * wcode[k];
                        gIx ++;
                }
            }
        }

        return fVal;
    }

// project beta to simplex ( N*log(N) ).
    void project_beta( double [] vec )
    {
        double[] mu_ = Arrays.copyOf(vec,vec.length);
        Arrays.sort(mu_,0,vec.length-1);
        // find rho.
        int rho = 0;
        double dsum = 0;
        for ( int i=0; i<vec.length; i++ ) {
                dsum += mu_[i];

                if ( mu_[i] - (dsum-1)/(i+1) > 0 )
                        rho = i;
        }

        double theta = 0;
        for ( int i=0; i<=rho; i++ ) {
                theta += mu_[i];
        }
        theta = (theta-1) / (rho+1);
      
        for ( int i=0; i<vec.length; i++ ) {
                vec[i] = Math.max(0.0, vec[i] - theta);
        }
    }
    // linear algorithm
    void project_beta2( double[]vec, double dZ, double epsilon )
    {
        double[] mu_ = new double[vec.length];
            ArrayList<Integer> U = new ArrayList(vec.length);
            for ( int i=0; i<vec.length; i++ ) {
                    mu_[i] = vec[i] - epsilon;
                    U.set(i, i+1); ;
            }
            double dZVal = dZ - epsilon * vec.length; // make sure dZVal > 0

            /* project to a simplex. */
            double s = 0;
            int p = 0;
            Random rand = new Random();
            while ( U.size()!=0 ) {
                    int nSize = U.size();
                    
                    int k = U.get(rand.nextInt(nSize));

                    /* partition U. */
                    ArrayList<Integer> G = new ArrayList<>();
                    ArrayList<Integer> L = new ArrayList<>();
                    int deltaP = 0;
                    double deltaS = 0;
                    for ( int i=0; i<nSize; i++ ) {
                            int j = U.get(i);

                            if ( mu_[j-1] >= mu_[k-1] ) {
                                    if ( j != k ) G.add( j );
                                    deltaP ++;
                                    deltaS += vec[j-1];
                            } else L.add( j );
                    }

                    if ( s + deltaS - (p + deltaP) * mu_[k-1] < dZ ) {
                            s += deltaS;
                            p += deltaP;
                            U = L;
                    } else {
                            U = G;
                    }
            }

            double theta = (s - dZ) / p;
            for ( int i=0; i<vec.length; i++ ) {
                    vec[i] = Math.max(mu_[i] - theta, 0.0) + epsilon;
            }
    }
    public void setOptimParam(double[] optimParam){
        int oi = 0;
        for(int wi = 0;wi<vocabularySize; wi++){
            for(int ai = 0;ai<numAspects; ai++){
                optimParam[oi] = beta[wi][ai];
                oi++;
            }
            
        }
    }
    public void setBeta(double[] optimParam){
        int oi = 0;
        for(int wi = 0;wi<vocabularySize; wi++){
            for(int ai = 0;ai<numAspects; ai++){
                beta[wi][ai] = optimParam[oi];
                oi++;
            }
            
        }
    }
}

