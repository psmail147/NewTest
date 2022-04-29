package algorithms;

import java.util.Vector; 
//import java.util.Collections;
import static utils.algorithms.Misc.generateRandomSolution;

import interfaces.Algorithm;
import interfaces.Problem;
import utils.RunAndStore.FTrend;
import utils.algorithms.Misc;
import utils.algorithms.PSOOp;
import utils.random.RandUtils;

public class PUSO extends Algorithm
{
	public FTrend execute(Problem problem, int maxEvaluations) throws Exception
	{
		FTrend FT = new FTrend(); 
		
		int problemDimensions = problem.getDimension(); 
		double[][] bounds = problem.getBounds();
		
		int colonySize = getParameter("p0").intValue();  
		double inertia = getParameter("p1").doubleValue();
		double pInfluence = getParameter("p2").doubleValue();
		double gInfluence = getParameter("p3").doubleValue();
		
		//double gamma1 = getParameter("p4").doubleValue(); //rarely used ==> = 1
		//double gamma2 = getParameter("p5").doubleValue(); //rarely used ==> = 1
		
		double[][] colony = new double[colonySize][problemDimensions]; 
		double[][]localBestAnt = new double[colonySize][problemDimensions]; 
		double[][] velocity = new double[colonySize][problemDimensions]; 
		double[] globalBestAnt = new double[problemDimensions]; 
		double[] localBestAntFitness = new double[problemDimensions]; 

		double globalBestFitVal = Double.NaN;
		double newFitVal = Double.NaN; 
		int counter = 0; 

		
		for(int i = 0; i < colonySize; i++)
		{
			double[] temp = generateRandomSolution(bounds, problemDimensions);
			for (int j = 0; j < problemDimensions; j++)
			{
				colony[i][j] = temp[j];
				localBestAnt[i][j] = temp[j];
				velocity[i][j] = RandUtils.uniform(0, 0.1);
			}
			localBestAntFitness[i] = problem.f(colony[i]);
//			pFitness[j] = fitness[j];
			
			if (i == 0 || localBestAntFitness[i] < globalBestFitVal)
			{
				globalBestFitVal = localBestAntFitness[i];
				for (int n = 0; n < problemDimensions; n++)
					globalBestAnt[n] = colony[i][n];
					FT.add(counter, globalBestFitVal);
					counter++;
			}
			
			counter++;
		}
		
		
		while(counter < maxEvaluations)
		{
			double rand = 0.0;
			double x = 0.0; 
			
			if (counter < (maxEvaluations*0.5))
			{
				//PSO
				for (int i = 0; i < colonySize /*&& i < maxEvaluations*/; i++)
				{
					// update velocity 
					velocity[i] = PSOOp.classicVelocityUpdate(velocity[i], colony[i], localBestAnt[i], globalBestAnt, inertia, pInfluence, gInfluence);
					for (int j = 0; j < problemDimensions; j++)
					{
						//velocity[i][j] = velocity[i][j] * inertia + pInfluence*Math.random() *(localBestAnt[i][j]-colony[i][j]) + gInfluence*Math.random()*(globalBestAnt[i]-colony[i][j]); 
						colony[i][j] = colony[i][j] + velocity[i][j];
					}
					
					
					//Move
					colony[i] = Misc.toro(colony[i], bounds); 
					newFitVal = problem.f(colony[i]);
					counter++; 

					
					// update personal best position
					if (newFitVal < localBestAntFitness[i])
					{
						for (int j = 0; j < problemDimensions; j++)
							localBestAnt[i][j] = colony[i][j];
						localBestAntFitness[i] = newFitVal;
						
						// best update
						if (newFitVal < globalBestFitVal)
						{
							globalBestFitVal = newFitVal;
							for (int n = 0; n < problemDimensions; n++)
								globalBestAnt[n] = colony[i][n];
							FT.add(counter, globalBestFitVal);
						}
					}
					
				}
			}
			
			else if(counter < (maxEvaluations*0.75))
			{
				for (int i = 0; i < colonySize; i++)
				{
					for (int j = 0; j < problemDimensions; j++)
					{
					    rand = utils.random.RandUtils.uniform(colony[i][j], bounds[0][1] * 0.5);
						x = Math.round(Math.random()); 
						if (x == 0.0)
							colony[i][j] = colony[i][j] - rand;
						else
							colony[i][j] = colony[i][j] + rand; 
					}
					
					colony[i] = Misc.toro(colony[i], bounds); 
					newFitVal = problem.f(colony[i]);
					counter++; 

					
					// update personal best position
					if (newFitVal < localBestAntFitness[i])
					{
						for (int j = 0; j < problemDimensions; j++)
							localBestAnt[i][j] = colony[i][j];
						localBestAntFitness[i] = newFitVal;
						
						// best update
						if (newFitVal < globalBestFitVal)
						{
							globalBestFitVal = newFitVal;
							for (int n = 0; n < problemDimensions; n++)
								globalBestAnt[n] = colony[i][n];
							FT.add(counter, globalBestFitVal);
						}
					}
					else 
					{
						if (x == 0.0)
						{
							for (int n = 0; n < problemDimensions; n++)
								colony[i][n] = colony[i][n] + rand*2;
						}
						else
							for (int n = 0; n < problemDimensions; n++)
								colony[i][n] = colony[i][n] - rand*2;
						
							
									
					}
				}
			}
			
			else
			{
				double beta = 0.0; 

				for (int i = 0; i < colonySize; i++)
				{
					for (int j = 0; j < problemDimensions; j++)
					{
						double u = Math.random(); 
						double p = beta/1; 
						if (u < p)
							colony[i][j] = globalBestAnt[j] + rand; 
						else
							colony[i][j] = localBestAnt[i][j] + rand; 
					}
					beta = beta + 0.001; 
					//Check fitness. 
					colony[i] = Misc.toro(colony[i], bounds); 
					newFitVal = problem.f(colony[i]);
					counter++; 

					
					// update personal best position
					if (newFitVal < localBestAntFitness[i])
					{
						for (int j = 0; j < problemDimensions; j++)
							localBestAnt[i][j] = colony[i][j];
						localBestAntFitness[i] = newFitVal;
						
						// best update
						if (newFitVal < globalBestFitVal)
						{
							globalBestFitVal = newFitVal;
							for (int n = 0; n < problemDimensions; n++)
								globalBestAnt[n] = colony[i][n];
							FT.add(counter, globalBestFitVal);
						}
					}
				}
				counter++; 
			}

			
		}
		
		finalBest = globalBestAnt; 
		FT.add(counter, globalBestFitVal);
		return FT; 
	}

}
