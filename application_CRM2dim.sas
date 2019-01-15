options nonotes=0; 
/**************************************************************************/
/* Simulation 															  */
/**************************************************************************/

* Scenario 1 with 2-parameter model (nb_param=2) and no start-up (start_up=0); 
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.03/0.05/0.08/0.13/0.20/0.29/ 
							0.05/0.08/0.13/0.20/0.29/0.40/
							0.08/0.13/0.20/0.29/0.40/0.53,
			start_up=0,
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1,
			g_size1=,
			histo=,
			a0=);

* Scenario 1 with 2-parameter model (nb_param=2) and start-up (start_up=1);
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.03/0.05/0.08/0.13/0.20/0.29/ 
							0.05/0.08/0.13/0.20/0.29/0.40/
							0.08/0.13/0.20/0.29/0.40/0.53,
			start_up=1,
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* Scenario 1 with 3-parameter model (nb_param=3) and start-up (start_up=1);
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.03/0.05/0.08/0.13/0.20/0.29/ 
							0.05/0.08/0.13/0.20/0.29/0.40/
							0.08/0.13/0.20/0.29/0.40/0.53,	
			start_up=1,	
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=3,
			prior_dist=E/E/E,
			val_param=1/1/1,
			nbi=1000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* Scenario 1 with 3-parameter model (nb_param=3) and start-up (start_up=1) and incorporation of historical data for first agent;
* For the incorporation of historical data, we need to specify the dataset including the historical data (histo=histo_data);
* and the weight assigned to the historical data which is controlled with a0. a0=0 corresponds to no inclusion and a0=1 complete inclusion;

* Fictive historical data for the first agent;
* We assume available information (n=10 observations) at dose level 5 for the first agent and dose level 0 for the second agent with 2 toxicities;
data histo_data(drop=i);
	do i=1 to 8;
		patid=i;si=5;tj=0;y=0;
		output;
	end;
	do i=1 to 2;
		patid=8+i;si=5;tj=0;y=1;
		output;
	end;
run;

proc print noobs;run;

dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.03/0.05/0.08/0.13/0.20/0.29/ 
							0.05/0.08/0.13/0.20/0.29/0.40/
							0.08/0.13/0.20/0.29/0.40/0.53,
			start_up=1,
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,	
			dist=0,
			nb_param=3,	
			prior_dist=E/E/E,
			val_param=1/1/1,
			nbi=1000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=histo_data,
			a0=0.5);

* Scenario 2 with start-up (start_up=1);
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.04/0.05/0.08/0.11/0.15/0.21/ 
							0.04/0.06/0.09/0.13/0.18/0.25/
							0.05/0.08/0.11/0.15/0.21/0.29,
			start_up=1,
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* Scenario 3 with start-up (start_up=1);
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.03/0.05/0.13/0.20/0.27/0.35/ 
							0.10/0.20/0.25/0.32/0.41/0.50/
							0.20/0.30/0.41/0.53/0.65/0.70,	
			start_up=1,		
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1,		
			g_size1=2,			
			histo=,
			a0=);

* Scenario 4 with start-up (start_up=1);
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.03/0.05/0.08/0.13/0.17/0.20/ 
							0.08/0.13/0.20/0.32/0.41/0.50/
							0.20/0.40/0.47/0.56/0.65/0.76,			
			start_up=1,	
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1,		
			g_size1=2,			
			histo=,
			a0=);

* Scenario 5 with start-up (start_up=1);
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.08/0.13/0.20/ 
							0.08/0.13/0.20/0.29/
							0.13/0.20/0.29/0.40/
							0.20/0.29/0.40/0.53,			
			start_up=1,
			n_patient=54,
			group_size=3,
			nb_simul=4000,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3,
			_b_=0.05/0.1/0.2/0.3,
			dist=0,
			nb_param=2,	
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* Scenario 6 which comes from scenario 1 of Riviere's paper (Table 1) for illustrating the different model specifications of the SAS macro CRM2dim;
* We assume a start-up (start_up=1) and no incorporation of historical data;

* With 2 parameter power model from Wang and Ivanova - exponential distribution of variance 1;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,								
		    start_up=1,
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=0,
			nb_param=2,	
			prior_dist=E/E,
			val_param=1/1,
			nbi=2000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 2 parameter power model from Wang and Ivanova - exponential distribution of variance 10;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,	
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=0,
			nb_param=2,	
			prior_dist=E/E,
			val_param=10/10,
			nbi=2000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 3 parameter power model from Wang and Ivanova - exponential distribution of variance 1;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,		
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=0,
			nb_param=3,	
			prior_dist=E/E/E,
			val_param=1/1/1,
			nbi=2000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 3 parameter power model from Wang and Ivanova - exponential distribution of variance 10;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,	
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=0,
			nb_param=3,	
			prior_dist=E/E/E,
			val_param=10/10/10,
			nbi=2000,
			nmc=5000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 3 parameter standard logistic regression model from Riviere et al -  normal distribution for intercept with variance 1 and ;
* exponential distribution for dose effect of variance 1;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,	
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=1,
			nb_param=3,	
			prior_dist=N/E/E,
			val_param=1/1/1,
			nbi=50000,
			nmc=50000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 3 parameter standard logistic regression model from Riviere et al -  normal distribution for intercept with variance 10 and ;
* exponential distribution for dose effect of variance 10;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=1,
			nb_param=3,	
			prior_dist=N/E/E,
			val_param=10/10/10,
			nbi=50000,
			nmc=50000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 3 parameter standard logistic regression model from Riviere et al -  normal distribution for intercept with variance 50 and ;
* exponential distribution for dose effect of variance 10;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,	
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=1,
			nb_param=3,	
			prior_dist=N/E/E,
			val_param=50/10/10,
			nbi=50000,
			nmc=50000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 4 parameter standard logistic regression model from Riviere et al -  normal distribution for intercept with variance 1 and ;
* exponential distribution for dose effect of variance 1;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=1,
			nb_param=4,	
			prior_dist=N/E/E/N,
			val_param=1/1/1/1,
			nbi=50000,
			nmc=50000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 4 parameter standard logistic regression model from Riviere et al -  normal distribution for intercept with variance 10 and ;
* exponential distribution for dose effect of variance 1;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=1,
			nb_param=4,	
			prior_dist=N/E/E/N,
			val_param=10/1/1/10,
			nbi=50000,
			nmc=50000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

* With 4 parameter standard logistic regression model from Riviere et al -  normal distribution for intercept with variance 50 and ;
* exponential distribution for dose effect of variance 10;
dm "log; clear; ";
%CRM2dim(	simul=1,
			scenario=       0.05/0.10/0.15/0.30/0.45/ 
							0.10/0.15/0.30/0.45/0.55/
							0.15/0.30/0.45/0.50/0.60,	
			start_up=1,	
			n_patient=60,
			group_size=3,
			nb_simul=2000,
			gamma=0.30,
			_a_=0.12/0.2/0.3/0.4/0.5,
			_b_=0.2/0.3/0.4,
			dist=1,
			nb_param=4,	
			prior_dist=N/E/E/N,
			val_param=50/10/10/50,
			nbi=50000,
			nmc=50000,
			seed=1,
			g_size1=2,
			histo=,
			a0=);

/**************************************************************************/
/* Conducting a dual-agent Bayesian CRM trial                             */
/**************************************************************************/

* Suppose we are running a two-agent phase I trial with 6 levels;
* of the first agent (max_level1=6) and 3 levels of the second agent (max_level2=3);
* with the target probability of toxicity of 0.20 (gamma=0.2), two patients per cohort in the start-up;
* and three patients per cohort in the two-dimensional CRM;

* Suppose that the start-up produced the following sequence of combinations;
* (1,1), (2,1), (3,1), (4,1), (5,1), (6,1), (4,2), (5,2), (3,3) and (4,3);
* and one toxicity was observed at the combination (5,2) and another at combination (4,3); 
* This information is collected in the dataset current_data;

dm "log; clear; ";
data current_data;
	input patid si tj y;
	cards;
	1 1	1	0
	2 1 1	0
	3 2	1	0
	4 2	1	0
	5 3 1	0
	6 3	1	0
	7 4 1   0
	8 4 1   0 
	9 5 1   0
	10 5 1  0
	11 6 1 0
	12 6 1 0
	13 4 2 0
	14 4 2 0
	15 5 2 1
	16 5 2 0
	17 3 3 0
	18 3 3 0
	19 4 3 0
	20 4 3 1
	;
run;

proc print noobs;run;

* The recommended dose combination based on information accumulated from the start-up data (simul=0) with a two-parameter model (nb_param=2) ;
* and no historical data (histo and a0 not specified) is returned by the call of CRM2dim;

%CRM2dim(	simul=0,
			start_up=1,
			n_patient=35,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1);

* -> Next dose combination recommendation: (5,1) with probability close to 0.20 (0.17639);
* Although there are dose levels with toxicity probability close to 0.20, we select the smallest level of agent 2 at the end of the startup ;

* The next step is the start of the two-dimensional dose escalation design;
* Assuming that three patients (#21,22,23) are treated at the combination (5,1) and 2 toxicity are observed;
* The accumulated data (from startup + 3 patients) is now updated in the current_data SAS dataset;

data current_data;
	input patid si tj y;
	cards;
	1 1	1	0
	2 1 1	0
	3 2	1	0
	4 2	1	0
	5 3 1	0
	6 3	1	0
	7 4 1   0
	8 4 1   0 
	9 5 1   0
	10 5 1  0
	11 6 1 0
	12 6 1 0
	13 4 2 0
	14 4 2 0
	15 5 2 1
	16 5 2 0
	17 3 3 0
	18 3 3 0
	19 4 3 0
	20 4 3 1
	21 5 1 0
	22 5 1 1
	23 5 1 1
	;
run;
proc print noobs;run;

* To define the next recommended dose, we call CRM2dim with simul=0 and start_up=0;
* The exploration of the next dose levels is performed around the current dose level with no skipping rule;
* The skipping rule allows to avoid to treat patients at toxic dose level;

%CRM2dim(	simul=0,
			start_up=0,
			n_patient=35,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1);
			
* -> Next dose combination recommendation: (4,2) with probability close to 0.20 (0.18793);

* Three new patients (#24,25,26) are treated at this recommended dose and assuming that one patient has a toxicity;
data current_data;
	input patid si tj y;
	cards;
	1 1	1	0
	2 1 1	0
	3 2	1	0
	4 2	1	0
	5 3 1	0
	6 3	1	0
	7 4 1   0
	8 4 1   0 
	9 5 1   0
	10 5 1  0
	11 6 1 0
	12 6 1 0
	13 4 2 0
	14 4 2 0
	15 5 2 1
	16 5 2 0
	17 3 3 0
	18 3 3 0
	19 4 3 0
	20 4 3 1
	21 5 1 0
	22 5 1 1
	23 5 1 1
	24 4 2 0
 	25 4 2 0
	26 4 2 1
	;
run;
proc print noobs;run;

* To identify the new recommended dose based on these 26 patients, we call CRM2dim;

%CRM2dim(	simul=0,
			start_up=0,
			n_patient=35,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1);

* -> Next dose combination recommendation: (4,2) with probability close to 0.20 (0.20192);

* Three new patients (#27,28,29) are treated at this recommended dose and we assume that no toxicity is observed;
			
data current_data;
	input patid si tj y;
	cards;
	1 1	1	0
	2 1 1	0
	3 2	1	0
	4 2	1	0
	5 3 1	0
	6 3	1	0
	7 4 1   0
	8 4 1   0 
	9 5 1   0
	10 5 1  0
	11 6 1 0
	12 6 1 0
	13 4 2 0
	14 4 2 0
	15 5 2 1
	16 5 2 0
	17 3 3 0
	18 3 3 0
	19 4 3 0
	20 4 3 1
	21 5 1 0
	22 5 1 1
	23 5 1 1
	24 4 2 0
 	25 4 2 0
	26 4 2 1
	27 4 2 0
	28 4 2 0
	29 4 2 0
	;
run;
proc print noobs;run;

%CRM2dim(	simul=0,
			start_up=0,
			n_patient=35,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1);

* -> Next dose combination recommendation: (3,3) with probability close to 0.20 (0.21135);

* Three new patients (#30,31,32) are treated at this recommended dose and no toxicity is observed;
			
data current_data;
	input patid si tj y;
	cards;
	1 1	1	0
	2 1 1	0
	3 2	1	0
	4 2	1	0
	5 3 1	0
	6 3	1	0
	7 4 1   0
	8 4 1   0 
	9 5 1   0
	10 5 1  0
	11 6 1 0
	12 6 1 0
	13 4 2 0
	14 4 2 0
	15 5 2 1
	16 5 2 0
	17 3 3 0
	18 3 3 0
	19 4 3 0
	20 4 3 1
	21 5 1 0
	22 5 1 1
	23 5 1 1
	24 4 2 0
 	25 4 2 0
	26 4 2 1
	27 4 2 0
	28 4 2 0
	29 4 2 0
	30 3 3 0
	31 3 3 0
	32 3 3 0
	;
run;
proc print noobs;run;

%CRM2dim(	simul=0,
			start_up=0,
			n_patient=35,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1);

* -> Next dose combination recommendation: (4,3) with probability close to 0.20 (0.21692);

* Three new patients (#33,34,35) are treated at this recommended dose and 3 toxicity are observed;
			
data current_data;
	input patid si tj y;
	cards;
	1 1	1  0
	2 1 1  0
	3 2	1  0
	4 2	1  0
	5 3 1  0
	6 3	1  0	
	7 4 1  0
	8 4 1  0 
	9 5 1  0
	10 5 1 0
	11 6 1 0
	12 6 1 0
	13 4 2 0
	14 4 2 0
	15 5 2 1
	16 5 2 0
	17 3 3 0
	18 3 3 0
	19 4 3 0
	20 4 3 1
	21 5 1 0
	22 5 1 1
	23 5 1 1
	24 4 2 0
 	25 4 2 0
	26 4 2 1
	27 4 2 0
	28 4 2 0
	29 4 2 0
	30 3 3 0
	31 3 3 0
	32 3 3 0
	33 4 3 1
	34 4 3 1
	35 4 3 1
	;
run;
proc print noobs;run;

* The total number of patients is included (n=35);
* The next execution of the macro allows to obtain the final recommended dose;

%CRM2dim(	simul=0,
			start_up=0,
			n_patient=35,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1);

* -> Next dose combination recommendation: (4,1), (4,2), (1,3);

* The code below generate a plot which summarizes the different escalation and de-escalation during the trial;
* (See figure 2 in the paper);

ods output list=list;
proc freq data=current_data;
	tables si*tj*y / list;
run;
ods output;

data list;
	set list;
	if si=2 and tj=1 then do;si_old=1;tj_old=1;end;
	if si=3 and tj=1 then do;si_old=2;tj_old=1;end;
	if si=4 and tj=1 then do;si_old=3;tj_old=1;end;
	if si=5 and tj=1 then do;si_old=4;tj_old=1;end;
	if si=6 and tj=1 then do;si_old=5;tj_old=1;end;
	if si=4 and tj=2 then do;si_old=6;tj_old=1;end;
	if si=5 and tj=2 then do;si_old=4;tj_old=2;end;
	if si=3 and tj=3 then do;si_old=5;tj_old=2;end;
	if si=4 and tj=3 and y=0 then do;si_old=3;tj_old=tj-0.025;tj=tj-0.025;end;

	if y=1 then do; 
			if si=5 and tj=2 then do;x_tox_st=si;y_tot_st=tj;end;
			if si=4 and tj=3 then do;x_tox_st=si;y_tot_st=tj-0.025;end;

			if si=5 and tj=1 then do;x_tox_st=si;y_tot_st=tj;end;
			
			if si=4 and tj=2 then do;x_tox_st=si;y_tot_st=tj;end;
		
	end;

	if si=4 and tj=2 then do;si_1=si;si_old_2=5;tj_1=tj;tj_old_2=1;end;
	if si=3 and tj=3 then do;si_1=si;si_old_2=4;tj_1=tj;tj_old_2=2;end;
	if si=4 and tj=3 and y=1 then do;si_1=si;si_old_2=3;tj_1=tj+0.025;tj_old_2=tj+0.025;end;

	if si=4 and tj=1 then do;si_RD=si;tj_RD=tj;end;
	if si=4 and tj=2 then do;si_RD=si;tj_RD=tj;end;
run;
data list2;
	set list end=eof;
	output;
	if eof then do;
       si_RD=1;tj_RD=3;output;
	   x_tox_st=5+0.1;y_tot_st=1;output;
	   x_tox_st=4;y_tot_st=3+0.025;output;
	   x_tox_st=4+0.1;y_tot_st=3+0.025;output;
	   x_tox_st=4+0.2;y_tot_st=3+0.025;output;
    end;
run;

ods graphics / reset noborder imagename="Flow_chart_running_phase1" imagefmt=jpeg /*BMP, GIF, JPEG, JPG, PDF, PS, TIFF*/;
ods listing gpath="X:\two_dimensional_CRM";
title "Flow chart of two-dimensional design with startup";
proc sgplot data=list2;
    * startup;
	vector x=si y=tj / xorigin=si_old yorigin=tj_old lineattrs=(thickness=1px pattern=dot color=black) legendlabel="startup";
	* two dim;
	vector x=si_1 y=tj_1 / xorigin=si_old_2 yorigin=tj_old_2 lineattrs=(thickness=2px pattern=dash color=black) legendlabel="two-dimensional";
	* DLT;
	scatter x=x_tox_st y=y_tot_st / MARKERATTRS=(color=red symbol=starfilled) legendlabel="dose limiting toxicity";
	* Recommended dose;
	scatter x=si_RD y=tj_RD / MARKERATTRS=(color=green symbol=circle size=30) legendlabel="recommended dose";

	xaxis label="Dose level of the first agent"  values=(1 to 6);
	yaxis label="Dose level of the second agent" values=(1 to 3);
	
run;
ods graphics off;
ods listing;

/**************************************************************************/
/* Conducting a dual-agent Bayesian CRM trial                             */
/* Assuming no startup                                                    */
/* This example is not reported in the manuscript                         */
/**************************************************************************/

* Suppose we are running a two-agent phase I trial with 6 levels of the first agent (max_level1=6) and 3 levels of the second agent (max_level2=3);
* with the target probability of toxicity of 0.20 (gamma=0.2) and three patients per cohort in the two-dimensional CRM;
* We assume no startup and the first 3 patients are treated at the current dose level (1,1) with no toxicity;
* The last information is collected in the current_data dataset;
* To specify no startup the macro variable start_up=0;

data current_data;
	input patid si tj y;
	cards;
	1 1	1	0
	2 1 1	0
	3 2	1	0
	;
run;

proc print noobs;run;

%CRM2dim(	simul=0,
			start_up=0,
			n_patient=35,
			gamma=0.20,
			_a_=0.05/0.1/0.2/0.3/0.5/0.7,
			_b_=0.05/0.1/0.2,		
			dist=0,
			nb_param=2,
			prior_dist=E/E,
			val_param=1/1,
			nbi=1000,
			nmc=5000,
			seed=1);

* -> Next dose combination recommendation: (1,2) with probability close to 0.20 (0.13181) among the possible dose levels;
* around the current dose level (1,1) since the rule of no skipping is implemented;
* This rule means that we cannot escalade from (1,1) to (2,2);

