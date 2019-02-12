*-----------------------------------------------------------------------------------------------------------------------;
*            SAS MACRO "CRM2dim"                                      													;
*                                                                     													;
* TITLE:     Implementation of Wang and Ivanova two-dimensional continual reassessment method in discrete dose space    ;       													     													;
*            for dose-finding in oncology Phase I trial (Biometrics Vol. 61, pp. 217-222 2005)                   		;	 
* AUTHOR:    Gwénaël Le Teuff                                         													;
*            Dept. of Biostatistics and Epidemiology                  													;
*            CESP INSERM U1018, OncoStat team                         													;
*            Gustave Roussy										      													;
*            114 rue Ed Vaillant 94805 Villejuif cedex - Paris        													;
* VERSION:   Version 1.1                                              													;
* DATE:      February 2019                                            												    ;
* COMMENT:   This program uses DATA steps and Base-SAS procedures and procedure MCMC for Bayesian analysis  			;
*-----------------------------------------------------------------------------------------------------------------------;
* HOW TO USE THE SAS MACRO "CRM2dim"                                 													;
* The main SAS macro CRM2dim calls 3 sub-macros %ESTIMATE, %DEFINE_DOSE and %REPORT                                     ;
*																														;
* The parameters of CRM2dim are described below                        													;
*																														;
* SIMUL       specifies if a simulation study is performed (simul=1) or not (simul=0). In the latter case, the data     ;
*             collected during a trial is continuously pooled by the user in a SAS data set call current_data           ;                                            ;
* SCENARIO    specifies the true toxicity probabilities at each dose combination                                    	;
* START_UP    flag to indicate whether a start-up step is run. Values are 0 or 1                                        ;
* N_PATIENT   total number of patients in the trial including patients from the start-up step is specified        		;
* GROUP_SIZE  number of patients per group assigned in the two-dimensional CRM											;
* NB_SIMUL    number of replications run for a simulation                                                        		;
* GAMMA       targeted probability of toxicity between 0 and 1. Usually, the targeted probability is 0.2 or 0.3         ;
* _A_         vector representing a set of constant ai(i=1,...J) separated by / with I being the number of dose levels 	;
*             of the first agent   																						;
* _B_         vector representing a set of constant bj(j=1,...J) separared by / with J being the number of dose levels  ;
*             of the second agent																						;
* DIST        Model specification to modeling the probability of toxicity (power or logistic)                           ;                               
* NB_PARAM    Number of model parameters																 				;
* PRIOR_DIST  Prior distribution of model parameters                                                                    ;
* VAL_PARAM   Parameters values of prior distribution                                                                   ;
* NBI         Number of burn-in iterations for Bayesian analysis														;
* NMC         Number of MCMC iterations, excluding the burn-in iterations												;
* SEED        Random seed for MCMC iterations																			;
* G_SIZE1     optional parameter specified when the macro variable start_up=1. This is the number of patients per cohort;
*             in the start-up 									                                                        ;
* HISTO       optional parameter indicating the name of SAS dataset containing historical data to include 				;
* A0          optional parameter between 0 and 1 representing the weight assigned to historical data  					;
*-----------------------------------------------------------------------------------------------------------------------;

*-----------------------------------------------------------------------------------------------------------------------;
/* The macro estimate estimates the toxicity posterior probability using the SAS proc MCMC                             */
%macro ESTIMATE(in=);

proc mcmc data=&in nbi=&nbi nmc=&nmc seed=&seed plots=none outpost=_out monitor=(p_post);
   array a(&max_level1);
   array b(&max_level2);

   /*******************************************************************************/
   /* Power model used by Wang and Ivanova */
  %if &dist=0 %then %do;
 
	   /* 2-parameters model */
	  %if &nb_param=2 %then %do;

	   	   parms alpha 1;parms beta 1;   
	   		
		   beginprior;
				array p_post(&max_level2,&max_level1);
				do _j_=1 to &max_level2;
		        	do _i_=1 to &max_level1;				    
		           		p_post[_j_,_i_]=1-((1-a[_i_])**alpha)*((1-b[_j_])**beta);					
		         	end;  
		        end;
		   endprior; 		   
	   	  
		   prior alpha ~ %if %scan(&prior_dist,1,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,1,"/"))))) %end;;
		   prior beta  ~ %if %scan(&prior_dist,2,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,2,"/"))))) %end;;		   

		   llike=log((1-((1-a[si])**alpha)*((1-b[tj])**beta))**y) + log((((1-a[si])**alpha)*((1-b[tj])**beta))**(1-y));
			%if &histo ne  %then %do;
				if (group = "D0") then	llike = &a0 * (log((1-((1-a[si])**alpha))**y) + log(((1-a[si])**alpha)**(1-y)));
			%end;

	   %end;

	   /* 3-parameters model */
	   %else %if &nb_param=3 %then %do; 

		   parms alpha 1; parms beta 1; parms gamma 1;

		   beginprior;
				array p_post(&max_level2,&max_level1);
				do _j_=1 to &max_level2;
		        	do _i_=1 to &max_level1;
					    p_post[_j_,_i_]=1-((1-a[_i_])**alpha)*((1-b[_j_])**(beta-gamma*log(1-a[_i_])));
		         	end;  
		        end;
		   endprior; 

		   prior alpha ~ %if %scan(&prior_dist,1,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,1,"/"))))) %end;;
		   prior beta  ~ %if %scan(&prior_dist,2,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,2,"/"))))) %end;;		   
           prior gamma ~ %if %scan(&prior_dist,3,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,3,"/"))))) %end;;
	   			        
		   llike=log((1-((1-a[si])**alpha)*((1-b[tj])**(beta-gamma*log(1-a[si]))))**y) + log((((1-a[si])**alpha)*((1-b[tj])**(beta-gamma*log(1-a[si]))))**(1-y));
			%if &histo ne %then %do;
				if (group = "D0") then llike = &a0 * (log((1-((1-a[si])**alpha))**y) + log(((1-a[si])**alpha)**(1-y)));
		    %end;

		%end;

   %end;

   /*******************************************************************************/
   /* Logistic model used by from Riviere et al */
  %else %if &dist=1 %then %do;

	   %if &nb_param=3 %then %do;

	   	   parms beta0 0 beta1 1 beta2 1; /* putting parameters in the same block improves the mixing of the chain when autocorrelations are observed */

		   array u(&max_level1);
		   array v(&max_level2);
		   do _i_=1 to &max_level1;
				u[_i_]=log(a[_i_]/(1-a[_i_]));
		   end;
		   do _j_=1 to &max_level2;
				v[_j_]=log(b[_j_]/(1-b[_j_]));
		   end;

		   beginprior;
				array p_post(&max_level2,&max_level1);
				do _j_=1 to &max_level2;
		        	do _i_=1 to &max_level1;			       						
						p_post[_j_,_i_]=exp(beta0 + beta1*u[_i_] + beta2*v[_j_])/(1+exp(beta0 + beta1*u[_i_] + beta2*v[_j_]));
		         	end;  
		        end;
		   endprior; 

		   prior beta0 ~ %if %scan(&prior_dist,1,"/")=N %then %do;normal(0,var=%scan(&val_param,1,"/"))   								%end;;
		   prior beta1 ~ %if %scan(&prior_dist,2,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,2,"/"))))) 	%end;;		   
           prior beta2 ~ %if %scan(&prior_dist,3,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,3,"/"))))) 	%end;;
		
		   num=exp(beta0+beta1*u[si]+beta2*v[tj]); den=1+exp(beta0+beta1*u[si]+beta2*v[tj]);
		   pij=num/den;
		   llike=log(pij**y) + log((1-pij)**(1-y));
		   %if &histo ne %then %do;
				if (group = "D0") then do;
						num_histo=exp(beta0+beta1*u[si]); den_histo=1+exp(beta0+beta1*u[si]);
		   				pij_histo=num_histo/den_histo;
		   				llike=&a0 * (log(pij_histo**y) + log((1-pij_histo)**(1-y)));						
				end;
		    %end;

	   %end;

	   %else %if &nb_param=4 %then %do;

		   parms beta0 0 beta1 1 beta2 1 beta3 0; /* putting parameters in the same block improves the mixing of the chain when autocorrelations are observed */

		   array u(&max_level1);
		   array v(&max_level2);
		   do _i_=1 to &max_level1;
				u[_i_]=log(a[_i_]/(1-a[_i_]));
		   end;
		   do _j_=1 to &max_level2;
				v[_j_]=log(b[_j_]/(1-b[_j_]));
		   end;

		   beginprior;
				array p_post(&max_level2,&max_level1);
				do _j_=1 to &max_level2;
		        	do _i_=1 to &max_level1;			       						
						p_post[_j_,_i_]=exp(beta0 + beta1*u[_i_] + beta2*v[_j_] + beta3*u[_i_]*v[_j_])/(1+exp(beta0 + beta1*u[_i_] + beta2*v[_j_] + beta3*u[_i_]*v[_j_]));
		         	end;  
		        end;
		   endprior; 

		   prior beta0 ~ %if %scan(&prior_dist,1,"/")=N %then %do;normal(0,var=%scan(&val_param,1,"/"))      							%end;;
		   prior beta1 ~ %if %scan(&prior_dist,2,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,2,"/"))))) 	%end;;		   
           prior beta2 ~ %if %scan(&prior_dist,3,"/")=E %then %do;expon(iscale = %sysevalf(1/%sysfunc(sqrt(%scan(&val_param,3,"/"))))) 	%end;;
		   prior beta3 ~ %if %scan(&prior_dist,4,"/")=N %then %do;normal(0,var=%scan(&val_param,4,"/"))      							%end;;
		  
		   num=exp(beta0+beta1*u[si]+beta2*v[tj]+beta3*u[si]*v[tj]); den=1+exp(beta0+beta1*u[si]+beta2*v[tj]+beta3*u[si]*v[tj]);
		   pij=num/den;
		   llike=log(pij**y) + log((1-pij)**(1-y));
		   %if &histo ne %then %do;
				if (group = "D0") then do;
						num_histo=exp(beta0+beta1*u[si]); den_histo=1+exp(beta0+beta1*u[si]);
		   				pij_histo=num_histo/den_histo;
		   				llike=&a0 * (log(pij_histo**y) + log((1-pij_histo)**(1-y)));						
				end;
		    %end;

	   %end;

   %end;
 
   model general(llike);  
run;


ods output summary=sum(keep=p_post:);
proc means data=_out mean;
run;
ods output close;

proc datasets nolist;
   delete _out;
run; quit;

%mend ESTIMATE;

*-----------------------------------------------------------------------------------------------------------------------;
/* The macro define_dose identifies the recommended dose combination                                                   */
%macro DEFINE_DOSE(in=,mode=);

data &in;
	row=1;
	set &in end=eof;

	array a(&max_level1);
	array b(&max_level2);
	array prob(&max_level2,&max_level1);
	array diff(&max_level2,&max_level1);

	if eof then do;

	set sum point=row;

	%let k=0;
   	%do _j_=1 %to &max_level2;
		%do _i_=1 %to &max_level1;
			%let k=%eval(&k+1);
		    prob(&_j_,&_i_)=p_post&k._mean;
			diff(&_j_,&_i_)=abs(prob(&_j_,&_i_)-&gamma);
		%end;  
	%end;

/*---------------------------------------------------------------*/
	%if &mode ne 2dim_design %then %do;
		%do j=1 %to &max_level2;
			init=diff(&j,1);indice&j=1;
			%let i=2;
			%do %while (&i <=&max_level1);
				if diff(&j,&i) < init then do;init=diff(&j,&i);indice&j=&i;end;
				%let i=%eval(&i+1);
			%end;
		%end;
		drop init;
	%end;
	end;
run;

proc datasets lib=work nolist;
	delete sum;
run;quit;

/*---------------------------------------------------------------*/
%if &mode=startup %then %do;
		/* the smallest level of agent 2 is choosed */
		%let current_tj=1;
		data _null_;
			set &in end=eof;
			if eof then call symput("current_si",indice1);
		run;		
%end; 
/*---------------------------------------------------------------*/
%if &mode=2dim_design %then %do;
		data &in;
			set &in end=eof;
			array diff(&max_level2,&max_level1);
			array near(&max_level2,&max_level1);
			array search(&max_level2,&max_level1);

			if eof then do;
			%do j=1 %to &max_level2;
				%do i=1 %to &max_level1;
					near(&j,&i)=99999; /* arbitrary value uses to identify the dose level which are not to explore from the current dose level */
				%end;
			%end;
			
			/* near=0 identify the dose level which can be explored from the current dose level */
			near(%eval(&current_tj),&current_si)=0;
			%if %eval(&current_tj-1) ge 1   %then %do;near(%eval(&current_tj-1),&current_si)=0;%end;
			%if %eval(&current_tj+1) le &max_level2 %then %do;near(%eval(&current_tj+1),&current_si)=0;%end;
			%if %eval(&current_si-1) ge 1 %then %do;near(&current_tj,%eval(&current_si-1))=0;%end;
			%if %eval(&current_tj+1) le &max_level2 and %eval(&current_si-1) ge 1 %then %do;near(%eval(&current_tj+1),%eval(&current_si-1))=0;%end;
			%if %eval(&current_si+1) le &max_level1 %then %do;near(&current_tj,%eval(&current_si+1))=0;%end;
			%if %eval(&current_tj-1) ge 1 and %eval(&current_si+1) le &max_level1 %then %do;near(%eval(&current_tj-1),%eval(&current_si+1))=0;%end;

			%do j=1 %to &max_level2;
				%do i=1 %to &max_level1;
					search(&j,&i)=diff(&j,&i)+near(&j,&i);
				%end;
			%end;

			init=99999; /* arbitrary value to initialize the research of the smallest difference */
			%do j=1 %to &max_level2;				
				%do i=1 %to &max_level1;
					if search(&j,&i) < init then do;init=search(&j,&i);indj=&j;indi=&i;end;
				%end;
			%end;
			call symput("current_tj",indj);	
			call symput("current_si",indi);	
			end;
			drop init;
		run;	
%end;
/*---------------------------------------------------------------*/
%if &mode=recommend_dose %then %do;
		data _null_;
			set &in end=eof;
			if eof then do;
				%do j=1 %to &max_level2;
					%global recom_dose&j;
					call symput("recom_dose&j",indice&j);
				%end;	
			end;
		run;		
%end;

%mend DEFINE_DOSE;

*-----------------------------------------------------------------------------------------------------------------------;
/* This macro report reports the results of simulation study                                                           */
%macro REPORT(in=,var=);

ods output summary=mean(keep=&var:);
ods html close; ods listing close;
proc means data=&in mean;
	var %do l=1 %to %sysevalf(&max_level2*&max_level1);
				&var.&l
		%end;
		;	
run;
ods html;ods listing;
ods output close;

data mean_&var.;
	set mean;
	
	%let k=0;
	%do j=1 %to &max_level2;
		%do i=1 %to &max_level1;
			%let k=%eval(&k+1);
			agentB=&j;
			%if &var=n_pat %then %do;agentA_level&i=((&var.&k._mean/&n_patient)*100);%end;
            %if &var=n_tox %then %do;agentA_level&i=&var.&k._mean;%end;
			%if &var=recommended_dose %then %do; agentA_level&i=&var.&k._mean*100;%end;
		%end;
		output;
	%end;
	keep agentB  %do i=1 %to &max_level1; agentA_level&i %end;;
run;

data mean_&var;
	length characteristic $16;
	set mean_&var;
	if _n_=1 then characteristic="&var";
run;

%mend REPORT;

*-----------------------------------------------------------------------------------------------------------------------;
/* The CRM2dim is the main macro for simulation or conducting a dual-agent Bayesian CRM phase 1 trial                  */
/* This macro calls the estimate, define_dose. The report macro is only executed for simulation study                  */ 
%macro CRM2DIM(	simul=,
				scenario=,
				start_up=,
				n_patient=, 
				group_size=,
				nb_simul=,
				gamma=,			
				_a_=,
				_b_=,
				dist=,
				nb_param=,
				prior_dist=,
				val_param=,				
				nbi=,
				nmc=,
				seed=,
				g_size1=,
				histo=,
				a0=													
			);

	%global current_tj current_si;

	%let max_level1=%sysfunc(countw(&_a_,"/"));
	%let max_level2=%sysfunc(countw(&_b_,"/"));

/****************************** conducting a CRM trial ********************************/
%if &simul=0 %then %do;

    %global first_run_startup first_run_2dim; /* Allow to identify a startup has be done and first 2dim cohort */

    %if &first_run_startup=%str() and &first_run_2dim=%str() %then %do;%let current_si=1;%let current_tj=1;%end; /* Initialisation to (1,1) the first cohort if no previous start-up */
   
	%put NOTE: 	conduction a dual agent CRM trial;

    ods _all_ close;

	data current_data ;
			set current_data;			
			/* Working model ai and bj */
			array a(&max_level1);
			array b(&max_level2);

			%do _i_=1 %to &max_level1;
				a(&_i_)=%scan(&_a_,&_i_,"/");		
			%end;

			%do _j_=1 %to &max_level2;
				b(&_j_)=%scan(&_b_,&_j_,"/");
			%end;
	run;

	* Compute the number of patients included in the actual phase 1;
	proc sql noprint;
		select count(*) into: nb_actual_patient
		from current_data;
	quit;

	%if &start_up=1 %then %do;
	    %let first_run_startup=1;
	    %put;
		%put Conducting a CRM trial - start-up;
		%ESTIMATE(in=current_data);
		%DEFINE_DOSE(in=current_data,mode=startup);
	%end;

	%else %if &start_up=0 %then %do;
		%let first_run_2dim=1;
        %if (&nb_actual_patient < &n_patient) %then %do;  			/* the total number of patient is not achieved. The dose escalation continue */
    		%put;
			%put Conducting a CRM trial - two dimensional design;
			%ESTIMATE(in=current_data);
			%DEFINE_DOSE(in=current_data,mode=2dim_design);
		%end;
		%else %do;													/* the total number of patients is achieved. The trial is finished */
			%put;
			%put Conducting a CRM trial - final recommended doses;
			%ESTIMATE(in=current_data);
			%DEFINE_DOSE(in=current_data,mode=recommend_dose);

			%let first_run_startup=;								/* Initialisation for the next trial of the current dose combination */
			%let first_run_2dim=;   
            %let current_si=;%let current_tj=;	
		%end;
	%end;
	
	data post_mean_current_data;
		set current_data end=eof;
		if eof then do;		
		%let k=0;
   		%do _j_=1 %to &max_level2;
			%do _i_=1 %to &max_level1;
				%let k=%eval(&k+1);
				agentB=&_j_;
				agentA_level&_i_=p_post&k._mean;  
			%end;
			output;
		%end;
		end;

	run;

	ods html; ods listing;
	%if (&nb_actual_patient < &n_patient) %then %do;
			footnote "Next dose combination recommendation: (%left(&current_si),%left(&current_tj))";
	%end;
	%else %do;
			footnote "Final dose combination recommendation: (%left(&recom_dose1),1), (%left(&recom_dose2),2), (%left(&recom_dose3),3)";
	%end;

	proc print data=post_mean_current_data noobs;
		var agentB agentA_:;
	run;
	footnote;

	dm 'ODSRESULTS' clear editor;

%end;

/****************************** Operating characteristics through simulations *********/
%else %if &simul=1 %then %do;

%put;
%put %str(NOTE	Parameters of the simulation study of a two-dimensional CRM);
%put;
%put %str(      Simulation under the scenario=&scenario);
%put %str(		Start-up phase included in the design=&start_up); 
%put %str(		Total patient number in the trial=&n_patient); 
%put %str(      Patient number per group assigned in the two-dimensional CRM=&group_size);
%put %str(		Number of replications in the simulation study=&nb_simul);
%put %str(		Targeted probability of toxicity=&gamma);
%put %str(      Model for drug combination-toxicity relationship=&dist);
%put %str(		Number of model parameters=&nb_param);
%put %str(		Prior distribution of parameters=&prior_dist);
%put %str(      Parameters values (variance) of prior distributions =&val_param);
%put %str(		Number of burn-in iterations for Bayesian analysis=&nbi);
%put %str(		Number of MCMC iterations, excluding the burn-in iterations=&nmc);
%put %str(		Random seed for MCMC iterations=&seed);
%put %str(		Patients number per cohort in the start-up=&g_size1);
%put %str(		SAS dataset containing historical data=&histo);
%put %str(		Weight assigns to the historical data=&a0);
%put;

%if &max_level1 < &max_level2 %then %put "Error: the number of dose for the first agent is lower than those for the second agent";
%let init_time = %sysfunc(datetime());%put Init time: &init_time;%put;

data pool;
		if (0);
run;

/*---------------------------------------------------------------*/
/* Starting of simulation -------------------------------------- */
%do rep=1 %to &nb_simul;
  
    ods _all_ close;
	
    %put 							SIMULATION: &rep;

	/* Important to define the current level of agent 1 (current_si) and agent 2 (current_tj) */
	/* for the 2 design escalation 															  */

/*---------------------------------------------------------------*/
/* Historical data   --------------------------------------------*/

	%if &histo ne %then %do; 
		data _temp_;
			set &histo;			
			/* Scenario_i_j contains the probability of true toxicity with J (t1, ...,tJ) row and I (s1,...,sI) column */
			array scenario(&max_level2,&max_level1);
			%let ind=0;
        	%do _j_=1 %to &max_level2;
				%do _i_=1 %to &max_level1;
			    	%let ind=%eval(&ind+1);
					scenario(&_j_,&_i_)=%scan(&scenario,&ind,"/");
				%end;  
			%end; 
			/* Working model ai and bj */
			array a(&max_level1);
			array b(&max_level2);

			%do _i_=1 %to &max_level1;
				a(&_i_)=%scan(&_a_,&_i_,"/");		
			%end;

			%do _j_=1 %to &max_level2;
				b(&_j_)=%scan(&_b_,&_j_,"/");
			%end;
		run;
	%end;
	
/*---------------------------------------------------------------*/
/* Start-up -----------------------------------------------------*/

	%if &start_up=0 %then %do;	
		data startup;
		    /* Scenario_i_j contains the probability of true toxicity with J (t1, ...,tJ) row and I (s1,...,sI) column */
			array scenario(&max_level2,&max_level1);
			%let ind=0;
        	%do _j_=1 %to &max_level2;
				%do _i_=1 %to &max_level1;
			    	%let ind=%eval(&ind+1);
					scenario(&_j_,&_i_)=%scan(&scenario,&ind,"/");
				%end;  
			%end; 
			/* Working model ai and bj */
			array a(&max_level1);
			array b(&max_level2);

			%do _i_=1 %to &max_level1;
				a(&_i_)=%scan(&_a_,&_i_,"/");		
			%end;

			%do _j_=1 %to &max_level2;
				b(&_j_)=%scan(&_b_,&_j_,"/");
			%end;

            /* Initialisation of toxicity probabiliy, patient number and toxicity number */
			%let size=%sysevalf(&max_level2*&max_level1);
			array n_pat(&max_level2,&max_level1);	retain n_pat1-n_pat&size 0;
			array n_tox(&max_level2,&max_level1);	retain n_tox1-n_tox&size 0;
			 
		run;		
			
		%let current_tj=1;%let current_si=1;%let nobs=0;

		%if &histo ne %then %do; data startup;
									set _temp_(in=a) startup(in=b);
									length group $7; 
									if a then group="D0";
									else if b then group="startup";																	 
                                 run;

								 proc datasets lib=work nolist;
								 	delete _temp_;
								 run;quit;
									
                      %end;
	%end;
/*---------------------------------------------------------------*/
	%else %if &start_up=1 %then %do;
		    
		data startup;			
			/* Scenario_i_j contains the probability of true toxicity with J (t1, ...,tJ) row and I (s1,...,sI) column */
			array scenario(&max_level2,&max_level1);
			%let ind=0;
        	%do _j_=1 %to &max_level2;
				%do _i_=1 %to &max_level1;
			    	%let ind=%eval(&ind+1);
					scenario(&_j_,&_i_)=%scan(&scenario,&ind,"/");
				%end;  
			%end; 
			/* Working model ai and bj */
			array a(&max_level1);
			array b(&max_level2);

			%do _i_=1 %to &max_level1;
				a(&_i_)=%scan(&_a_,&_i_,"/");		
			%end;
 
			%do _j_=1 %to &max_level2;
				b(&_j_)=%scan(&_b_,&_j_,"/");
			%end;
	
			/* Initialisation of toxicity probabiliy, patient number and toxicity number */
			%let size=%sysevalf(&max_level2*&max_level1);
			array n_pat(&max_level2,&max_level1);	retain n_pat1-n_pat&size 0;
			array n_tox(&max_level2,&max_level1);	retain n_tox1-n_tox&size 0;
			 
			/* Initialisation of tj agent 2 and si agent 1 */
			tj=1;si=1;

			/* Use to assign an idenfifiant for each observation and compute the number of toxicities for each cohort */
			retain patid sum 0;
			stop=0;

			/* Seed for the generation of bernoulli variable */
			call streaminit(&rep);
			do while (tj<=&max_level2 and stop=0);
				end_level=0;tox=0;
				do while (si<=&max_level1 and end_level=0 and tox=0);			
					do l=1 to &g_size1; /* Cohort of g_size1 patients */
						patid+1;
						p=scenario(tj,si);
	                    y=rand("bernoulli",p);						
						sum+y; /* Allow to compute the total number of toxicities */
						n_pat(tj,si)+1;
						n_tox(tj,si)+y;
						output;
					end;
					if sum ge 1 then do; /* We observe at least one toxicity. So the start-up is terminated at the current level of agent 2 */
								tox=1; 
								if si ge 3 then si=si-2;
								else stop=1; /* Stop because early toxicity but not end of level si */
					end;
					else if sum=0	then do;
							if si ne &max_level1 then si=si+1; /* Escalation of first agent i+1 */
							else if si=&max_level1 then do;end_level=1;si=si-2;end; /* End of level first agent with no toxicity */					
					end;		
					sum=0;	/* Useful to set at 0 the number of toxicities for the next cohort */
				end;
				tj+1;
			end;			
		run;	
		
		%if &histo ne %then %do; data startup;
									set _temp_(in=a) startup(in=b);
									length group $7; 
									if a then group="D0";
									else if b then group="startup";									
								 run;								 
							%end;		
/*---------------------------------------------------------------*/
/* Estimate the parameters and the probability of toxicity on data from the start-up */
		%ESTIMATE(in=startup);
		%DEFINE_DOSE(in=startup,mode=startup);

		proc sql noprint;	
				select count(*) into: nobs
				from startup
				%if &histo ne %then %do;
					where group ne "D0" /* We did not consider historical data in the sample size of two dimensional design */
				%end;;
		quit;
	%end;
/*---------------------------------------------------------------*/
/* End of start-up--> 1: acquire some data and dose level chosen for the starting combination for the 2 dimensional design */

	data simul&rep;
		set startup;
	run;
/*---------------------------------------------------------------*/
/* dose escalation ----------------------------------------------*/;
	
    %let  group_size__=&group_size;
    %do %while (&nobs lt &n_patient);
        %if %eval(&nobs+&group_size__) >&n_patient %then %do;%let group_size__=%eval(&n_patient-&nobs);%end;
		
		data simul&rep;
			
			retain patid;
			array scenario(&max_level2,&max_level1);

            %let size=%sysevalf(&max_level2*&max_level1);
			array n_pat(&max_level2,&max_level1); retain n_pat1-n_pat&size;
			array n_tox(&max_level2,&max_level1); retain n_tox1-n_tox&size;

			set simul&rep end=eof;
            output;
			
			if eof then do; 
			array prob(&max_level2,&max_level1);array diff(&max_level2,&max_level1);array search(&max_level2,&max_level1);
		   	%do _j_=1 %to &max_level2;
				%do _i_=1 %to &max_level1;
					prob(&_j_,&_i_)=.;
					diff(&_j_,&_i_)=.;
					search(&_j_,&_i_)=.;
				%end;  
				indice&_j_=.;
			%end;

				call streaminit(%eval(&rep+&nobs));				
				do l=1 to &group_size__; /* Cohort of n_design patients for 2 dimensional design */
							%if &histo ne %then %do;group="Current";%end;
							patid+1;
							tj=&current_tj;si=&current_si;
							p=scenario(tj,si);
	                        y=rand("bernoulli",p);						
							n_pat(tj,si)+1;
							n_tox(tj,si)+y;
							output;
				end;
			end;
		run;				

		%if  &start_up=0 %then %do;
		data simul&rep;
			set simul&rep;
			if patid=. then delete;
		run;
		%end;
	
		/* For the final analysis, we recommend one level of agent 1 by level of agent 2 */
		%let nobs=%eval(&nobs+&group_size__); 

		%if  (&nobs lt &n_patient) %then %do; 
			%ESTIMATE(in=simul&rep);
			%DEFINE_DOSE(in=simul&rep,mode=2dim_design);
		%end;

		%else %do;
			%ESTIMATE(in=simul&rep);
			%DEFINE_DOSE(in=simul&rep,mode=recommend_dose);
		%end;

	%end;
	
	data simul&rep;
		set simul&rep end=eof;

		array recommended_dose(&max_level2,&max_level1);

		if eof then do; 
		%do j=1 %to &max_level2;
			%do i=1 %to &max_level1;
				recommended_dose(&j,&i)=0;
			%end;
		%end;
		end;
		if eof then do;
			%do j=1 %to &max_level2;
				recommended_dose(&j,&&recom_dose&j)=1;
			%end;
		end;
	run;

data pool;
	set pool simul&rep(in=a);
	if a then simul=&rep;
run;

proc datasets lib=work noprint;
	delete startup postSumInt simul&rep;
run;
quit;

ods listing;

dm 'ODSRESULTS' clear editor;
%end; 	

%let end_time = %sysfunc(datetime());%put End Time: &end_time;
%let duration = %sysevalf((&end_time.-&init_time.));
%put Execution time : &duration. seconds;

/*---------------------------------------------------------------*/
/* End of simulation --------------------------------------------*/

/*---------------------------------------------------------------*/
/* Summary of results -------------------------------------------*/

proc sort data=pool;by simul;run;

data summary;
	set pool(keep=simul n_pat: n_tox: recommended_dose:);
	by simul;
	if last.simul;
run;

%REPORT(in=summary,var=n_pat);%REPORT(in=summary,var=n_tox);%REPORT(in=summary,var=recommended_dose);

data empty_row;
run;

data results_final;
	set mean_n_pat empty_row mean_n_tox empty_row mean_recommended_dose;
run;

footnote;
proc print data=results_final noobs;
run;

proc datasets lib=work noprint;
	delete pool summary empty_row mean mean_n_pat mean_n_tox mean_recommended_dose;
run;
quit;

%end;

%mend CRM2DIM;

