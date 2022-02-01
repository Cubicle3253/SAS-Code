* Calculate Korn and Graubard confidence intervals from SUDAAN Proc Descript output           *;
* and apply National Center for Health Statistics Data Presentation Standards for Proportions *;
* For use with National Health and Nutrition Examination Survey (NHANES) data                 *;
* Programmer: Crescent Martin, NCHS                                                           *;
* v1.0 7/11/2019 for release to DHANES Analysis Branch                                        *;
* v1.1 - v1.2.5 : additional releases for enhancements and bug fixes                          *;
* v1.3   updated 1/8/2020 - was shared across NCHS Divisions                                  *;
* v1.3.1 updated 1/29/2021 - no change to calculations, changes only to clean up comments     *;
* v2.0 updated 3/10/2021 - added option df_adjust to support analysis of NHIS data            *;
*                        - added support for special case of 0 obs in an age-adjustment cell  *;
*                        - added option to output additional data presentation flags in same  *;
*                          format as HUS macro (for backward compatability with existing code)*;

%* Write macro version to the log (when code is '%include'-ed) for traceability *;
%put NOTE: Using macro calculate_KornGraubard_CIs version 2.0;

%macro calculate_KornGraubard_CIs(dsIn = &syslast. , 
                                  dsOut= &dsIn._KG ,                                   
                                  proportionTotal=100,
                                  roundingIncrement=0.1,
                                  df_adjust=yes,
                                  age_adjusted=no,
                                  ageVar=ageCat,
                                  ageWtShares=&ageShares,
                                  dropAgeSpecific=yes,
                                  HUS_output=no
                                  );
  %* ----------------------------------------------------------------------------------------------------------------------------------------------*;
  %* Macro to calculate Clopper-Pearson (exact) confidence limits for proportions, modified for complex surveys by the method of Korn and Graubard *;
  %* accounts for the degrees of freedom for each subpopulation (domain level)                                                                     *;
  %* i.e. DF = # of PSUs represented IN THE SUBPOPULATION minus # of strata represented IN THE SUBPOPULATION                                       *;
  
  %* Call this macro after a SUDAAN Proc Descript statement.                                                                                       *;
  %* ----------------------------------------------------------------------------------------------------------------------------------------------*;
  %* Macro parameters: All are optional (default values are provided, and user can override if desired)                                            *;
  %* dsIn =  Input dataset. Defaults to using the last SAS dataset created (using automatic macro variable syslast.)                               *;
  %*         Must contain values for the percent or the mean (of a 0/1 or 0/100 indicator), standard error, number of observations                 *;  
  %*         For NHANES data, must also contain atlev1 (nStrata per subdomain), atlev2 (nPSU per subdomain)                                        *;
  %* dsOut = Output dataset name. Defaults to appending _KG to the input dataset name.                                                             *;
  %* proportionTotal= Value meaning '100%'. Allowable values: 100 (default) or 1. Must match the coding on your indicator variable(s).             *;
  %*                  Set to 100 for 0/100 indicator variable (proportion shown as a percent e.g. 30) or 1 for 0/1 indicator (decimals e.g. 0.3)   *;
  %* roundingIncrement = rounding unit for Percent, SE, and confidence intervals, as applied to a Percent. Defaults to 0.1                         *;
  %*                     i.e. a calculated percent 12.345689 is rounded to 12.3 if proportionTotal=100 or .123 if proportionTotal=1                *;
  %* df_adjust = option to apply a degrees of freedom adjustment to the effective sample size. Yes (alias NHANES, is default) or No (alias NHIS)   *;
  %*             NHANES and NHIS allowable values reflect the analytic guidelines for each survey (as of February 2021)                            *;
  %*             If Yes or NHANES, the input dataset MUST contain variables atlev1 (nStrata per subdomain) and atlev2 (nPSU per subdomain)         *;
  %* age_adjusted = yes/no, whether to calculate confidence intervals for an age-adjusted proportion                                               *;
  %*                if yes - dataset must include age-specific estimates (for each age group used in age adjustment)                               *;
  %*                         as well as the age-adjusted estimates                                                                                 *;
  %* --- Parameters that are relevant only if age_adjusted=yes: ---------------------------------------------------------------------------------- *;
  %* ageVar = Name of the variable that contains the age groups for age adjustment (must match the STDVAR statement in the Proc Descript)          *;
  %*          Defaults to variable ageCat                                                                                                          *;
  %* ageWtShares = the age-adjustment weights containing the proportion of the standardization population in each of the age groups.               *;
  %*               Must match the STDWGT statement in the Proc Descript call                                                                       *;
  %*               Weights must add up to 1                                                                                                        *;
  %*               Defaults to macro variable ageShares                                                                                            *;
  %* dropAgeSpecific = yes/no, whether to drop the age-specific estimates from the output file with KG confidence intervals                        *;
  %* ----------------------------------------------------------------------------------------------------------------------------------------------*;
  %* -- For internal NCHS use --                                                                                                                   *;
  %* HUS_output = option to create additional output variables that map to the output from the Health, United States KG macro                      *;
  %*              kg_l, kg_u kg_l, kg_wdth, p_reliable, q_reliable, df_clerical, p_clerical                                                        *;
  
  
  %* Save users existing settings for some SAS options - so can reinstate them at the end of the macro *;
  %local sasOptionSettings ;
  %let sasOptionSettings = %sysfunc(getoption(MPRINT));
  
  %* Turn on mprint to ensure this data step code is printed to log *;
  option mprint;

  %* Save resolved names of input/output datasets - because syslast macro variable will get reset  *;
  %let dsIn_tmp = &dsIn;
  %let dsOut_tmp = %sysfunc(compress(&dsOut));
  %let dsIn=&dsIn_tmp;
  %let dsOut = &dsOut_tmp;

  %************************************************************;
  %* Check that required variables exist on the input dataset *;
  %************************************************************;
  %let varcheck_flg=0;

  %* open the input dataset and determine if variables MEAN and PERCENT are present *;
  %let dsid = %sysfunc(open(&dsIn));
  %let val_mean = %sysfunc(varnum(&dsid,mean));
  %let val_percent = %sysfunc(varnum(&dsid,percent));

  %* set the variable names to match the variables in the provided dataset - MEAN or PERCENT *;
  %* then check if the appropropriate SE variable is present *;
  %if %eval(&val_mean >0) %then %do;
    %let varname_mean=mean;
    %let varname_stderr=semean;

    %if %eval(%sysfunc(varnum(&dsid,semean)) = 0) %then %do;
      %put %SYSFUNC(COMPRESS(ER ROR:)) Input dataset &dsIn contains MEAN but not SEMEAN. ;
      %let varcheck_flg=1;
    %end;
  %end;
  %else %if %eval(&val_percent >0) %then %do;
    %let varname_mean=percent;
    %let varname_stderr=sepercent;

    %if %eval(%sysfunc(varnum(&dsid,sepercent)) = 0) %then %do;
      %put %SYSFUNC(COMPRESS(ER ROR:)) Input dataset &dsIn contains PERCENT but not SEPERCENT. ;
      %let varcheck_flg=1;
    %end;
  %end;
  %else %do;
      %put %SYSFUNC(COMPRESS(ER ROR:)) Input dataset &dsIn does not contain a variable named MEAN or PERCENT. ;
      %let varcheck_flg=1;
  %end;

  %* Check for NSUM, ATLEV1, ATLEV2 *;
  %let varname_n=nsum;
  %if %eval(%sysfunc(varnum(&dsid,nsum)) = 0) %then %do;
      %put %SYSFUNC(COMPRESS(ER ROR:)) Input dataset &dsIn does not contain NSUM. ;
      %let varcheck_flg=1;
  %end; 

  %* check allowable values for option df_adjust and check for ATLEV1, ATLEV2 variables on file if df_adjust is yes *;
  %if %upcase(%substr(&df_adjust,1,1)) = Y or %upcase(&df_adjust.) = NHANES %then %do;
    %if %eval(%sysfunc(varnum(&dsid,atlev1)) = 0) %then %do;
        %put %SYSFUNC(COMPRESS(ER ROR:)) Input dataset &dsIn does not contain ATLEV1 (number of strata.) ;
        %put %SYSFUNC(COMPRESS(ER ROR-)) ATLEV1 is required when option df_adjust is YES or NHANES. ;
        %let varcheck_flg=1;
    %end; 
    %if %eval(%sysfunc(varnum(&dsid,atlev2)) = 0) %then %do;
        %put %SYSFUNC(COMPRESS(ER ROR:)) Input dataset &dsIn does not contain ATLEV2 (number of PSUs.) ;
        %put %SYSFUNC(COMPRESS(ER ROR-)) ATLEV2 is required when option df_adjust is YES or NHANES. ;
        %let varcheck_flg=1;
    %end;
  %end;
  %else %if %upcase(&df_adjust.) = NO or %upcase(&df_adjust.) = NHIS %then %do;
    %* no further checks required *;
    %* standardize allowable value of df_adjust for checks later in the program ;
    %let df_adjust=NO;
  %end;
  %else %do; %* invalid value provided for df_adjust option *;
    %put %SYSFUNC(COMPRESS(ER ROR:)) Invalid macro option df_adjust (&df_adjust.); 
    %put %SYSFUNC(COMPRESS(ER ROR-)) Allowable values are YES (alias NHANES) or NO (alias NHIS) ;
    %let varcheck_flg=1;
  %end;

  %* close dataset *;
  %let rc = %sysfunc(close(&dsid));

  %* exit macro if problems were found *;
  %if %eval(&varcheck_flg=1) %then %return;

  %************************************************************;
  %* Check that the proportionTotal option is consistent with data values provided (decimal vs. percent)*;  
  %let varcheck_decimal=0;
  proc means data = &dsIN noprint;
    var &varname_mean;
    output out=_checkRange max=max;
  run;
  
  data _null_;
    set _checkRange;
    %if &proportionTotal=100 %then %do;
      if max<=1 then do;
        putlog "WAR" "NING: Macro option proportionTotal=100 indicating the outcomes are specified as percents (0/100 coding).";
        putlog "WAR" "NING- But all values of &varname_mean are <=1. Your outcomes may be specified as decimals (0/1 coding) instead.";
        putlog "WAR" "NING- Double-check your proportion values, and if they are decimals then set macro option proportionTotal=1.";  
        putlog "WAR" "NING- If you are analyzing an outcome with very low prevalence (<1 percent), you may ignore this message.";
      end;
    %end;
    %else %if &proportionTotal=1 %then %do;
      if max>1 then do;
        putlog "ER" "ROR: Macro option proportionTotal=1 indicating the outcomes are specified as decimals (0/1 coding).";
        putlog "ER" "ROR- But some values of &varname_mean are >1, an invalid proportion.";
        putlog "ER" "ROR- Double-check your proportion values, and if they are percents then set macro option proportionTotal=100.";         
        call symput("varcheck_decimal", 1);
      end;
    %end;
  run;  

  proc datasets library=work nolist ;
    delete _checkRange;
  quit;

  %* exit macro if problems were found *;
  %if %eval(&varcheck_decimal=1) %then %return;  

  %*************************************************************************;
  %** For age-adjusted estimates: calculate denominator of design effect  **;
  %*************************************************************************;

  %if %upcase(%substr(&age_adjusted,1,1)) = Y %then %do;
    ** Prepare design effect for age-adjusted estimates **;
    %let calc_age_adjusted=Y;

    %* Confirm that the age variable exists on the input dataset (i.e. there are age-specific estimates), else exit *;
    %local rc dsid ;
    %let dsid = %sysfunc(open(&dsIN));
 
    %if %eval(%sysfunc(varnum(&dsid, &ageVar)) = 0) %then %do;
      %let rc = %sysfunc(close(&dsid));
      %put %SYSFUNC(COMPRESS(ER ROR:)) Error in calculating age-adjusted confidence intervals. ;
      %put %SYSFUNC(COMPRESS(ER ROR-)) Specified ageVar &ageVar. does not exist on input dataset &dsIn.. ;
      %put %SYSFUNC(COMPRESS(ER ROR-)) You must also provide age-specific estimates (needed for the calculation of the design effect.) ;
      %return;
    %end;
    %else %let rc = %sysfunc(close(&dsid));

    * Prep files of age-adjusted and non-age-adjusted estimates *;
    data _ageAdjusted 
         _ageSpecific;
      set &dsIN;
      origSort=_n_;
      %* Values indicating age-adjusted estimates: *;
      %* Coded value 0 in SUDAAN output indicates the corresponding value of the class variable is a TOTAL. *;
      %* Coded value -2 in SUDAAN output indicates the corresponding value of the class variable is not used in that table. *;
      %* Missing value from SAS output indicates the corresponding value of the class variable is not used in that domain. *;
      if &ageVar in (., 0, -2 ) then output _ageAdjusted ;
      else output _ageSpecific;      
    run;     

    * Put list of age groups into a macro variable *;
    proc sort nodupkey data = _ageSpecific out=_ages;
      by &ageVar;
    run;
    
    %* remove any formats that may be on ageVar (i.e. a label) *;
    data _ages;
      set _ages;
      format &ageVar.;
    run;

    proc sql noprint;
      select  &ageVar. into: ageList separated by " "
      from _ages;
    quit;  

    %* count number of unique age groups and number of age weights provided *;
    %let nAgeShares = %sysfunc(countw("&ageWtShares",,s));
    %let nAges = %sysfunc(countw("&ageList",,s));

    %if (&nAgeShares. ne &nAges.) %then %do;
      %put %SYSFUNC(COMPRESS(ER ROR:)) Number of age-adjustment weights does not match number of age groups in the input dataset.;
      %put %SYSFUNC(COMPRESS(ER ROR-)) Number of age-adjustment weights: &nAgeShares.. ageWtShare macro variable = &ageWtShares;
      %put %SYSFUNC(COMPRESS(ER ROR-)) Number of age groups on input dataset: &nAges.. Levels of variable &ageVar.: &ageList;
      %return;
    %end;

    %* check the sum of age weights provided *;
    %let WtSumCheck=0; 
    %do i = 1 %to %sysfunc(countw("&ageWtShares.", ' ') );                                                                                                              
        %let WtSumCheck= %sysevalf(&WtSumCheck. + %scan(&ageWtShares, &i, ' '));                                                                                                                                                                         
    %end;                

    %if (%sysfunc(fuzz(&WtSumCheck.)) ne 1 ) %then %do;
      %put %SYSFUNC(COMPRESS(ER ROR:)) Age-adjustment weights do not sum to 1. Sum=&WtSumCheck.;
      %put %SYSFUNC(COMPRESS(ER ROR-)) Correct macro variable ageWtShares and resubmit.;
      %put %SYSFUNC(COMPRESS(ER ROR-)) ageWtShares was set to: &ageWtShares.;
      
      %return;
    %end; 

    ** Calculate denominator of design effect (variance of simple random sample) **;
    data _calc_denom_deff;
      set _ageSpecific end=eof ;
      varsrs= &varname_mean. * (&proportionTotal. - &varname_mean.) / &varname_n.;
      * specify the age-adjustment weight for each age group *;
      select(&ageVar);
      %do i = 1 %to &nAges;
        when (%sysfunc(scan("&ageList", &i, " "))) ageWt= %sysfunc(scan("&ageWtShares", &i, " "));
      %end;
      otherwise ageWt=0;
      end; 
      prop=varsrs*(ageWt**2);      
      * check for any standardization cells with 0 observations *;
      retain reweightFlag;
      if nsum=0 then reweightFlag=1;
      if eof then call symputx("reweightFlag", reweightFlag);
      drop reweightFlag;
      reweighted=0;
    run;

    %* Get list of analysis ("table") variables *;
    proc contents noprint data = &dsIN out= _cont;
    run;

    data _analysisVars ;
      set _cont;      
      %* delete variables that are standard output of proc descript *;
      %* except: keep "variable" (in case user analyzes more than one variable or more than one catlevel in same proc descript) *;
      if upcase(name) in ( %upcase("&ageVar"), "ATLEV1", "ATLEV2",  %upcase("&varname_mean."), 
                           %upcase("&varname_stderr."), %upcase("&varname_n."), "PROCNUM", "TABLENO" ) then delete;
      if upcase(name) in ("SMLCELL", "NSUM", "WSUM", "TOTAL", "SETOTAL", "LOWTOTAL", "UPTOTAL", "MEAN", "SEMEAN", "LOWMEAN", "UPMEAN", 
                          "T_MEAN", "P_MEAN", "PERCENT", "SEPERCENT", "LOWPCT", "UPPCT", "T_PCT", "P_PCT", "GEOMEAN", "SEGEOMEAN", "DEFFMEAN",
                          "DEFFTOTAL", "DEFFPCT", "ATLEV1", "ATLEV2", "VARMEAN", "VARPCT", "VARTOTAL", "DDFMEAN", "DDFPCT", "DDFTOTAL") then delete;
      if name=:"_C" then delete;
    run;

    proc sql noprint;
      select name into: analysisVarList separated by " "
      from _analysisVars;
    quit;

    %* If any standardization cell within a subgroup has 0 observations  - renormalize weights over remaining cells for that subgroup *;
    %* to match SUDAAN Proc Descript behavior (v10 and v11) *;
    %if &reweightFlag.=1 %then %do;

      * get sum of standardization weights with observations *;
      proc means noprint data = _calc_denom_deff (where=(nsum>0)) nway missing;
        var ageWt;
        class &analysisVarList;
        output out= _wtSums sum(ageWt) = sum_ageWts;
      run;    

      * merge back to age-specific estimates *;
      proc sort data = _calc_denom_deff ;
        by &analysisVarList.;
      run;

      data _calc_denom_deff;
        merge _calc_denom_deff (in = a)
              _wtSums (keep = &analysisVarList. sum_ageWts);
        by &analysisVarList.;
        * renormalize age shares *;
        if nsum>0 then do;
          ageWt=ageWt/sum_ageWts;
          prop=varsrs*(ageWt**2);
        end;
        else if nsum=0 then ageWt=0;
        * flag estimates with weights that are renormalized  *;
        reweighted=(fuzz(sum_ageWts)<1);
      run;

    %end;

    %* Sum across the age groups to get the SRS variance *;
    %* check sums: age weights should add up to 1 and weighted average of age-specific means should equal age-adusted mean *;
    proc means noprint data = _calc_denom_deff (where=(nsum>0)) nway missing;
      var prop ageWt reweighted;
      var &varname_mean. / weight=ageWt;
      class &analysisVarList;
      output out= _denom_deff sum(prop)=srsvar sum(ageWt) = check_ageShares mean(&varname_mean.)=check_ageAdjMean max(reweighted)=reweighted;
    run;

    * Merge to age-adjusted estimates;
    proc sort data = _ageAdjusted ;
      by &analysisVarList;
    run;
    
    data _ageAdjusted _aaProblems _aaReweight;
      merge _ageAdjusted (in = a)
            _denom_deff (in = b) end=eof;
      by &analysisVarList ;
      if a;
      deffmean=&varname_stderr.**2/srsvar;
      * check for problems with dataset merge or the calculation of the denominator of design effect *;
      * and flag if standardization weights were renormalized *;
      problemFlag=0;
      if a and not b then problemFlag=1;
      else if abs(check_ageAdjMean - &varname_mean) > .001 then problemFlag=2;      
      if problemFlag ne 0 then output _aaProblems;
      if reweighted=1 then output _aaReweight;
      output _ageAdjusted;
    run;

    %* print to log if problems *; 
    %* check if any records were flagged *;
    proc sql noprint;
      select nobs into : nProblems
      from dictionary.tables
      where libname='WORK' and memname='_AAPROBLEMS';
    quit;

    %if (&nProblems > 0  or &reweightFlag.=1) %then %do;

      %* Increment highest title number so do not clobber users existing title statements *;    
      data _savedTitles;
        number=0;
        output;
      run;

      data _savedTitles;
        set _savedTitles
            sashelp.vtitle (where = (type="T"));
        call symputx("newTitleNum", number+1);
      run;

      %if &reweightFlag.=1 %then %do;

        title&newTitleNum. "Estimates where one or more age standardization cells has 0 observations";

        proc print data = _aaReweight noobs;
          var &analysisVarList.;
        run;
 
        title&newTitleNum.;

      %end;

      %if &nProblems > 0 %then %do;
        title&newTitleNum. "Problems in calculation of the denominator of the design effect.";
        proc format;
          value de_msg
          1="Corresponding age-specific estimates not found on input dataset."
          2="Option ageWtShares does not match the age weights used for the estimate."
          ;
        run;  
              
        proc print data = _ageAdjusted (where = (problemFlag>0));
          format problemflag de_msg.;
        run;      

        title&newTitleNum.;

        %put %SYSFUNC(COMPRESS(ER ROR:)) Problems in calculation of the denominator of the design effect for some or all estimates.;
        %put %SYSFUNC(COMPRESS(ER ROR-)) Check output / list file or dataset _aaProblems for records with issues. ;

        %return;
      %end;
    %end;

    %* Reassemble age-adjusted and age-specific estimates, OR just use age-adjusted if option dropAgeSpecific= yes *;
    data _combined;
      %if %upcase(%substr(&dropAgeSpecific,1,1)) = Y %then %do;
      set _ageAdjusted (in = a drop = problemFlag check_ageShares check_ageAdjMean _type_ _freq_ srsvar);
      %end;
      %else %do;
      set _ageAdjusted (in = a drop = problemFlag check_ageShares check_ageAdjMean _type_ _freq_ srsvar)
          _ageSpecific (in = b);
      %end;
      if a then ageAdjusted=1;
      else ageAdjusted=0;
    run;

    proc sort data = _combined;
      by origSort;
    run;

* comment out clean-up for testing *;
    %* Clean up intermediate datasets *;
    proc datasets library=work nolist ;
      delete _ageAdjusted _ageSpecific _ages _calc_denom_deff  _cont _analysisVars _denom_deff _aaProblems _savedTitles;
    quit;

  %end;
  %else %do;
    %let calc_age_adjusted=N;
  %end;

  %*************************************************************************;
  %** Calculate Korn and Graubard confidence intervals                    **;
  %*************************************************************************;
  data &dsOut;
    %if &calc_age_adjusted. = Y %then %do; set _combined; %end;
    %else %do; set &dsIn ; %end;

    * save original variable, if input dataset contains percent (i.e. SUDAAN Proc Descript with catlevel statement) *;
    %if &varname_mean = percent %then %do;
      percent_orig=&varname_mean.;
      &varname_stderr._orig= &varname_stderr.;
    %end;

    * create decimal values for calculations *;
    p= fuzz(&varname_mean / &proportionTotal); * fuzz statement sets to integer (0 or 1) if the value is within 1E-12 of an integer *;
    sep= &varname_stderr / &proportionTotal;
    * for display as percent *;
    percent=round(p*100,&roundingIncrement);    
    sepercent=round(sep*100, &roundingIncrement);    

    %if &calc_age_adjusted ne Y %then %do;
      * specify that all calculations are crude (not age-adjusted) *;
      ageAdjusted=0;
    %end;
  
    if p in (0,1) then do;    
      * Special treatment for confidence interval limits if estimated proportion = 0 or 1 *;
      n_eff=&varname_n.;
      n_eff_df= &varname_n.;
      x=p*n_eff_df;

      if p=0 then do;
        lowerCL=0;
        upperCL = (1 + ((n_eff_df - x    ) / ((x+1) *finv(0.975, 2*(x+1), 2*(n_eff_df-x  )))))**(-1);
      end;
      else if p=1 then do;
        lowerCL = (1 + ((n_eff_df - x + 1) / ( x    *finv(0.025, 2*x    , 2*(n_eff_df-x+1)))))**(-1);
        upperCL = 1;
      end;
    end;
    else do;  
      * calculate modified Clopper-Pearson confidence intervals for a proportion, according to the approach of Korn and Graubard  *;
      * n_eff = effective sample size *;
      if ageAdjusted=1 then n_eff = &varname_n./deffmean;      
      else do;
        n_eff=p*(1-p)/sep**2;
        deffmean= (sep**2)/ (p*(1-p)/nsum);
      end;

      * Adjustment for degrees of freedom (DF);
      %if &df_adjust=NO %then %do;
        * Option df_adjust=NO or NHIS was specified, so do not apply degrees-of-freedom adjustment to effective sample size. *;
        t_adj=1;        
      %end;
      %else %do;
        * Option df_adjust=Yes or NHANES was specified, so apply a degrees-of-freedom adjustment to effective sample size. ;
        * DF = degrees of freedom = # of clusters - # of strata *;
        df=atlev2-atlev1;
        t_adj=(tinv(.975,&varname_n. -1)/tinv(.975,df))**2;
      %end;

      * Effective sample size adjusted for DF, capped at actual sample size (replaces sample size in the Clopper-Pearson computation);
      n_eff_df=min(&varname_n. ,n_eff*t_adj);

      * proportion x DF-adjusted effective sample size (replaces the number of positive responses in the Clopper-Pearson computation) *;    
      x=p*n_eff_df;
      %* use quantile function instead of finv (b/c the finv function breaks if DF < 1 in extreme cases) *;
      %*lowerCL = (1 + ((n_eff_df - x + 1) / ( x    *finv(0.025, 2*x    , 2*(n_eff_df-x+1)))))**(-1);
      %*upperCL = (1 + ((n_eff_df - x    ) / ((x+1) *finv(0.975, 2*(x+1), 2*(n_eff_df-x  )))))**(-1);
      lowerCL = (1 + ((n_eff_df - x + 1) / ( x    *quantile('F',0.025, 2*x    , 2*(n_eff_df-x+1)))))**(-1);
      upperCL = (1 + ((n_eff_df - x    ) / ((x+1) *quantile('F',0.975, 2*(x+1), 2*(n_eff_df-x  )))))**(-1);

    end;
    
    * CIW=confidence interval width (in decimal) *;
    CIW=upperCL-lowerCL;
    
    * RCIW=relative confidence interval width *;
    if p>0 then RCIW=round((CIW/p)*100, &roundingIncrement);
    if p<1 then RCIW_complement = round((CIW/(1-p))*100, &roundingIncrement);
    
    %if %upcase(%substr(&HUS_output,1,1)) = Y %then %do;
      * create output variables to match Health, United States macro output *;
      %* always output as decimals with full computer precision *;
      kg_l=lowerCL;
      kg_u=upperCL;
      kg_wdth=CIW;
    %end;  
    
    %* Convert CIW, lower and upper bounds to percent *;
    CIW = round(CIW*100, &roundingIncrement);
    lowerCL=round(lowerCL*100, &roundingIncrement);
    upperCL=round(upperCL*100, &roundingIncrement);
        
    %if %eval(&proportionTotal.=1) %then %do;
      ** convert results back to decimal, if proportions were provided as decimals **;
      array toDecimal (*) percent sepercent lowerCL upperCL;
      do i = 1 to dim(toDecimal);
        toDecimal{i} = toDecimal{i}/100;
      end;
      drop i;
    %end;

    ** Add formatted columns for presentation in data table **;
    %* number of digits post decimal (when displayed as a percent) *;
    postDecimal=log10(1/&roundingIncrement);
    %if &proportionTotal=100 %then %do;
      fmt = compress(put(postDecimal+4, best8.) || "." || put(postDecimal, best8.));
    %end;
    %else %if &proportionTotal=1 %then %do;
      fmt = compress(put(postDecimal+4, best8.) || "." || put(postDecimal+2, best8.));
    %end;      
    drop postDecimal fmt;
    
    length pct_se $50. pct_ci $50.;
    * percent and SE *;
    pct_se=putn(percent, fmt) ||  " (" || strip(putn( sepercent , fmt)) || ")";    
    label pct_se="Percent (std err)";
    * percent and CI *;
    pct_ci=putn(percent, fmt) ||  " (" || strip(putn( lowerCL , fmt)) || ", " || strip(putn( upperCL ,  fmt))  || ")";    
    label pct_ci="Percent (K-G CI)";

    ** apply criteria of NCHS Data Presentation Standards for Proportions **;
    suppress1=(min(n_eff, &varname_n.)<30 );
    suppress2=(CIW>=30);
    suppress3=( not(CIW<=5) and (RCIW>130));
    * if fuzzed proportion=0 or 1 --> flag for consideration as if number of events (or its compliment) = 0 ;    
    statReview1= (p in (0,1) );
    %if &df_adjust=NO %then %do;
      %* Assume large degrees of freedom *;
      statReview2=0;
    %end;
    %else %do;
      statReview2=(df<8);
    %end;

    * Result of data presentation standards checks: suppress, statistical review, or present *;
    suppress=(suppress1 or suppress2 or suppress3);
    %* Apply hierarchy of outcomes -- if estimate requires suppression then stat review is moot *;
    statReview=not suppress and (statReview1 or statReview2);
    present = (not suppress and not statReview);

    * for proportions that may be presented: check if need to include footnote that the complementary proportion (1-p) may be unreliable *;
    fn_complement=( present and CIW>5 and (RCIW_complement>130)) ;

    label suppress="Suppress"
          statReview="Statistical review"
          present="Present"
          suppress1="Suppress: effective sample size <=30"
          suppress2="Supress: CIW>=30"
          suppress3="Suppress: CIW>5 and RCIW>130"
          statReview1="Stat review: Events=0 or rounded Pct=0"
          statReview2="Stat review: DF<8"
          fn_complement="Complement may be unreliable"
          percent="Percent"
          sepercent="SE Percent"
          lowerCL="Lower CL"
          upperCL="Upper CL"
          n_eff="Effective sample size"
          n_eff_df="Effective sample size adjusted for DF"
          deffmean="Design effect"; 

    %if %upcase(%substr(&HUS_output,1,1)) = Y %then %do;
      * create output variables to match Health, United States macro output *;     
      p_reliable = (present or statReview);
      if p_reliable then q_reliable=1-fn_complement;
      df_clerical=statReview2;
      p_clerical=statReview;      
    %end;  

    drop p sep  ;
    drop t_adj x RCIW_complement;
    %if &calc_age_adjusted eq Y %then %do; drop origsort; %end;
  run;


  %if &calc_age_adjusted eq Y %then %do;  

    %* if any estimates are reweighted - write info to log at end of macro processing *;
    %if &reweightFlag.=1 %then %do;
      %put %SYSFUNC(COMPRESS(WA RNING:)) One or more standardization cells are 0 for one or more estimates.;
      %put %SYSFUNC(COMPRESS(WA RNING-)) Standardization weights have been renormalized over remaining cells;
      %put %SYSFUNC(COMPRESS(WA RNING-)) in accordance with SUDAAN Proc Descript calculations.;
      %put %SYSFUNC(COMPRESS(WA RNING-)) Check output / list file for affected estimates and consider the analytic implications.;
    %end; 

    * delete intermediate dataset *;
    proc datasets library=work nolist ;
        delete _combined;
    quit;
  %end;

  %* Restore users SAS options (mprint) *;
  options &sasOptionSettings;

%mend calculate_KornGraubard_CIs;
