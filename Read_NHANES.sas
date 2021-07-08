options msglevel=n;

%macro CreateDS(myDS);
  %let i = 1;
  %let DS = %scan(&myDS, &i);

  %do %until(&DS = %nrstr());
    %let Prefix = %lowcase(%substr(&DS, 1, 2));
    %let Suffix = %lowcase(%substr(&DS, %eval(%length(&DS)-1)));

    %if (&Prefix = p_) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2017-2018/&DS..xpt"; %end;
    %else %if (&Suffix = _j) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2017-2018/&DS..xpt"; %end;
    %else %if (&Suffix = _i) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2015-2016/&DS..xpt"; %end;
    %else %if (&Suffix = _h) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2013-2014/&DS..xpt"; %end;
    %else %if (&Suffix = _g) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2011-2012/&DS..xpt"; %end;
    %else %if (&Suffix = _f) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2009-2010/&DS..xpt"; %end;
    %else %if (&Suffix = _e) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2007-2008/&DS..xpt"; %end;
    %else %if (&Suffix = _d) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2005-2006/&DS..xpt"; %end;
    %else %if (&Suffix = _c) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2003-2004/&DS..xpt"; %end;
    %else %if (&Suffix = _b) %then %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/2001-2002/&DS..xpt"; %end;
    %else %do; filename &DS url "https://wwwn.cdc.gov/nchs/nhanes/1999-2000/&DS..xpt"; %end;

    libname &DS xport;
    data &DS;
      set &DS..&DS;
    run;

    %let i = %eval(&i+1);
    %let DS = %scan(&myDS, &i);
  %end;
%mend CreateDS;
