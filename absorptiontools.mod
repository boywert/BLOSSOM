	  ã(  l   k820309    Å
          11.1        nÁQ                                                                                                           
       /home1/01937/cs390/BLOSSOM/src/absorptiontools.f90 ABSORPTIONTOOLS                                                    
                                                                             `                                                                                           
                                                        
                                                        
                 
                    `û!	@        3.14159265359                                                 
                                                        
                                                   	     
                                                   
     
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                    `                                                                                                             #MPIPRIV1%MPI_BOTTOM    #MPIPRIV1%MPI_IN_PLACE    #MPIPRIV1%MPI_STATUS_IGNORE                                                                                                                                                                                                      p          p            p                                                                                                     #MPIPRIV2%MPI_STATUSES_IGNORE    #MPIPRIV2%MPI_ERRCODES_IGNORE                                                                              p          p          p            p          p                                                                                                           p          p            p                                                                                                      #MPIPRIVC%MPI_ARGVS_NULL !   #MPIPRIVC%MPI_ARGV_NULL "   -                                         !                                 p          p          p            p          p                                  -                                         "                                p          p            p                                  #         @                                  #                  #N_CAL%TRIM $   #N %   #R &   #RHO '                                              $     TRIM           D                                 %                     D                                &                    
     p           & p         5 r          5 r    p         p                                   D                                '                    
     p           & p         5 r          5 r    p         p                          %         @                                (                   
       #FIND_T_K%SQRT )   #M0 *   #ZCOLL +                                              )     SQRT                                           *     
                                                 +     
       #         @                                  ,               	   #TAU_CAL%LOG10 -   #TAU_CAL%SQRT .   #M0 /   #ZCOLL 0   #IMPACT_PARAM 1   #N 2   #R 3   #RHO 4   #SIGMA_V 5   #AREATAU0 6   #TAU0 7                                                                    -     LOG10                                            .     SQRT                                           /     
                                                 0     
                                                 1     
                                                 2                                                    3                    
     p           & p         5 r          5 r    p         p                                                                   4                    
     p           & p         5 r          5 r    p         p                                    D                                5     
                 D                                6     
                 D                                7     
       #         @                                 8                  #LORENTZ_TERM_CAL%REAL 9   #N :   #F0 ;   #DELTA_F <   #GAMMA =   #LORENTZ_TERM >                               @              9     REAL                                            :                                                      ;     
                                                 <     
                                                 =     
                D                                >                    
      5  p        r :     &  5  p        r :   5  p        r :         5  p        r :    5  p        r :   p                          #         @                                 ?                  #DOPPLER_TERM_CAL%EXP @   #DOPPLER_TERM_CAL%SQRT A   #DOPPLER_TERM_CAL%REAL B   #N C   #F0 D   #DELTA_F E   #SIGMA F   #DOPPLER_TERM G                                              @     EXP                                            A     SQRT                             @              B     REAL                                            C                                                      D     
                                                 E     
                                                 F     
                D                                G                    
      5  p        r C     &  5  p        r C   5  p        r C         5  p        r C    5  p        r C   p                          #         @                                  H                  #VOIGT_TERM_CAL%REAL I   #N J   #F0 K   #DELTA_F L   #GAMMA M   #SIGMA N   #VOIGT_TERM O                               @              I     REAL           D @                               J                      D @                              K     
                 D @                              L     
                 D @                              M     
                 D @                              N     
                D                                O                    
      5  p        r J     &  5  p        r J   5  p        r J         5  p        r J    5  p        r J   p                          %         @                                P                   
       #VOIGT_FWHM%ALOG Q   #VOIGT_FWHM%SQRT R   #GAMMA S   #SIGMA T   #NU0 U                                              Q     ALOG                                            R     SQRT                                           S     
                                                 T     
                                                 U     
       (        `                                V                                   
    #ABSORPTION_LINE%CEILING W   #ABSORPTION_LINE%EXP X   #ABSORPTION_LINE%REAL Y   #NU_MIN Z   #NU_MAX [   #NU_BINSIZE \   #TOTALBIN ]   #NU_CENTER ^   #FWHM _   #TAU `   p          5 O p            5 O p                                                                     W     CEILING                                            X     EXP                             @              Y     REAL                                           Z     
                                                 [     
                                                 \     
                                                  ]                                                      ^     
                                                 _     
                                                 `     
       #         @                                  a                  #SAVE_ABSORPTIONLINE%LEN_TRIM b   #SAVE_ABSORPTIONLINE%ADJUSTL c   #SAVE_ABSORPTIONLINE%TRIM d   #SAVE_ABSORPTIONLINE%REAL e   #ABSORPLINE f   #N g   #NU_MIN h   #NU_BINSIZE i   #Z_S j   #FILENAME k                                              b     LEN_TRIM                                            c     ADJUSTL                                            d     TRIM                             @              e     REAL                                          f                    
     p          5  p        r g       5  p        r g                                                                g                                                      h     
                  @                              i     
                 D @                              j     d                                 D @                              k                                   K      fn#fn    ë   @   j   COMMON_VARS %   +  @       MAX_SIZE+COMMON_VARS -   k  @       DEN_PROFILE_FILE+COMMON_VARS #   «  @       ZETA_T+COMMON_VARS     ë  @       NU0+COMMON_VARS    +  }       PI+COMMON_VARS '   ¨  @       RHO_CRIT_0+COMMON_VARS $   è  @       LAMBDA0+COMMON_VARS $   (  @       OMEGA_0+COMMON_VARS #   h  @       ETASUS+COMMON_VARS #   ¨  @       ETATIS+COMMON_VARS    è  @       H+COMMON_VARS $   (  @       OMEGA_B+COMMON_VARS "   h  @       M_SOL+COMMON_VARS "   ¨  @       MTTIL+COMMON_VARS '   è  @       GRAV_CONST+COMMON_VARS "   (  @       AMASS+COMMON_VARS #   h  @       BOLTZK+COMMON_VARS    ¨  @       NH+COMMON_VARS    è  @       MU+COMMON_VARS    (  @       C+COMMON_VARS     h  @       A10+COMMON_VARS #   ¨  @       T_STAR+COMMON_VARS (   è  @       RESULT_PATH+COMMON_VARS >   (  ¤      MPI_CONSTANTS!MPIPRIV1+MPI_CONSTANTS=MPIPRIV1 2   Ì  H      MPIPRIV1%MPI_BOTTOM+MPI_CONSTANTS 4     H      MPIPRIV1%MPI_IN_PLACE+MPI_CONSTANTS 9   \  ¤      MPIPRIV1%MPI_STATUS_IGNORE+MPI_CONSTANTS >    	        MPI_CONSTANTS!MPIPRIV2+MPI_CONSTANTS=MPIPRIV2 ;   	  Ä      MPIPRIV2%MPI_STATUSES_IGNORE+MPI_CONSTANTS ;   X
  ¤      MPIPRIV2%MPI_ERRCODES_IGNORE+MPI_CONSTANTS >   ü
        MPI_CONSTANTS!MPIPRIVC+MPI_CONSTANTS=MPIPRIVC 6     Ä      MPIPRIVC%MPI_ARGVS_NULL+MPI_CONSTANTS 5   I  ¤      MPIPRIVC%MPI_ARGV_NULL+MPI_CONSTANTS    í  o       N_CAL    \  =      N_CAL%TRIM      @   a   N_CAL%N    Ù  Ä   a   N_CAL%R      Ä   a   N_CAL%RHO    a  v       FIND_T_K    ×  =      FIND_T_K%SQRT      @   a   FIND_T_K%M0    T  @   a   FIND_T_K%ZCOLL      ä       TAU_CAL    x  >      TAU_CAL%LOG10    ¶  =      TAU_CAL%SQRT    ó  @   a   TAU_CAL%M0    3  @   a   TAU_CAL%ZCOLL %   s  @   a   TAU_CAL%IMPACT_PARAM    ³  @   a   TAU_CAL%N    ó  Ä   a   TAU_CAL%R    ·  Ä   a   TAU_CAL%RHO     {  @   a   TAU_CAL%SIGMA_V !   »  @   a   TAU_CAL%AREATAU0    û  @   a   TAU_CAL%TAU0 !   ;         LORENTZ_TERM_CAL &   ×  =      LORENTZ_TERM_CAL%REAL #     @   a   LORENTZ_TERM_CAL%N $   T  @   a   LORENTZ_TERM_CAL%F0 )     @   a   LORENTZ_TERM_CAL%DELTA_F '   Ô  @   a   LORENTZ_TERM_CAL%GAMMA .        a   LORENTZ_TERM_CAL%LORENTZ_TERM !   4  Ñ       DOPPLER_TERM_CAL %     <      DOPPLER_TERM_CAL%EXP &   A  =      DOPPLER_TERM_CAL%SQRT &   ~  =      DOPPLER_TERM_CAL%REAL #   »  @   a   DOPPLER_TERM_CAL%N $   û  @   a   DOPPLER_TERM_CAL%F0 )   ;  @   a   DOPPLER_TERM_CAL%DELTA_F '   {  @   a   DOPPLER_TERM_CAL%SIGMA .   »     a   DOPPLER_TERM_CAL%DOPPLER_TERM    Û  £       VOIGT_TERM_CAL $   ~  =      VOIGT_TERM_CAL%REAL !   »  @   a   VOIGT_TERM_CAL%N "   û  @   a   VOIGT_TERM_CAL%F0 '   ;  @   a   VOIGT_TERM_CAL%DELTA_F %   {  @   a   VOIGT_TERM_CAL%GAMMA %   »  @   a   VOIGT_TERM_CAL%SIGMA *   û     a   VOIGT_TERM_CAL%VOIGT_TERM             VOIGT_FWHM     ´  =      VOIGT_FWHM%ALOG     ñ  =      VOIGT_FWHM%SQRT !   .   @   a   VOIGT_FWHM%GAMMA !   n   @   a   VOIGT_FWHM%SIGMA    ®   @   a   VOIGT_FWHM%NU0     î   \      ABSORPTION_LINE (   J"  @      ABSORPTION_LINE%CEILING $   "  <      ABSORPTION_LINE%EXP %   Æ"  =      ABSORPTION_LINE%REAL '   #  @   a   ABSORPTION_LINE%NU_MIN '   C#  @   a   ABSORPTION_LINE%NU_MAX +   #  @   a   ABSORPTION_LINE%NU_BINSIZE )   Ã#  @   a   ABSORPTION_LINE%TOTALBIN *   $  @   a   ABSORPTION_LINE%NU_CENTER %   C$  @   a   ABSORPTION_LINE%FWHM $   $  @   a   ABSORPTION_LINE%TAU $   Ã$        SAVE_ABSORPTIONLINE -   Ô%  A      SAVE_ABSORPTIONLINE%LEN_TRIM ,   &  @      SAVE_ABSORPTIONLINE%ADJUSTL )   U&  =      SAVE_ABSORPTIONLINE%TRIM )   &  =      SAVE_ABSORPTIONLINE%REAL /   Ï&  ´   a   SAVE_ABSORPTIONLINE%ABSORPLINE &   '  @   a   SAVE_ABSORPTIONLINE%N +   Ã'  @   a   SAVE_ABSORPTIONLINE%NU_MIN /   (  @   a   SAVE_ABSORPTIONLINE%NU_BINSIZE (   C(  P   a   SAVE_ABSORPTIONLINE%Z_S -   (  P   a   SAVE_ABSORPTIONLINE%FILENAME 