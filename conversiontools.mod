	  n  .   k820309    Å
          11.1        nÁQ                                                                                                           
       /home1/01937/cs390/BLOSSOM/src/conversiontools.f90 CONVERSIONTOOLS                                                    
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                        
                                                   	     
                                                   
     
                                                        
                                                        
                    @                                   	                                                        
                                                                          #MPIPRIV1%MPI_BOTTOM    #MPIPRIV1%MPI_IN_PLACE    #MPIPRIV1%MPI_STATUS_IGNORE                                                                                                                                                                                                      p          p            p                                                                                                     #MPIPRIV2%MPI_STATUSES_IGNORE    #MPIPRIV2%MPI_ERRCODES_IGNORE                                                                              p          p          p            p          p                                                                                                           p          p            p                                                                                                     #MPIPRIVC%MPI_ARGVS_NULL    #MPIPRIVC%MPI_ARGV_NULL    -                                                                          p          p          p            p          p                                  -                                                                         p          p            p                                  %         @                                                    
       #Z                                                   
       %         @                                                   
       #D_TO_Z%SQRT    #D                                                    SQRT                                                
       %         @                                                   
       #D_TO_NU%SQRT    #D                                                     SQRT                                                 
       %         @                                !                    
       #NU "                                             "     
       %         @                               #                    
       #GRID_MASS $                                             $     
       %         @                               %                    
       #LENGTH &   #Z '                                             &     
                                                 '     
       %         @                               (                    
       #VEL )   #Z *                                             )     
                                                 *     
       #         @                                  +                  #RUN_UNIT_TEST%INT ,   #Z -                               @              ,     INT           D @                              -     
              K      fn#fn    ë   @   J   COMMON_VARS    +  @       C+COMMON_VARS    k  @       H0+COMMON_VARS     «  @       NU0+COMMON_VARS #   ë  @       M_GRID+COMMON_VARS #   +  @       LSCALE+COMMON_VARS #   k  @       TSCALE+COMMON_VARS $   «  @       OMEGA_0+COMMON_VARS $   ë  @       OMEGA_B+COMMON_VARS    +  @       H+COMMON_VARS $   k  @       LAMBDA0+COMMON_VARS (   «  @       CO_BOXWIDTH+COMMON_VARS $   ë  @       BOXSIZE+COMMON_VARS "   +  @       M_SOL+COMMON_VARS >   k  ¤      MPI_CONSTANTS!MPIPRIV1+MPI_CONSTANTS=MPIPRIV1 2     H      MPIPRIV1%MPI_BOTTOM+MPI_CONSTANTS 4   W  H      MPIPRIV1%MPI_IN_PLACE+MPI_CONSTANTS 9     ¤      MPIPRIV1%MPI_STATUS_IGNORE+MPI_CONSTANTS >   C        MPI_CONSTANTS!MPIPRIV2+MPI_CONSTANTS=MPIPRIV2 ;   ×  Ä      MPIPRIV2%MPI_STATUSES_IGNORE+MPI_CONSTANTS ;     ¤      MPIPRIV2%MPI_ERRCODES_IGNORE+MPI_CONSTANTS >   ?        MPI_CONSTANTS!MPIPRIVC+MPI_CONSTANTS=MPIPRIVC 6   È  Ä      MPIPRIVC%MPI_ARGVS_NULL+MPI_CONSTANTS 5   	  ¤      MPIPRIVC%MPI_ARGV_NULL+MPI_CONSTANTS    0
  W       Z_TO_D    
  @   a   Z_TO_D%Z    Ç
  h       D_TO_Z    /  =      D_TO_Z%SQRT    l  @   a   D_TO_Z%D    ¬  i       D_TO_NU      =      D_TO_NU%SQRT    R  @   a   D_TO_NU%D      X       NU_TO_D    ê  @   a   NU_TO_D%NU &   *  _       CONVERT_MASS2PHYSICAL 0     @   a   CONVERT_MASS2PHYSICAL%GRID_MASS (   É  c       CONVERT_LENGTH2PHYSICAL /   ,  @   a   CONVERT_LENGTH2PHYSICAL%LENGTH *   l  @   a   CONVERT_LENGTH2PHYSICAL%Z %   ¬  `       CONVERT_VEL2PHYSICAL )     @   a   CONVERT_VEL2PHYSICAL%VEL '   L  @   a   CONVERT_VEL2PHYSICAL%Z      f       RUN_UNIT_TEST "   ò  <      RUN_UNIT_TEST%INT     .  @   a   RUN_UNIT_TEST%Z 