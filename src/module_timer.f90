MODULE TIMER  

  IMPLICIT NONE 
  
  INTEGER,SAVE :: t_srt  
  LOGICAL,SAVE :: t_act
  
  CONTAINS

    SUBROUTINE TIMER_START()  
      INTEGER :: dum
      INTEGER :: t_rate,t_max 
      IF ( t_act.EQV. .FALSE. ) THEN 
         t_act = .TRUE. 
         CALL SYSTEM_CLOCK(t_srt,t_rate,t_max)  
      ENDIF
      RETURN 

    END SUBROUTINE TIMER_START

    SUBROUTINE TIMER_FINISH(msg) 
      CHARACTER(LEN=*) :: msg
      INTEGER :: t_rate,t_max,t_fin
      
      IF ( t_act.EQV. .TRUE. ) THEN 
         CALL SYSTEM_CLOCK(t_fin,t_rate,t_max)  
         t_act=.FALSE.
         WRITE(*,*) msg,(1000.*(t_fin-t_srt)/t_rate)/1000.,'s'
101      FORMAT(a30,':',1x,f5.3,'s')
      ELSE 
         STOP 'TIMER ERROR -- TIMER NOT ACTIVE' 
      ENDIF

      RETURN 

    END SUBROUTINE TIMER_FINISH

END MODULE TIMER
