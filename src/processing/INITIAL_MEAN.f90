    SUBROUTINE INITIAL_MEAN_TG
        use init_info
        use mesh_info
        use flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE 
       
        INTEGER(4) :: IDR, I, J, K, IP, KP, JJ
        REAL(WP)   :: UUJJ, VELENTER
       
        UU(:,:,1) = 0.0_WP
        DO IDR=1,NDV     
       
            IF (IDR.EQ.1) THEN 
                 
                DO J=1,N2DO(MYID) 
                    JJ=JCL2G(J)               
                    UUJJ    =0.0_WP
                    VELENTER=0.0_WP    
                            
                    DO I=1,NCL1_TG
                        IP=IPV_TG(I)        
                        DO K=1,NCL3                                
                            VELENTER= ( Q_tg(IP,J,K,IDR) + Q_tg(I,J,K,IDR) ) * 0.5_WP     
                            UUJJ=VELENTER+UUJJ
                        ENDDO
                   
                    ENDDO                 
                    UU(J,IDR,1)  =UUJJ/DBLE(NCL1_TG*NCL3)
                ENDDO
              
            ELSE IF (IDR.EQ.3) THEN
           
                DO J=1,N2DO(MYID) 
                    JJ=JCL2G(J)        
                    UUJJ    =0.0_WP
                    VELENTER=0.0_WP    
                            
                    DO I=1,NCL1_TG
                         
                        DO K=1,NCL3    
                            KP=KPV(K)                                   
                            !VELENTER= ( Q_tg(I,J,KP,IDR)*RCCI1(JJ) + Q_tg(I,J,K,IDR)*RCCI1(JJ) ) * 0.5_WP 
                            VELENTER= ( Q_tg(I,J,KP,IDR) + Q_tg(I,J,K,IDR) ) * 0.5_WP    
                            UUJJ=VELENTER+UUJJ
                        ENDDO
                   
                    ENDDO                 
                    UU(J,IDR,1)  =UUJJ/DBLE(NCL1_TG*NCL3)
                ENDDO
           
            ELSE IF (IDR.EQ.2) THEN
!                NYI = 1
!                IF(MYID==0 .AND. ICASE==IPIPEC) THEN
!                    J = 1
!                    JJ=  JCL2G(J)
!                    JP=  JLPV(J)
!                    JJP= JGPV(JJ)
!                    UUJJ    =0.0_WP
!                    VELENTER=0.0_WP           
!                    DO I=1,NCL1_TG      
!                        DO K=1,NCL3  
!                            KS = KSYM(K)
!                            !V1=(Q_tg(I,JP,K,IDR) - Q_TG(I,JP,KS,IDR))*RNDI1(JJP)/2.0_WP
!                            !VELENTER= ( Q_tg(I,JP,K,IDR)*RNDI1(JJP) + V1 ) * 0.5_WP
                            
!                            UUJJ=VELENTER+UUJJ
!                        ENDDO
!                    ENDDO                  
!                    UU(J,IDR,1)  =UUJJ/DBLE(NCL1_TG*NCL3) 
!                    NYI = 2
!                END IF
                
!                DO J=NYI,N2DO(MYID)  
!                    JJ=JCL2G(J) 
!                    JP=JLPV(J)
!                    UUJJ    =0.0_WP
!                    VELENTER=0.0_WP           
!                    DO I=1,NCL1_TG      
!                        DO K=1,NCL3                               
!                            VELENTER= ( Q_tg(I,JP,K,IDR)*RNDI1(JJ) + Q_tg(I,J,K,IDR)*RNDI1(JJ) ) * 0.5_WP
!                            UUJJ=VELENTER+UUJJ
!                        ENDDO
!                    ENDDO                  
!                    UU(J,IDR,1)  =UUJJ/DBLE(NCL1_TG*NCL3) 
!                ENDDO  
                
                DO J=1,N2DO(MYID)             
                    UUJJ    =0.0_WP
                    VELENTER=0.0_WP           
                    DO I=1,NCL1_TG      
                        DO K=1,NCL3                               
                            VELENTER= Q_TG(I,J,K,IDR)    
                            UUJJ=VELENTER+UUJJ
                        ENDDO
                    ENDDO                  
                    UU(J,IDR,1)  =UUJJ/DBLE(NCL1_TG*NCL3) 
                ENDDO      
                
            ELSE
                CALL ERRHDL('DIRECTIONS SHOULD BE 1, 2 or 3, and no others.',myid)  
            END IF
        END DO
       
      
        RETURN
      
    END SUBROUTINE
      
!***********************************************************************************
    SUBROUTINE INITIAL_MEAN_IO
        use init_info
        use mesh_info
        use flow_info
        USE POSTPROCESS_INFO
        IMPLICIT NONE 
       
        INTEGER(4) :: IDR, I, J, K, IP, KP, JJ
        REAL(WP)   :: UUJJ, VELENTER
       
        UU(:,:,2) = 0.0_WP
        DO IDR=1,NDV     
       
            IF (IDR.EQ.1) THEN 
                 
                DO J=1,N2DO(MYID)              
                    UUJJ    =0.0_WP
                    VELENTER=0.0_WP    
                            
                    DO I=1,NCL1_IO
                        IP=IPV_IO(I)        
                        DO K=1,NCL3                                
                            VELENTER= ( Q_IO(IP,J,K,IDR) + Q_IO(I,J,K,IDR) ) * 0.5_WP     
                            UUJJ=VELENTER+UUJJ
                        ENDDO
                   
                    ENDDO                 
                    UU(J,IDR,2)  =UUJJ/DBLE(NCL1_IO*NCL3)
                ENDDO
              
            ELSE IF (IDR.EQ.3) THEN
           
                DO J=1,N2DO(MYID) 
                    JJ=JCL2G(J)
                    UUJJ    =0.0_WP
                    VELENTER=0.0_WP    
                    DO I=1,NCL1_IO
                        DO K=1,NCL3    
                            KP=KPV(K)                                   
                            !VELENTER= ( Q_IO(I,J,KP,IDR)*RCCI1(JJ) + Q_IO(I,J,K,IDR)*RCCI1(JJ) ) * 0.5_WP     
                            VELENTER= ( Q_IO(I,J,KP,IDR) + Q_IO(I,J,K,IDR) ) * 0.5_WP     
                            UUJJ=VELENTER+UUJJ
                        ENDDO
                   
                    ENDDO                 
                    UU(J,IDR,2)  =UUJJ/DBLE(NCL1_IO*NCL3)
                ENDDO
           
           ELSE IF (IDR.EQ.2) THEN
      
!                NYI = 1
                
!                IF(MYID==0 .AND. ICASE==IPIPEC) THEN
!                    J = 1
!                    JJ=  JCL2G(J)
!                    JP=  JLPV(J)
!                    JJP= JGPV(JJ)
!                    UUJJ    =0.0_WP
!                    VELENTER=0.0_WP           
!                    DO I=1,NCL1_IO      
!                        DO K=1,NCL3  
!                            KS = KSYM(K)
!                            V1=(Q_IO(I,JP,K,IDR) - Q_IO(I,JP,KS,IDR))*RNDI1(JJP)/2.0_WP
!                            VELENTER= ( Q_IO(I,JP,K,IDR)*RNDI1(JJP) + V1 ) * 0.5_WP
!                            UUJJ=VELENTER+UUJJ
!                        ENDDO
!                    ENDDO                  
!                    UU(J,IDR,2)  =UUJJ/DBLE(NCL1_IO*NCL3) 
!                    NYI =  2
!                END IF 
                
!                DO J=NYI,N2DO(MYID)  
!                    JJ=JCL2G(J)    
!                    JP =JLPV(J)   
!                    UUJJ    =0.0_WP
!                    VELENTER=0.0_WP           
!                    DO I=1,NCL1_IO    
!                        DO K=1,NCL3                               
!                            VELENTER= ( Q_IO(I,JP,K,IDR)*RNDI1(JJ) + Q_IO(I,J,K,IDR)*RNDI1(JJ) ) * 0.5_WP
!                            UUJJ=VELENTER+UUJJ
!                        ENDDO
!                    ENDDO                  
!                    UU(J,IDR,2)  =UUJJ/DBLE(NCL1_IO*NCL3) 
!                ENDDO     

                DO J=1,N2DO(MYID)             
                    UUJJ    =0.0_WP
                    VELENTER=0.0_WP           
                    DO I=1,NCL1_IO     
                        DO K=1,NCL3                               
                            VELENTER= Q_IO(I,J,K,IDR)    
                            UUJJ=VELENTER+UUJJ
                        ENDDO
                    ENDDO                  
                    UU(J,IDR,2)  =UUJJ/DBLE(NCL1_IO*NCL3) 
                ENDDO   
                
                
                
           ELSE
              CALL ERRHDL('DIRECTIONS SHOULD BE 1, 2 or 3, and no others.',myid)  
           END IF
       END DO
       
      
      RETURN
      
      END SUBROUTINE
