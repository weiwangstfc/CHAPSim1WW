!***********************************************************************
!> @author 
!> Mehdi Seddighi-Moornani, Wei Wang @ HETF, ME, The University of Sheffield.
!
! DESCRIPTION: 
!> Brief description of subroutine. 
!> @details 
!>  CALCULATION OF COEFF. OF THE POISSON EQ. (Tridiagonal System of Equation)
!> @param AMPH(1:NND2M)   only used in  SUBROUTINE PRCALC  
!> @param APPH(1:NND2M)   only used in  SUBROUTINE PRCALC
!> @param ACPH(1:NND2M)   only used in  SUBROUTINE PRCALC
!> @param AMVR(1:NND2M)   \phi_{i-1} Coefficient of Possion Eq. in Uniform and non-uniform Eq.3.13 Mehdi thesis
!> @param APVR(1:NND2M)   \phi_{i}   Coefficient of Possion Eq. in Uniform and non-uniform Eq.3.13 Mehdi thesis
!> @param ACVR(1:NND2M)   \phi_{i+1} Coefficient of Possion Eq. in Uniform and non-uniform Eq.3.13 Mehdi thesis
!
!> @todo
!> Nothing left to do.
!> @note: For 3D  case with two periodic directons and with the third non-homogenous
!>        dierction, a very good efficiency for the solution of the elliptic equation
!>        requries two consecutive FFTs in the homogeneous direction and a tridiagonal 
!>        solver in the third direction.

!> @note: this subroutine calculates the coefficients for the tridiagonal solver 
!>        for the non-homogenous direction.  See Page 29 of Orlandi Book
!
! REVISION HISTORY:
! ??/??/20??- Initial Version, by Mehdi Seddighi
! 06/02/2014- Code structured is optimized by Wei Wang
!**********************************************************************
    SUBROUTINE TDMA_COEF
!>      TEST the same as Ke Kui version in both channel and pipe.
!>       For WALL b.c.: u, v, w   on the wall
!>                    : thermal properties on the wall
!>                    : pressure on the ghost cells
        use mesh_info
        use init_info
        IMPLICIT NONE     
      
        INTEGER(4) :: JC, JM, JP, IDR
        INTEGER(4) :: NQ
        INTEGER(4) :: NFIL
        REAL(WP)    :: A22ICC
        REAL(WP)    :: A22ICP
        REAL(WP)    :: A22
        REAL(WP)    :: A22DP
        REAL(WP)    :: UCAJ11, UCAJ
        REAL(WP)    :: UGMMM
        REAL(WP)    :: AC2
        
        IF(MYID.NE.0) RETURN

        !================CARTESIAN============================================
        IF(ICASE==ICHANL) THEN
            DO JC=1,NCL2      
                JP=JC+1
                
                AMPH(JC)= DYFI(JC)*DYCI(JC)
                APPH(JC)= DYFI(JC)*DYCI(JP)
                
                IF(JC==1)    AMPH(JC)= DYFI(JC)*DYCI(JC)*0.5_WP ! for b.c. phi in ghost cell
                IF(JP==NND2) APPH(JC)= DYFI(JC)*DYCI(JP)*0.5_WP
                
            ENDDO
            ! set b.c. for phi=====
            !AMPH0 = AMPH(1)
            !APPH0 = APPH(NCL2)
            ACPH(1:NCL2)= -( AMPH(1:NCL2)+APPH(1:NCL2) )
            !================
            !pressure gradient equal to zero at wall
            ACPH(1)    = AMPH(1) + ACPH(1)
            AMPH(1)    = 0.0_WP
          
            ACPH(NCL2) = APPH(NCL2) + ACPH(NCL2)
            APPH(NCL2) = 0.0_WP
          
            !IF(TGFLOWflg .or. (visthemflg==visimplicit)) THEN
            DO IDR=1,3
                IF( (IDR.EQ.1) .OR. (IDR.EQ.3) ) THEN
                    DO JC=2,NCL2-1
                        JP=JC+1
                        AMVR(JC,IDR)= DYFI(JC)*DYCI(JC)
                        APVR(JC,IDR)= DYFI(JC)*DYCI(JP)
                        ACVR(JC,IDR)=-( AMVR(JC,IDR) + APVR(JC,IDR) )
                    ENDDO
    !>              @note : FOR JC=1   ! to do for u,w on the wall!!! check
                    UCAJ=4.0_WP/( 2.0_WP/DYCI(2)+1.0_WP/DYFI(1) )
                    AMVR(1,IDR)=   0.0_WP
                    APVR(1,IDR)=   UCAJ*DYCI(2)
                    ACVR(1,IDR)=  -UCAJ*( DYCI(2)+2.0_WP*DYFI(1) )
                    AMVR1 = UCAJ*(2.0_WP*DYFI(1))
    !>               @note : FOR JC=NND2M ! to do for u,w on the wall!!!  check
                    UCAJ=4.0_WP/( 2.0_WP/DYCI(NCL2)+1.0_WP/DYFI(NCL2) )
                    AMVR(NCL2,IDR)= UCAJ*DYCI(NCL2)
                    APVR(NCL2,IDR)= 0.0_WP
                    ACVR(NCL2,IDR)=-UCAJ*( 2.0_WP*DYFI(NCL2)+DYCI(NCL2) )
                    APVRN = UCAJ*(2.0_WP*DYFI(NCL2))
    
                ELSE
    !>          @note: clustered grids in non-homogenous direction (y) , the same methods as above.
    !>                 coefficients in Eq.3.13 Page 44 of Mehdi thesis.
                    DO JC=2,NCL2
                        JM=JC-1
                        AMVR(JC,IDR)=  DYFI(JM)*DYCI(JC)
                        APVR(JC,IDR)=  DYFI(JC)*DYCI(JC)
                        ACVR(JC,IDR)= -( AMVR(JC,IDR)+APVR(JC,IDR) )
                    ENDDO
                    AMVR(1,IDR)=0.0_WP
                    APVR(1,IDR)=0.0_WP
                    ACVR(1,IDR)=1.0_WP
                    AMVR(NND2,IDR)=0.0_WP
                    APVR(NND2,IDR)=0.0_WP
                    ACVR(NND2,IDR)=1.0_WP
                
                ENDIF
             
            ENDDO
            !END IF
            
        END IF
      
        
        !================CYLINDRICAL============================================
        IF(ICASE==IPIPEC .OR. ICASE==IANNUL) THEN
            ! for pressure correction==========
            do jc=2,NCL2-1
                jp=jc+1
                A22ICC=DYCI(jc)/RNDI1(jc) !rc(jC)
                A22ICP=DYCI(jp)/RNDI1(JP) !rc(jp)
                AC2=-(A22ICC+A22ICP)
                UGMMM=DYFI(JC)/RCCI1(JC)       !rm(jc)*     
                AMPH(JC)=(A22ICC)*UGMMM 
                APPH(JC)=(A22ICP)*UGMMM 
                ACPH(JC)=AC2*UGMMM 
            enddo
          
            JC=1
            JP=JC+1
            A22ICP=DBLE(jp-jc)*DYCI(jp)/RNDI1(JP)     !rc(jp)*      
            UGMMM=DYFI(JC)/RCCI1(JC) !rm(jc)*
            AMPH(JC)=0.0_WP
            APPH(JC)= (A22ICP)*UGMMM         
            ACPH(JC)=-(A22ICP)*UGMMM        
             
            JC=NCL2
            A22ICC=DYCI(JC)/RNDI1(JC)  !rc(jc)*
            UGMMM =DYFI(JC)/RCCI1(JC)  !rm(jc)*
            AMPH(JC)= (A22ICC)*UGMMM        
            APPH(JC)=0.0_WP
            ACPH(JC)=-(A22ICC)*UGMMM      
                
            !IF(TGFLOWflg .or. (visthemflg==visimplicit)) THEN
            do nq=1,3
                if(nq.eq.3) then                   
                    do jc=2,NCL2-1
                        jp=jc+1
                        jm=jc-1     !@
                        A22ICC=DYCI(JC)/RNDI1(JC)          !rc(jc)*
                        A22ICP=DYCI(JP)/RNDI1(JP)          !rc(jp)* 
                        AMVR(JC,NQ)=  A22ICC*DYFI(JC)*RCCI1(JM)      !/rm(jm)            
                        APVR(JC,NQ)=  A22ICP*DYFI(JC)*RCCI1(JP)      !/rm(jp)             
                        ACVR(JC,NQ)=-(A22ICC+A22ICP)*DYFI(JC)*RCCI1(jc)-RCCI2(JC)!1.0_WP/rm(jc)**2    
                    enddo
    
                    JC=1
                    JP=JC+1                 
                    UCAJ11=DYFI(JC)*DYCI(JP)/RNDI1(JP)      !*rc(jp)          
                    AMVR(JC,NQ)=0.0_WP                   
                    APVR(JC,NQ)= UCAJ11*RCCI1(JP)!/rm(jp)               
                    ACVR(JC,NQ)=-UCAJ11*RCCI1(jc)-RCCI2(JC)!1.0_WP/rm(jc)**2       
    
                    JC=NCL2
                    JP=NND2
                    jm=jc-1 
                    UGMMM=DYFI(JC)/RNDI1(JC)!*rm(jc)       
                                           
                    AMVR(JC,NQ)=UGMMM*DYCI(jc)*RNDI1(JC)    !/rc(jc)
                    APVR(JC,NQ)=0. !-UGMMM/rc(jp)/cac(jp)*2.0_WP                    
                    !ACVR(JC,NQ)=-UGMMM*DYCI(jc)*RNDI1(JC)-UGMMM*DYCI(jp)*2.0_WP*RNDI1(JP)   ! W on the ghost cell  
                    ACVR(JC,NQ)=-UGMMM*DYCI(jc)*RNDI1(JC)-UGMMM*DYCI(jp)*RNDI1(JP)     ! W on the wall.   
                endif
             
                if(nq.eq.1) then
                    do jc=2,NCL2-1
                        jp=jc+1
                        ugmmm=DYFI(JC)*RCCI1(JC)!/rm(jc)        
                        amvr(JC,NQ)=+ugmmm*DYCI(jc)/RNDI1(JC)   !*rc(jc)        
                        apvr(JC,NQ)=+ugmmm*DYCI(jp)/RNDI1(JP)   !*rc(jp)       
                        acvr(JC,NQ)=-(amvr(JC,NQ)+apvr(JC,NQ))                    
                    enddo
     
                    jc=1
                    jp=jc+1
                    ugmmm=DYFI(jc)*RCCI1(JC)!/rm(jc)            
                    amvr(jc,nq)=0.0_WP                
                    apvr(jc,nq)= ugmmm*(DYCI(jp))/RNDI1(JP)    !rc(jp)*    
                    acvr(jc,nq)=-ugmmm*(DYCI(jp))/RNDI1(JP)      !rc(jp)*       
    
                    jc=NCL2
                    jp=NND2
                    ugmmm=DYFI(jc)*RCCI1(JC)!/rm(jc)                                                                  
                    amvr(jc,nq)=+ugmmm*DYCI(jc)/RNDI1(JC)      !rc(jc)*
                    apvr(jc,nq)=0.0_wp !ugmmm*rc(jp)/cac(jp)*2.0_WP    
                    !acvr(jc,nq)=-ugmmm*DYCI(jc)/RNDI1(jc)-ugmmm*DYCI(jp)*2.0_WP/RNDI1(jp)      ! U on the ghost cell 
                    acvr(jc,nq)=-ugmmm*DYCI(jc)/RNDI1(jc)-ugmmm*DYCI(jp)/RNDI1(jp)  ! U on the wall.         
                endif
             
                if(nq.eq.2) then
                    DO jc=2,NCL2
                        jm=jc-1
                        a22  = DYCI(jc)
                        a22dp= DYCI(jc)/(2.0_WP/RNDI1(jc))
                        apvr(jc,nq)=  a22*DYFI(jc)-a22dp
                        acvr(jc,nq)=-(a22*DYFI(jc)+a22*DYFI(jm))
                        amvr(jc,nq)=  a22*DYFI(jm)+a22dp
                    ENDDO
                     
                    do jc=1,NND2,NCL2
                        apvr(jc,nq)=0.0_WP
                        acvr(jc,nq)=0.0_WP
                        amvr(jc,nq)=0.0_WP
                    enddo
                endif
                                               
            enddo   
            !end if
        END IF
        
        
        
        !================CARTESIAN============================================
        IF(ICASE==IBOX3P) THEN
        
            DO JC=1,NCL2      
                JP=JC+1
                
                AMPH(JC)= DYFI(JC)*DYCI(JC)
                APPH(JC)= DYFI(JC)*DYCI(JP)
                
                IF(JC==1)    AMPH(JC)= DYFI(JC)*DYCI(JC) ! DYCI is whole length, not half at wall
                IF(JP==NND2) APPH(JC)= DYFI(JC)*DYCI(JP) ! DYCI is whole length, not half at wall
                
            ENDDO
            ! set b.c. for phi=====
            ! = AMPH(1)
            !APPH0 = APPH(NCL2)
            ACPH(1:NCL2)= -( AMPH(1:NCL2)+APPH(1:NCL2) )
            
          
            !IF(TGFLOWflg .or. (visthemflg==visimplicit)) THEN
            DO IDR=1,3
                
                IF(IDR==1 .OR. IDR==3) THEN
                    DO JC=2,NCL2-1
                        JP=JC+1
                        AMVR(JC,IDR)= DYFI(JC)*DYCI(JC)
                        APVR(JC,IDR)= DYFI(JC)*DYCI(JP)
                        ACVR(JC,IDR)=-( AMVR(JC,IDR) + APVR(JC,IDR) )
                    ENDDO
                    JC=1
                    AMVR(1,IDR)= DYFI(1)*DYCI(1)
                    APVR(1,IDR)= DYFI(1)*DYCI(2)
                    ACVR(1,IDR)=-( AMVR(1,IDR) + APVR(1,IDR) )
                    JC=NCL2
                    AMVR(NCL2,IDR)= DYFI(NCL2)*DYCI(NCL2)
                    APVR(NCL2,IDR)= DYFI(NCL2)*DYCI(NND2)
                    ACVR(NCL2,IDR)=-( AMVR(NCL2,IDR) + APVR(NCL2,IDR) )
                ELSE
                    DO JC=2,NCL2
                        JM=JC-1
                        AMVR(JC,IDR)=  DYFI(JM)*DYCI(JC)
                        APVR(JC,IDR)=  DYFI(JC)*DYCI(JC)
                        ACVR(JC,IDR)= -( AMVR(JC,IDR)+APVR(JC,IDR) )
                    ENDDO
                    JC=1
                    AMVR(1,IDR)= DYFI(NCL2)*DYCI(1)
                    APVR(1,IDR)= DYFI(1)*DYCI(1)
                    ACVR(1,IDR)=-( AMVR(1,IDR) + APVR(1,IDR) )
                    JC=NND2
                    AMVR(NND2,IDR)=  DYFI(NCL2)*DYCI(NND2)
                    APVR(NND2,IDR)=  DYFI(NCL2)*DYCI(NND2)
                    ACVR(NND2,IDR)= -( AMVR(1,IDR)+APVR(1,IDR) )
                END IF

         
            ENDDO
            !END IF
            
            
!            ACPH(:) = ACVR(:,1)
!            AMPH(:) = AMVR(:,1)
!            APPH(:) = APVR(:,1)
            
        END IF
        
!!**********************************************************************    
        IF(MYID==0) THEN
            NFIL=16
            OPEN(NFIL,FILE='CHK_COEF_TDMA.dat')
            !REWIND NFIL
            WRITE(NFIL,'(A)') &
                   '#AMPH(JC),  ACPH(JC), APPH(JC)'
            DO JC=1,NCL2
                WRITE(NFIL,'(3ES15.7)') &
                    AMPH(JC),  ACPH(JC), APPH(JC)
            END DO
            
            !IF(TGFLOWflg .or. (visthemflg==visimplicit)) THEN
            WRITE(NFIL,'(A)') &
                   '#AMVR(JC,1), ACVR(JC,1), APVR(JC,1)'
            DO JC=1,NCL2
                WRITE(NFIL,'(3ES15.7)') &
                    AMVR(JC,1), ACVR(JC,1), APVR(JC,1)
            END DO
            
            WRITE(NFIL,'(A)') &
                   '#AMVR(JC,2), ACVR(JC,2), APVR(JC,2)'
            DO JC=1,NND2
                WRITE(NFIL,'(3ES15.7)') &
                    AMVR(JC,2), ACVR(JC,2), APVR(JC,2)
            END DO
            
            WRITE(NFIL,'(A)') &
                   '#AMVR(JC,3), ACVR(JC,3), APVR(JC,3)'
            DO JC=1,NCL2
                WRITE(NFIL,'(3ES15.7)') &
                    AMVR(JC,3), ACVR(JC,3), APVR(JC,3)
            END DO
            
!            !IF(TGFLOWflg .or. (visthemflg==visimplicit)) THEN
!            WRITE(NFIL,'(A)') &
!                   '#AMVR(JC,1), ACVR(JC,1), APVR(JC,1), ',  &
!                   '#AMVR(JC,2), ACVR(JC,2), APVR(JC,2), ',  &
!                   '#AMVR(JC,3), ACVR(JC,3), APVR(JC,3)  '
!            DO JC=1,NCL2
!                WRITE(NFIL,'(9ES15.7)') &
!                    AMVR(JC,1), ACVR(JC,1), APVR(JC,1),  &
!                    AMVR(JC,2), ACVR(JC,2), APVR(JC,2),  &
!                    AMVR(JC,3), ACVR(JC,3), APVR(JC,3)
!            END DO
            !END IF
            CLOSE(NFIL)
        END IF
        
        RETURN
      
    END SUBROUTINE

!**********************************************************************

    SUBROUTINE BCAST_TDMA_COEF
        use mesh_info
        use init_info
        IMPLICIT NONE 
        
        CALL MPI_BCAST( AMPH, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ACPH, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( APPH, NCL2, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        !IF(TGFLOWflg .or. (visthemflg==visimplicit)) THEN
        CALL MPI_BCAST( AMVR, NND2*3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( ACVR, NND2*3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        CALL MPI_BCAST( APVR, NND2*3, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        
        !CALL MPI_BCAST( AMPH0, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        !CALL MPI_BCAST( APPH0, 1, MPI_DOUBLE_PRECISION, 0, ICOMM, IERROR )
        !END IF
        
 
        RETURN
    END SUBROUTINE
    
    
