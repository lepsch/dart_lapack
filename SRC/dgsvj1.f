      SUBROUTINE DGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             EPS, SFMIN, TOL;
      int                INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( N ), SVA( N ), V( LDV, * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Local Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      double             AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, LARGE, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN;
      int                BLSKIP, EMPTSW, i, ibr, igl, IERR, IJBLSK, ISWROT, jbc, jgl, KBL, MVL, NOTROT, nblc, nblr, p, PSKIPPED, q, ROWSKIP, SWBAND;
      bool               APPLV, ROTOK, RSVEC;
      // ..
      // .. Local Arrays ..
      double             FASTR( 5 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DABS, MAX, DBLE, MIN, DSIGN, DSQRT
      // ..
      // .. External Functions ..
      double             DDOT, DNRM2;
      int                IDAMAX;
      bool               LSAME;
      // EXTERNAL IDAMAX, LSAME, DDOT, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DLASCL, DLASSQ, DROTM, DSWAP, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      APPLV = LSAME( JOBV, 'A' )
      RSVEC = LSAME( JOBV, 'V' )
      IF( .NOT.( RSVEC .OR. APPLV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( N.LT.0 ) .OR. ( N.GT.M ) ) THEN
         INFO = -3
      ELSE IF( N1.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.M ) THEN
         INFO = -6
      ELSE IF( ( RSVEC.OR.APPLV ) .AND. ( MV.LT.0 ) ) THEN
         INFO = -9
      ELSE IF( ( RSVEC.AND.( LDV.LT.N ) ).OR. ( APPLV.AND.( LDV.LT.MV ) )  ) THEN
         INFO = -11
      ELSE IF( TOL.LE.EPS ) THEN
         INFO = -14
      ELSE IF( NSWEEP.LT.0 ) THEN
         INFO = -15
      ELSE IF( LWORK.LT.M ) THEN
         INFO = -17
      ELSE
         INFO = 0
      END IF

      // #:(
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGSVJ1', -INFO )
         RETURN
      END IF

      IF( RSVEC ) THEN
         MVL = N
      ELSE IF( APPLV ) THEN
         MVL = MV
      END IF
      RSVEC = RSVEC .OR. APPLV

      ROOTEPS = DSQRT( EPS )
      ROOTSFMIN = DSQRT( SFMIN )
      SMALL = SFMIN / EPS
      BIG = ONE / SFMIN
      ROOTBIG = ONE / ROOTSFMIN
      LARGE = BIG / DSQRT( DBLE( M*N ) )
      BIGTHETA = ONE / ROOTEPS
      ROOTTOL = DSQRT( TOL )

      // .. Initialize the right singular vector matrix ..

      // RSVEC = LSAME( JOBV, 'Y' )

      EMPTSW = N1*( N-N1 )
      NOTROT = 0
      FASTR( 1 ) = ZERO

      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

      KBL = MIN( 8, N )
      NBLR = N1 / KBL
      IF( ( NBLR*KBL ).NE.N1 )NBLR = NBLR + 1

      // .. the tiling is nblr-by-nblc [tiles]

      NBLC = ( N-N1 ) / KBL
      IF( ( NBLC*KBL ).NE.( N-N1 ) )NBLC = NBLC + 1
      BLSKIP = ( KBL**2 ) + 1
*[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = MIN( 5, KBL )
*[TP] ROWSKIP is a tuning parameter.
      SWBAND = 0
*[TP] SWBAND is a tuning parameter. It is meaningful and effective
      // if SGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm SGESVJ.


      // | *   *   * [x] [x] [x]|
      // | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
      // | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
      // |[x] [x] [x] *   *   * |
      // |[x] [x] [x] *   *   * |
      // |[x] [x] [x] *   *   * |


      DO 1993 i = 1, NSWEEP
      // .. go go go ...

         MXAAPQ = ZERO
         MXSINJ = ZERO
         ISWROT = 0

         NOTROT = 0
         PSKIPPED = 0

         DO 2000 ibr = 1, NBLR

            igl = ( ibr-1 )*KBL + 1


*........................................................
* ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1

            DO 2010 jbc = 1, NBLC

               jgl = N1 + ( jbc-1 )*KBL + 1

         // doing the block at ( ibr, jbc )

               IJBLSK = 0
               DO 2100 p = igl, MIN( igl+KBL-1, N1 )

                  AAPP = SVA( p )

                  IF( AAPP.GT.ZERO ) THEN

                     PSKIPPED = 0

                     DO 2200 q = jgl, MIN( jgl+KBL-1, N )

                        AAQQ = SVA( q )

                        IF( AAQQ.GT.ZERO ) THEN
                           AAPP0 = AAPP

      // .. M x 2 Jacobi SVD ..

         // .. Safe Gram matrix computation ..

                           IF( AAQQ.GE.ONE ) THEN
                              IF( AAPP.GE.AAQQ ) THEN
                                 ROTOK = ( SMALL*AAPP ).LE.AAQQ
                              ELSE
                                 ROTOK = ( SMALL*AAQQ ).LE.AAPP
                              END IF
                              IF( AAPP.LT.( BIG / AAQQ ) ) THEN
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP
                              ELSE
                                 CALL DCOPY( M, A( 1, p ), 1, WORK, 1 )
                                 CALL DLASCL( 'G', 0, 0, AAPP, D( p ), M, 1, WORK, LDA, IERR )                                  AAPQ = DDOT( M, WORK, 1, A( 1, q ), 1 )*D( q ) / AAQQ
                              END IF
                           ELSE
                              IF( AAPP.GE.AAQQ ) THEN
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              ELSE
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              END IF
                              IF( AAPP.GT.( SMALL / AAQQ ) ) THEN
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP
                              ELSE
                                 CALL DCOPY( M, A( 1, q ), 1, WORK, 1 )
                                 CALL DLASCL( 'G', 0, 0, AAQQ, D( q ), M, 1, WORK, LDA, IERR )                                  AAPQ = DDOT( M, WORK, 1, A( 1, p ), 1 )*D( p ) / AAPP
                              END IF
                           END IF

                           MXAAPQ = MAX( MXAAPQ, DABS( AAPQ ) )

         // TO rotate or NOT to rotate, THAT is the question ...

                           IF( DABS( AAPQ ).GT.TOL ) THEN
                              NOTROT = 0
            // ROTATED  = ROTATED + 1
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1

                              IF( ROTOK ) THEN

                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*DABS(AQOAP-APOAQ) / AAPQ
                                 IF( AAQQ.GT.AAPP0 )THETA = -THETA

                                 IF( DABS( THETA ).GT.BIGTHETA ) THEN
                                    T = HALF / THETA
                                    FASTR( 3 ) = T*D( p ) / D( q )
                                    FASTR( 4 ) = -T*D( q ) / D( p )
                                    CALL DROTM( M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, DABS( T ) )
                                 ELSE

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -DSIGN( ONE, AAPQ )
                                    IF( AAQQ.GT.AAPP0 )THSIGN = -THSIGN
                                    T = ONE / ( THETA+THSIGN* DSQRT( ONE+THETA*THETA ) )
                                    CS = DSQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
                                    MXSINJ = MAX( MXSINJ, DABS( SN ) )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )

                                    APOAQ = D( p ) / D( q )
                                    AQOAP = D( q ) / D( p )
                                    IF( D( p ).GE.ONE ) THEN

                                       IF( D( q ).GE.ONE ) THEN
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q )*CS
                                          CALL DROTM( M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR )
                                       ELSE
                                          CALL DAXPY( M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                           CALL DAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )
                                          IF( RSVEC ) THEN
                                             CALL DAXPY( MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                              CALL DAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )
                                          END IF
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q ) / CS
                                       END IF
                                    ELSE
                                       IF( D( q ).GE.ONE ) THEN
                                          CALL DAXPY( M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                           CALL DAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )
                                          IF( RSVEC ) THEN
                                             CALL DAXPY( MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                              CALL DAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )
                                          END IF
                                          D( p ) = D( p ) / CS
                                          D( q ) = D( q )*CS
                                       ELSE
                                          IF( D( p ).GE.D( q ) ) THEN
                                             CALL DAXPY( M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                              CALL DAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )
                                             D( p ) = D( p )*CS
                                             D( q ) = D( q ) / CS
                                             IF( RSVEC ) THEN
                                                CALL DAXPY( MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                                 CALL DAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )
                                             END IF
                                          ELSE
                                             CALL DAXPY( M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                              CALL DAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )
                                             D( p ) = D( p ) / CS
                                             D( q ) = D( q )*CS
                                             IF( RSVEC ) THEN
                                                CALL DAXPY( MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                                 CALL DAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )
                                             END IF
                                          END IF
                                       END IF
                                    END IF
                                 END IF

                              ELSE
                                 IF( AAPP.GT.AAQQ ) THEN
                                    CALL DCOPY( M, A( 1, p ), 1, WORK, 1 )                                     CALL DLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR )                                     CALL DLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )
                                    TEMP1 = -AAPQ*D( p ) / D( q )
                                    CALL DAXPY( M, TEMP1, WORK, 1, A( 1, q ), 1 )                                     CALL DLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 ELSE
                                    CALL DCOPY( M, A( 1, q ), 1, WORK, 1 )                                     CALL DLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR )                                     CALL DLASCL( 'G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR )
                                    TEMP1 = -AAPQ*D( q ) / D( p )
                                    CALL DAXPY( M, TEMP1, WORK, 1, A( 1, p ), 1 )                                     CALL DLASCL( 'G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR )
                                    SVA( p ) = AAPP*DSQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 END IF
                              END IF
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q)
            // .. recompute SVA(q)
                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) THEN                                     SVA( q ) = DNRM2( M, A( 1, q ), 1 )* D( q )
                                 ELSE
                                    T = ZERO
                                    AAQQ = ONE
                                    CALL DLASSQ( M, A( 1, q ), 1, T, AAQQ )
                                    SVA( q ) = T*DSQRT( AAQQ )*D( q )
                                 END IF
                              END IF
                              IF( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) THEN
                                 IF( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) THEN                                     AAPP = DNRM2( M, A( 1, p ), 1 )* D( p )
                                 ELSE
                                    T = ZERO
                                    AAPP = ONE
                                    CALL DLASSQ( M, A( 1, p ), 1, T, AAPP )
                                    AAPP = T*DSQRT( AAPP )*D( p )
                                 END IF
                                 SVA( p ) = AAPP
                              END IF
               // end of OK rotation
                           ELSE
                              NOTROT = NOTROT + 1
            // SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1
                              IJBLSK = IJBLSK + 1
                           END IF
                        ELSE
                           NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                           IJBLSK = IJBLSK + 1
                        END IF

       // IF ( NOTROT .GE. EMPTSW )  GO TO 2011
                        IF( ( i.LE.SWBAND ) .AND. ( IJBLSK.GE.BLSKIP ) ) THEN
                           SVA( p ) = AAPP
                           NOTROT = 0
                           GO TO 2011
                        END IF
                        IF( ( i.LE.SWBAND ) .AND. ( PSKIPPED.GT.ROWSKIP ) ) THEN
                           AAPP = -AAPP
                           NOTROT = 0
                           GO TO 2203
                        END IF


 2200                CONTINUE
         // end of the q-loop
 2203                CONTINUE

                     SVA( p ) = AAPP

                  ELSE
                     IF( AAPP.EQ.ZERO )NOTROT = NOTROT + MIN( jgl+KBL-1, N ) - jgl + 1
                     IF( AAPP.LT.ZERO )NOTROT = 0
***      IF ( NOTROT .GE. EMPTSW )  GO TO 2011
                  END IF

 2100          CONTINUE
      // end of the p-loop
 2010       CONTINUE
      // end of the jbc-loop
 2011       CONTINUE
*2011 bailed out of the jbc-loop
            DO 2012 p = igl, MIN( igl+KBL-1, N )
               SVA( p ) = DABS( SVA( p ) )
 2012       CONTINUE
***   IF ( NOTROT .GE. EMPTSW ) GO TO 1994
 2000    CONTINUE
*2000 :: end of the ibr-loop

      // .. update SVA(N)
         IF( ( SVA( N ).LT.ROOTBIG ) .AND. ( SVA( N ).GT.ROOTSFMIN ) ) THEN
            SVA( N ) = DNRM2( M, A( 1, N ), 1 )*D( N )
         ELSE
            T = ZERO
            AAPP = ONE
            CALL DLASSQ( M, A( 1, N ), 1, T, AAPP )
            SVA( N ) = T*DSQRT( AAPP )*D( N )
         END IF

      // Additional steering devices

         IF( ( i.LT.SWBAND ) .AND. ( ( MXAAPQ.LE.ROOTTOL ) .OR. ( ISWROT.LE.N ) ) )SWBAND = i
          IF( ( i.GT.SWBAND+1 ) .AND. ( MXAAPQ.LT.DBLE( N )*TOL ) .AND. ( DBLE( N )*MXAAPQ*MXSINJ.LT.TOL ) ) THEN
            GO TO 1994
         END IF


         IF( NOTROT.GE.EMPTSW )GO TO 1994

 1993 CONTINUE
      // end i=1:NSWEEP loop
* #:) Reaching this point means that the procedure has completed the given
      // number of sweeps.
      INFO = NSWEEP - 1
      GO TO 1995
 1994 CONTINUE
* #:) Reaching this point means that during the i-th sweep all pivots were
      // below the given threshold, causing early exit.

      INFO = 0
* #:) INFO = 0 confirms successful iterations.
 1995 CONTINUE

      // Sort the vector D

      DO 5991 p = 1, N - 1
         q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         IF( p.NE.q ) THEN
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            TEMP1 = D( p )
            D( p ) = D( q )
            D( q ) = TEMP1
            CALL DSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
            IF( RSVEC )CALL DSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
         END IF
 5991 CONTINUE

      RETURN
      // ..
      // .. END OF DGSVJ1
      // ..
      }
