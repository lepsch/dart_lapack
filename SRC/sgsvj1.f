      SUBROUTINE SGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               EPS, SFMIN, TOL
      int                INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), D( N ), SVA( N ), V( LDV, * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Local Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0;
      // ..
      // .. Local Scalars ..
      REAL               AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, LARGE, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN
      int                BLSKIP, EMPTSW, i, ibr, igl, IERR, IJBLSK, ISWROT, jbc, jgl, KBL, MVL, NOTROT, nblc, nblr, p, PSKIPPED, q, ROWSKIP, SWBAND;
      bool               APPLV, ROTOK, RSVEC;
      // ..
      // .. Local Arrays ..
      REAL               FASTR( 5 )
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, FLOAT, MIN, SIGN, SQRT
      // ..
      // .. External Functions ..
      REAL               SDOT, SNRM2
      int                ISAMAX;
      bool               LSAME;
      // EXTERNAL ISAMAX, LSAME, SDOT, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SLASCL, SLASSQ, SROTM, SSWAP, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      APPLV = LSAME( JOBV, 'A' )
      RSVEC = LSAME( JOBV, 'V' )
      if ( .NOT.( RSVEC .OR. APPLV .OR. LSAME( JOBV, 'N' ) ) ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -2
      } else if ( ( N.LT.0 ) .OR. ( N.GT.M ) ) {
         INFO = -3
      } else if ( N1.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.M ) {
         INFO = -6
      } else if ( ( RSVEC.OR.APPLV ) .AND. ( MV.LT.0 ) ) {
         INFO = -9
      } else if ( ( RSVEC.AND.( LDV.LT.N ) ).OR. ( APPLV.AND.( LDV.LT.MV ) )  ) {
         INFO = -11
      } else if ( TOL.LE.EPS ) {
         INFO = -14
      } else if ( NSWEEP.LT.0 ) {
         INFO = -15
      } else if ( LWORK.LT.M ) {
         INFO = -17
      } else {
         INFO = 0
      }

      // #:(
      if ( INFO.NE.0 ) {
         xerbla('SGSVJ1', -INFO );
         RETURN
      }

      if ( RSVEC ) {
         MVL = N
      } else if ( APPLV ) {
         MVL = MV
      }
      RSVEC = RSVEC .OR. APPLV

      ROOTEPS = SQRT( EPS )
      ROOTSFMIN = SQRT( SFMIN )
      SMALL = SFMIN / EPS
      BIG = ONE / SFMIN
      ROOTBIG = ONE / ROOTSFMIN
      LARGE = BIG / SQRT( FLOAT( M*N ) )
      BIGTHETA = ONE / ROOTEPS
      ROOTTOL = SQRT( TOL )

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

                  if ( AAPP.GT.ZERO ) {

                     PSKIPPED = 0

                     DO 2200 q = jgl, MIN( jgl+KBL-1, N )

                        AAQQ = SVA( q )

                        if ( AAQQ.GT.ZERO ) {
                           AAPP0 = AAPP

      // .. M x 2 Jacobi SVD ..

         // .. Safe Gram matrix computation ..

                           if ( AAQQ.GE.ONE ) {
                              if ( AAPP.GE.AAQQ ) {
                                 ROTOK = ( SMALL*AAPP ).LE.AAQQ
                              } else {
                                 ROTOK = ( SMALL*AAQQ ).LE.AAPP
                              }
                              if ( AAPP.LT.( BIG / AAQQ ) ) {
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP
                              } else {
                                 scopy(M, A( 1, p ), 1, WORK, 1 );
                                 CALL SLASCL( 'G', 0, 0, AAPP, D( p ), M, 1, WORK, LDA, IERR )                                  AAPQ = SDOT( M, WORK, 1, A( 1, q ), 1 )*D( q ) / AAQQ
                              }
                           } else {
                              if ( AAPP.GE.AAQQ ) {
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              } else {
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              }
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP
                              } else {
                                 scopy(M, A( 1, q ), 1, WORK, 1 );
                                 CALL SLASCL( 'G', 0, 0, AAQQ, D( q ), M, 1, WORK, LDA, IERR )                                  AAPQ = SDOT( M, WORK, 1, A( 1, p ), 1 )*D( p ) / AAPP
                              }
                           }

                           MXAAPQ = MAX( MXAAPQ, ABS( AAPQ ) )

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ABS( AAPQ ).GT.TOL ) {
                              NOTROT = 0
            // ROTATED  = ROTATED + 1
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*ABS( AQOAP-APOAQ ) / AAPQ
                                 IF( AAQQ.GT.AAPP0 )THETA = -THETA

                                 if ( ABS( THETA ).GT.BIGTHETA ) {
                                    T = HALF / THETA
                                    FASTR( 3 ) = T*D( p ) / D( q )
                                    FASTR( 4 ) = -T*D( q ) / D( p )
                                    srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, ABS( T ) )
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ )
                                    IF( AAQQ.GT.AAPP0 )THSIGN = -THSIGN
                                    T = ONE / ( THETA+THSIGN* SQRT( ONE+THETA*THETA ) )
                                    CS = SQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
                                    MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )

                                    APOAQ = D( p ) / D( q )
                                    AQOAP = D( q ) / D( p )
                                    if ( D( p ).GE.ONE ) {

                                       if ( D( q ).GE.ONE ) {
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q )*CS
                                          srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                       } else {
                                          saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                           CALL SAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          if ( RSVEC ) {
                                             saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                              CALL SAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                          }
                                          D( p ) = D( p )*CS
                                          D( q ) = D( q ) / CS
                                       }
                                    } else {
                                       if ( D( q ).GE.ONE ) {
                                          saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                           CALL SAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          if ( RSVEC ) {
                                             saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                              CALL SAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                          }
                                          D( p ) = D( p ) / CS
                                          D( q ) = D( q )*CS
                                       } else {
                                          if ( D( p ).GE.D( q ) ) {
                                             saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                              CALL SAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             D( p ) = D( p )*CS
                                             D( q ) = D( q ) / CS
                                             if ( RSVEC ) {
                                                saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                                 CALL SAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             }
                                          } else {
                                             saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                              CALL SAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             D( p ) = D( p ) / CS
                                             D( q ) = D( q )*CS
                                             if ( RSVEC ) {
                                                saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                                 CALL SAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             }
                                          }
                                       }
                                    }
                                 }

                              } else {
                                 if ( AAPP.GT.AAQQ ) {
                                    scopy(M, A( 1, p ), 1, WORK, 1 )                                     CALL SLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR )                                     CALL SLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                    TEMP1 = -AAPQ*D( p ) / D( q )
                                    saxpy(M, TEMP1, WORK, 1, A( 1, q ), 1 )                                     CALL SLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 } else {
                                    scopy(M, A( 1, q ), 1, WORK, 1 )                                     CALL SLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR )                                     CALL SLASCL( 'G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR );
                                    TEMP1 = -AAPQ*D( q ) / D( p )
                                    saxpy(M, TEMP1, WORK, 1, A( 1, p ), 1 )                                     CALL SLASCL( 'G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA( p ) = AAPP*SQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q)
            // .. recompute SVA(q)
                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) THEN                                     SVA( q ) = SNRM2( M, A( 1, q ), 1 )* D( q )
                                 } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    slassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ )*D( q )
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) {
                                 IF( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) THEN                                     AAPP = SNRM2( M, A( 1, p ), 1 )* D( p )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    slassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*SQRT( AAPP )*D( p )
                                 }
                                 SVA( p ) = AAPP
                              }
               // end of OK rotation
                           } else {
                              NOTROT = NOTROT + 1
            // SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1
                              IJBLSK = IJBLSK + 1
                           }
                        } else {
                           NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                           IJBLSK = IJBLSK + 1
                        }

       // IF ( NOTROT .GE. EMPTSW )  GO TO 2011
                        if ( ( i.LE.SWBAND ) .AND. ( IJBLSK.GE.BLSKIP ) ) {
                           SVA( p ) = AAPP
                           NOTROT = 0
                           GO TO 2011
                        }
                        if ( ( i.LE.SWBAND ) .AND. ( PSKIPPED.GT.ROWSKIP ) ) {
                           AAPP = -AAPP
                           NOTROT = 0
                           GO TO 2203
                        }


 2200                CONTINUE
         // end of the q-loop
 2203                CONTINUE

                     SVA( p ) = AAPP

                  } else {
                     IF( AAPP.EQ.ZERO )NOTROT = NOTROT + MIN( jgl+KBL-1, N ) - jgl + 1
                     IF( AAPP.LT.ZERO )NOTROT = 0
***      IF ( NOTROT .GE. EMPTSW )  GO TO 2011
                  }

 2100          CONTINUE
      // end of the p-loop
 2010       CONTINUE
      // end of the jbc-loop
 2011       CONTINUE
*2011 bailed out of the jbc-loop
            DO 2012 p = igl, MIN( igl+KBL-1, N )
               SVA( p ) = ABS( SVA( p ) )
 2012       CONTINUE
***   IF ( NOTROT .GE. EMPTSW ) GO TO 1994
 2000    CONTINUE
*2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ).LT.ROOTBIG ) .AND. ( SVA( N ).GT.ROOTSFMIN ) ) {
            SVA( N ) = SNRM2( M, A( 1, N ), 1 )*D( N )
         } else {
            T = ZERO
            AAPP = ONE
            slassq(M, A( 1, N ), 1, T, AAPP );
            SVA( N ) = T*SQRT( AAPP )*D( N )
         }

      // Additional steering devices

         IF( ( i.LT.SWBAND ) .AND. ( ( MXAAPQ.LE.ROOTTOL ) .OR. ( ISWROT.LE.N ) ) )SWBAND = i
          if ( ( i.GT.SWBAND+1 ) .AND. ( MXAAPQ.LT.FLOAT( N )*TOL ) .AND. ( FLOAT( N )*MXAAPQ*MXSINJ.LT.TOL ) ) {
            GO TO 1994
         }


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
         q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         if ( p.NE.q ) {
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            TEMP1 = D( p )
            D( p ) = D( q )
            D( q ) = TEMP1
            sswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            IF( RSVEC )CALL SSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
         }
 5991 CONTINUE

      RETURN
      // ..
      // .. END OF SGSVJ1
      // ..
      }
