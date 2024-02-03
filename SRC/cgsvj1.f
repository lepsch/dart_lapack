      SUBROUTINE CGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               EPS, SFMIN, TOL
      int                INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK )
      REAL               SVA( N )
      // ..

*  =====================================================================

      // .. Local Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0;
      // ..
      // .. Local Scalars ..
      COMPLEX            AAPQ, OMPQ
      REAL               AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN
      int                BLSKIP, EMPTSW, i, ibr, igl, IERR, IJBLSK, ISWROT, jbc, jgl, KBL, MVL, NOTROT, nblc, nblr, p, PSKIPPED, q, ROWSKIP, SWBAND;
      bool               APPLV, ROTOK, RSVEC;
      // ..
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, CONJG, REAL, MIN, SIGN, SQRT
      // ..
      // .. External Functions ..
      REAL               SCNRM2
      COMPLEX            CDOTC
      int                ISAMAX;
      bool               LSAME;
      // EXTERNAL ISAMAX, LSAME, CDOTC, SCNRM2
      // ..
      // .. External Subroutines ..
      // .. from BLAS
      // EXTERNAL CCOPY, CROT, CSWAP, CAXPY
      // .. from LAPACK
      // EXTERNAL CLASCL, CLASSQ, XERBLA
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
         CALL XERBLA( 'CGSVJ1', -INFO )
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
      // LARGE = BIG / SQRT( REAL( M*N ) )
      BIGTHETA = ONE / ROOTEPS
      ROOTTOL = SQRT( TOL )

      // .. Initialize the right singular vector matrix ..

      // RSVEC = LSAME( JOBV, 'Y' )

      EMPTSW = N1*( N-N1 )
      NOTROT = 0

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
      // if CGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm CGEJSV.


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

      // Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
      // 1 <= p < q <= N. This is the first step toward a blocked implementation
      // of the rotations. New implementation, based on block transformations,
      // is under development.

         DO 2000 ibr = 1, NBLR

            igl = ( ibr-1 )*KBL + 1



* ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1

             // DO 2010 jbc = ibr + 1, NBL
            DO 2010 jbc = 1, NBLC

               jgl = ( jbc-1 )*KBL + N1 + 1

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

         // Safe Gram matrix computation

                           if ( AAQQ.GE.ONE ) {
                              if ( AAPP.GE.AAQQ ) {
                                 ROTOK = ( SMALL*AAPP ).LE.AAQQ
                              } else {
                                 ROTOK = ( SMALL*AAQQ ).LE.AAPP
                              }
                              if ( AAPP.LT.( BIG / AAQQ ) ) {
                                 AAPQ = ( CDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP
                              } else {
                                 CALL CCOPY( M, A( 1, p ), 1, WORK, 1 )                                  CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR )
                                 AAPQ = CDOTC( M, WORK, 1, A( 1, q ), 1 ) / AAQQ
                              }
                           } else {
                              if ( AAPP.GE.AAQQ ) {
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              } else {
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              }
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( CDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / MAX(AAQQ,AAPP) ) / MIN(AAQQ,AAPP)
                              } else {
                                 CALL CCOPY( M, A( 1, q ), 1, WORK, 1 )                                  CALL CLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR )
                                 AAPQ = CDOTC( M, A( 1, p ), 1, WORK, 1 ) / AAPP
                              }
                           }

                            // AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                           AAPQ1  = -ABS(AAPQ)
                           MXAAPQ = MAX( MXAAPQ, -AAPQ1 )

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ABS( AAPQ1 ).GT.TOL ) {
                              OMPQ = AAPQ / ABS(AAPQ)
                              NOTROT = 0
*[RTD]      ROTATED  = ROTATED + 1
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*ABS( AQOAP-APOAQ )/ AAPQ1
                                 IF( AAQQ.GT.AAPP0 )THETA = -THETA

                                 if ( ABS( THETA ).GT.BIGTHETA ) {
                                    T  = HALF / THETA
                                    CS = ONE
                                    CALL CROT( M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T )
                                    if ( RSVEC ) {
                                        CALL CROT( MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T )
                                    }
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, ABS( T ) )
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ1 )
                                    IF( AAQQ.GT.AAPP0 )THSIGN = -THSIGN
                                    T = ONE / ( THETA+THSIGN* SQRT( ONE+THETA*THETA ) )
                                    CS = SQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
                                    MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) )

                                    CALL CROT( M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN )
                                    if ( RSVEC ) {
                                        CALL CROT( MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN )
                                    }
                                 }
                                 D(p) = -D(q) * OMPQ

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                               if ( AAPP.GT.AAQQ ) {
                                    CALL CCOPY( M, A( 1, p ), 1, WORK, 1 )                                     CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK,LDA, IERR )                                     CALL CLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )                                     CALL CAXPY( M, -AAPQ, WORK, 1, A( 1, q ), 1 )                                     CALL CLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                               } else {
                                   CALL CCOPY( M, A( 1, q ), 1, WORK, 1 )                                     CALL CLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK,LDA, IERR )                                     CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR )                                     CALL CAXPY( M, -CONJG(AAPQ), WORK, 1, A( 1, p ), 1 )                                     CALL CLASCL( 'G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR )
                                    SVA( p ) = AAPP*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                               }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // .. recompute SVA(q), SVA(p)
                              if ( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) {
                                    SVA( q ) = SCNRM2( M, A( 1, q ), 1)
                                  } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    CALL CLASSQ( M, A( 1, q ), 1, T, AAQQ )
                                    SVA( q ) = T*SQRT( AAQQ )
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) {
                                 if ( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) {
                                    AAPP = SCNRM2( M, A( 1, p ), 1 )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    CALL CLASSQ( M, A( 1, p ), 1, T, AAPP )
                                    AAPP = T*SQRT( AAPP )
                                 }
                                 SVA( p ) = AAPP
                              }
               // end of OK rotation
                           } else {
                              NOTROT = NOTROT + 1
*[RTD]      SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1
                              IJBLSK = IJBLSK + 1
                           }
                        } else {
                           NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                           IJBLSK = IJBLSK + 1
                        }

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
***
 2000    CONTINUE
*2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ).LT.ROOTBIG ) .AND. ( SVA( N ).GT.ROOTSFMIN ) ) {
            SVA( N ) = SCNRM2( M, A( 1, N ), 1 )
         } else {
            T = ZERO
            AAPP = ONE
            CALL CLASSQ( M, A( 1, N ), 1, T, AAPP )
            SVA( N ) = T*SQRT( AAPP )
         }

      // Additional steering devices

         IF( ( i.LT.SWBAND ) .AND. ( ( MXAAPQ.LE.ROOTTOL ) .OR. ( ISWROT.LE.N ) ) )SWBAND = i

         if ( ( i.GT.SWBAND+1 ) .AND. ( MXAAPQ.LT.SQRT( REAL( N ) )* TOL ) .AND. ( REAL( N )*MXAAPQ*MXSINJ.LT.TOL ) ) {
            GO TO 1994
         }

         IF( NOTROT.GE.EMPTSW )GO TO 1994

 1993 CONTINUE
      // end i=1:NSWEEP loop

* #:( Reaching this point means that the procedure has not converged.
      INFO = NSWEEP - 1
      GO TO 1995

 1994 CONTINUE
* #:) Reaching this point means numerical convergence after the i-th
      // sweep.

      INFO = 0
* #:) INFO = 0 confirms successful iterations.
 1995 CONTINUE

      // Sort the vector SVA() of column norms.
      DO 5991 p = 1, N - 1
         q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         if ( p.NE.q ) {
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            AAPQ = D( p )
            D( p ) = D( q )
            D( q ) = AAPQ
            CALL CSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
            IF( RSVEC )CALL CSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
         }
 5991 CONTINUE


      RETURN
      // ..
      // .. END OF CGSVJ1
      // ..
      }
