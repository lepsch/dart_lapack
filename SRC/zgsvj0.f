      SUBROUTINE ZGSVJ0( JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE
      // .. Scalar Arguments ..
      int                INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP;
      double             EPS, SFMIN, TOL;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK )
      double             SVA( N );
      // ..

*  =====================================================================

      // .. Local Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0;
      COMPLEX*16   CZERO,                  CONE
      const      CZERO = (0.0D0, 0.0D0), CONE = (1.0D0, 0.0D0) ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         AAPQ, OMPQ
      double             AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN;
      int                BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND;
      bool               APPLV, ROTOK, RSVEC;
      // ..
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, CONJG, DBLE, MIN, SIGN, SQRT
      // ..
      // .. External Functions ..
      double             DZNRM2;
      COMPLEX*16         ZDOTC
      int                IDAMAX;
      bool               LSAME;
      // EXTERNAL IDAMAX, LSAME, ZDOTC, DZNRM2
      // ..
      // ..
      // .. External Subroutines ..
      // ..
      // from BLAS
      // EXTERNAL ZCOPY, ZROT, ZSWAP, ZAXPY
      // from LAPACK
      // EXTERNAL ZLASCL, ZLASSQ, XERBLA
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
      } else if ( LDA.LT.M ) {
         INFO = -5
      } else if ( ( RSVEC.OR.APPLV ) .AND. ( MV.LT.0 ) ) {
         INFO = -8
      } else if ( ( RSVEC.AND.( LDV.LT.N ) ).OR. ( APPLV.AND.( LDV.LT.MV ) ) ) {
         INFO = -10
      } else if ( TOL.LE.EPS ) {
         INFO = -13
      } else if ( NSWEEP.LT.0 ) {
         INFO = -14
      } else if ( LWORK.LT.M ) {
         INFO = -16
      } else {
         INFO = 0
      }

      // #:(
      if ( INFO.NE.0 ) {
         xerbla('ZGSVJ0', -INFO );
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
      BIGTHETA = ONE / ROOTEPS
      ROOTTOL = SQRT( TOL )

      // .. Row-cyclic Jacobi SVD algorithm with column pivoting ..

      EMPTSW = ( N*( N-1 ) ) / 2
      NOTROT = 0

      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..


      SWBAND = 0
*[TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
      // if ZGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure
      // works on pivots inside a band-like region around the diagonal.
      // The boundaries are determined dynamically, based on the number of
      // pivots above a threshold.

      KBL = MIN( 8, N )
*[TP] KBL is a tuning parameter that defines the tile size in the
      // tiling of the p-q loops of pivot pairs. In general, an optimal
      // value of KBL depends on the matrix dimensions and on the
      // parameters of the computer's memory.

      NBL = N / KBL
      IF( ( NBL*KBL ).NE.N )NBL = NBL + 1

      BLSKIP = KBL**2
*[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = MIN( 5, KBL )
*[TP] ROWSKIP is a tuning parameter.

      LKAHEAD = 1
*[TP] LKAHEAD is a tuning parameter.

      // Quasi block transformations, using the lower (upper) triangular
      // structure of the input matrix. The quasi-block-cycling usually
      // invokes cubic convergence. Big part of this cycle is done inside
      // canonical subspaces of dimensions less than M.


      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

      for (i = 1; i <= NSWEEP; i++) { // 1993

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

         for (ibr = 1; ibr <= NBL; ibr++) { // 2000

            igl = ( ibr-1 )*KBL + 1

            DO 1002 ir1 = 0, MIN( LKAHEAD, NBL-ibr )

               igl = igl + ir1*KBL

               DO 2001 p = igl, MIN( igl+KBL-1, N-1 )

      // .. de Rijk's pivoting

                  q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
                  if ( p.NE.q ) {
                     zswap(M, A( 1, p ), 1, A( 1, q ), 1 );
                     IF( RSVEC )CALL ZSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
                     TEMP1 = SVA( p )
                     SVA( p ) = SVA( q )
                     SVA( q ) = TEMP1
                     AAPQ = D(p)
                     D(p) = D(q)
                     D(q) = AAPQ
                  }

                  if ( ir1.EQ.0 ) {

         // Column norms are periodically updated by explicit
         // norm computation.
         // Caveat:
         // Unfortunately, some BLAS implementations compute SNCRM2(M,A(1,p),1)
         // as SQRT(S=ZDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to
         // overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to
         // underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
         // Hence, DZNRM2 cannot be trusted, not even in the case when
         // the true norm is far from the under(over)flow boundaries.
         // If properly implemented DZNRM2 is available, the IF-THEN-ELSE-END IF
         // below should be replaced with "AAPP = DZNRM2( M, A(1,p), 1 )".

                     if ( ( SVA( p ).LT.ROOTBIG ) .AND. ( SVA( p ).GT.ROOTSFMIN ) ) {
                        SVA( p ) = DZNRM2( M, A( 1, p ), 1 )
                     } else {
                        TEMP1 = ZERO
                        AAPP = ONE
                        zlassq(M, A( 1, p ), 1, TEMP1, AAPP );
                        SVA( p ) = TEMP1*SQRT( AAPP )
                     }
                     AAPP = SVA( p )
                  } else {
                     AAPP = SVA( p )
                  }

                  if ( AAPP.GT.ZERO ) {

                     PSKIPPED = 0

                     DO 2002 q = p + 1, MIN( igl+KBL-1, N )

                        AAQQ = SVA( q )

                        if ( AAQQ.GT.ZERO ) {

                           AAPP0 = AAPP
                           if ( AAQQ.GE.ONE ) {
                              ROTOK = ( SMALL*AAPP ).LE.AAQQ
                              if ( AAPP.LT.( BIG / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP
                              } else {
                                 CALL ZCOPY( M, A( 1, p ), 1, WORK, 1 )                                  CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR )                                  AAPQ = ZDOTC( M, WORK, 1, A( 1, q ), 1 ) / AAQQ
                              }
                           } else {
                              ROTOK = AAPP.LE.( AAQQ / SMALL )
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAPP ) / AAQQ
                              } else {
                                 zcopy(M, A( 1, q ), 1, WORK, 1 )                                  CALL ZLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR );
                                 AAPQ = ZDOTC( M, A( 1, p ), 1, WORK, 1 ) / AAPP
                              }
                           }

                            // AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q)
                           AAPQ1  = -ABS(AAPQ)
                           MXAAPQ = MAX( MXAAPQ, -AAPQ1 )

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ABS( AAPQ1 ).GT.TOL ) {
                              OMPQ = AAPQ / ABS(AAPQ)

            // .. rotate
*[RTD]      ROTATED = ROTATED + ONE

                              if ( ir1.EQ.0 ) {
                                 NOTROT = 0
                                 PSKIPPED = 0
                                 ISWROT = ISWROT + 1
                              }

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*ABS( AQOAP-APOAQ )/AAPQ1

                                 if ( ABS( THETA ).GT.BIGTHETA ) {

                                    T  = HALF / THETA
                                    CS = ONE
                                     zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
                                    }
                                     SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, ABS( T ) )

                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ1 )
                                    T = ONE / ( THETA+THSIGN* SQRT( ONE+THETA*THETA ) )
                                    CS = SQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS

                                    MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) )

                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 D(p) = -D(q) * OMPQ

                                 } else {
               // .. have to use modified Gram-Schmidt like transformation
                                 zcopy(M, A( 1, p ), 1, WORK, 1 )                                  CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR );
                                 zlascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )                                  CALL ZAXPY( M, -AAPQ, WORK, 1, A( 1, q ), 1 )                                  CALL ZLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )                                  SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // recompute SVA(q), SVA(p).

                              if ( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) {
                                    SVA( q ) = DZNRM2( M, A( 1, q ), 1 )
                                 } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    zlassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ )
                                 }
                              }
                              if ( ( AAPP / AAPP0 ).LE.ROOTEPS ) {
                                 if ( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) {
                                    AAPP = DZNRM2( M, A( 1, p ), 1 )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    zlassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*SQRT( AAPP )
                                 }
                                 SVA( p ) = AAPP
                              }

                           } else {
         // A(:,p) and A(:,q) already numerically orthogonal
                              IF( ir1.EQ.0 )NOTROT = NOTROT + 1
*[RTD]      SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1
                           }
                        } else {
         // A(:,q) is zero column
                           IF( ir1.EQ.0 )NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                        }

                        if ( ( i.LE.SWBAND ) .AND. ( PSKIPPED.GT.ROWSKIP ) ) {
                           IF( ir1.EQ.0 )AAPP = -AAPP
                           NOTROT = 0
                           GO TO 2103
                        }

 2002                CONTINUE
      // END q-LOOP

 2103                CONTINUE
      // bailed out of q-loop

                     SVA( p ) = AAPP

                  } else {
                     SVA( p ) = AAPP
                     IF( ( ir1.EQ.0 ) .AND. ( AAPP.EQ.ZERO ) ) NOTROT = NOTROT + MIN( igl+KBL-1, N ) - p
                  }

 2001          CONTINUE
      // end of the p-loop
      // end of doing the block ( ibr, ibr )
 1002       CONTINUE
      // end of ir1-loop

* ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1

            DO 2010 jbc = ibr + 1, NBL

               jgl = ( jbc-1 )*KBL + 1

         // doing the block at ( ibr, jbc )

               IJBLSK = 0
               DO 2100 p = igl, MIN( igl+KBL-1, N )

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
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP
                              } else {
                                 zcopy(M, A( 1, p ), 1, WORK, 1 )                                  CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR );
                                 AAPQ = ZDOTC( M, WORK, 1, A( 1, q ), 1 ) / AAQQ
                              }
                           } else {
                              if ( AAPP.GE.AAQQ ) {
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              } else {
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              }
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / MAX(AAQQ,AAPP) ) / MIN(AAQQ,AAPP)
                              } else {
                                 zcopy(M, A( 1, q ), 1, WORK, 1 )                                  CALL ZLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR );
                                 AAPQ = ZDOTC( M, A( 1, p ), 1, WORK, 1 ) / AAPP
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
                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
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

                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 D(p) = -D(q) * OMPQ

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                               if ( AAPP.GT.AAQQ ) {
                                    zcopy(M, A( 1, p ), 1, WORK, 1 )                                     CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK,LDA, IERR )                                     CALL ZLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )                                     CALL ZAXPY( M, -AAPQ, WORK, 1, A( 1, q ), 1 )                                     CALL ZLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                               } else {
                                   zcopy(M, A( 1, q ), 1, WORK, 1 )                                     CALL ZLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK,LDA, IERR )                                     CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR )                                     CALL ZAXPY( M, -CONJG(AAPQ), WORK, 1, A( 1, p ), 1 )                                     CALL ZLASCL( 'G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA( p ) = AAPP*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                               }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // .. recompute SVA(q), SVA(p)
                              if ( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) {
                                    SVA( q ) = DZNRM2( M, A( 1, q ), 1)
                                  } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    zlassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ )
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) {
                                 if ( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) {
                                    AAPP = DZNRM2( M, A( 1, p ), 1 )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    zlassq(M, A( 1, p ), 1, T, AAPP );
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
            SVA( N ) = DZNRM2( M, A( 1, N ), 1 )
         } else {
            T = ZERO
            AAPP = ONE
            zlassq(M, A( 1, N ), 1, T, AAPP );
            SVA( N ) = T*SQRT( AAPP )
         }

      // Additional steering devices

         IF( ( i.LT.SWBAND ) .AND. ( ( MXAAPQ.LE.ROOTTOL ) .OR. ( ISWROT.LE.N ) ) )SWBAND = i

         if ( ( i.GT.SWBAND+1 ) .AND. ( MXAAPQ.LT.SQRT( DBLE( N ) )* TOL ) .AND. ( DBLE( N )*MXAAPQ*MXSINJ.LT.TOL ) ) {
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
         q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         if ( p.NE.q ) {
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            AAPQ = D( p )
            D( p ) = D( q )
            D( q ) = AAPQ
            zswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            IF( RSVEC )CALL ZSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
         }
 5991 CONTINUE

      RETURN
      // ..
      // .. END OF ZGSVJ0
      // ..
      }
