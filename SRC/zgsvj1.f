      SUBROUTINE ZGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE
      // .. Scalar Arguments ..
      double             EPS, SFMIN, TOL;
      int                INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK );
      double             SVA( N );
      // ..

*  =====================================================================

      // .. Local Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         AAPQ, OMPQ;
      double             AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN;
      int                BLSKIP, EMPTSW, i, ibr, igl, IERR, IJBLSK, ISWROT, jbc, jgl, KBL, MVL, NOTROT, nblc, nblr, p, PSKIPPED, q, ROWSKIP, SWBAND;
      bool               APPLV, ROTOK, RSVEC;
      // ..
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, DBLE, MIN, SIGN, SQRT
      // ..
      // .. External Functions ..
      double             DZNRM2;
      COMPLEX*16         ZDOTC;
      int                IDAMAX;
      bool               LSAME;
      // EXTERNAL IDAMAX, LSAME, ZDOTC, DZNRM2
      // ..
      // .. External Subroutines ..
      // .. from BLAS
      // EXTERNAL ZCOPY, ZROT, ZSWAP, ZAXPY
      // .. from LAPACK
      // EXTERNAL ZLASCL, ZLASSQ, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      APPLV = LSAME( JOBV, 'A' );
      RSVEC = LSAME( JOBV, 'V' );
      if ( !( RSVEC || APPLV || LSAME( JOBV, 'N' ) ) ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( ( N < 0 ) || ( N > M ) ) {
         INFO = -3;
      } else if ( N1 < 0 ) {
         INFO = -4;
      } else if ( LDA < M ) {
         INFO = -6;
      } else if ( ( RSVEC || APPLV ) && ( MV < 0 ) ) {
         INFO = -9;
      } else if ( ( RSVEC && ( LDV < N ) ) || ( APPLV && ( LDV < MV ) )  ) {
         INFO = -11;
      } else if ( TOL <= EPS ) {
         INFO = -14;
      } else if ( NSWEEP < 0 ) {
         INFO = -15;
      } else if ( LWORK < M ) {
         INFO = -17;
      } else {
         INFO = 0;
      }

      // #:(
      if ( INFO != 0 ) {
         xerbla('ZGSVJ1', -INFO );
         return;
      }

      if ( RSVEC ) {
         MVL = N;
      } else if ( APPLV ) {
         MVL = MV;
      }
      RSVEC = RSVEC || APPLV;

      ROOTEPS = SQRT( EPS );
      ROOTSFMIN = SQRT( SFMIN );
      SMALL = SFMIN / EPS;
      BIG = ONE / SFMIN;
      ROOTBIG = ONE / ROOTSFMIN;
      // LARGE = BIG / SQRT( DBLE( M*N ) )
      BIGTHETA = ONE / ROOTEPS;
      ROOTTOL = SQRT( TOL );

      // .. Initialize the right singular vector matrix ..

      // RSVEC = LSAME( JOBV, 'Y' )

      EMPTSW = N1*( N-N1 );
      NOTROT = 0;

      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

      KBL = MIN( 8, N );
      NBLR = N1 / KBL;
      IF( ( NBLR*KBL ) != N1 )NBLR = NBLR + 1;

      // .. the tiling is nblr-by-nblc [tiles]

      NBLC = ( N-N1 ) / KBL;
      IF( ( NBLC*KBL ) != ( N-N1 ) )NBLC = NBLC + 1;
      BLSKIP = ( KBL**2 ) + 1;
*[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = MIN( 5, KBL );
*[TP] ROWSKIP is a tuning parameter.
      SWBAND = 0;
*[TP] SWBAND is a tuning parameter. It is meaningful and effective
      // if ZGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm ZGEJSV.


      // | *   *   * [x] [x] [x]|
      // | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
      // | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
      // |[x] [x] [x] *   *   * |
      // |[x] [x] [x] *   *   * |
      // |[x] [x] [x] *   *   * |


      for (i = 1; i <= NSWEEP; i++) { // 1993

      // .. go go go ...

         MXAAPQ = ZERO;
         MXSINJ = ZERO;
         ISWROT = 0;

         NOTROT = 0;
         PSKIPPED = 0;

      // Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
      // 1 <= p < q <= N. This is the first step toward a blocked implementation
      // of the rotations. New implementation, based on block transformations,
      // is under development.

         for (ibr = 1; ibr <= NBLR; ibr++) { // 2000

            igl = ( ibr-1 )*KBL + 1;



* ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1;

             // DO 2010 jbc = ibr + 1, NBL
            for (jbc = 1; jbc <= NBLC; jbc++) { // 2010

               jgl = ( jbc-1 )*KBL + N1 + 1;

         // doing the block at ( ibr, jbc )

               IJBLSK = 0;
               DO 2100 p = igl, MIN( igl+KBL-1, N1 );

                  AAPP = SVA( p );
                  if ( AAPP > ZERO ) {

                     PSKIPPED = 0;

                     DO 2200 q = jgl, MIN( jgl+KBL-1, N );

                        AAQQ = SVA( q );
                        if ( AAQQ > ZERO ) {
                           AAPP0 = AAPP;

      // .. M x 2 Jacobi SVD ..

         // Safe Gram matrix computation

                           if ( AAQQ >= ONE ) {
                              if ( AAPP >= AAQQ ) {
                                 ROTOK = ( SMALL*AAPP ) <= AAQQ;
                              } else {
                                 ROTOK = ( SMALL*AAQQ ) <= AAPP;
                              }
                              if ( AAPP < ( BIG / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP;
                              } else {
                                 zcopy(M, A( 1, p ), 1, WORK, 1 );
                                 zlascl('G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR );
                                 AAPQ = ZDOTC( M, WORK, 1, A( 1, q ), 1 ) / AAQQ;
                              }
                           } else {
                              if ( AAPP >= AAQQ ) {
                                 ROTOK = AAPP <= ( AAQQ / SMALL );
                              } else {
                                 ROTOK = AAQQ <= ( AAPP / SMALL );
                              }
                              if ( AAPP > ( SMALL / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / MAX(AAQQ,AAPP) ) / MIN(AAQQ,AAPP);
                              } else {
                                 zcopy(M, A( 1, q ), 1, WORK, 1 );
                                 zlascl('G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR );
                                 AAPQ = ZDOTC( M, A( 1, p ), 1, WORK, 1 ) / AAPP;
                              }
                           }

                            // AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                           AAPQ1  = -ABS(AAPQ);
                           MXAAPQ = MAX( MXAAPQ, -AAPQ1 );

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ABS( AAPQ1 ) > TOL ) {
                              OMPQ = AAPQ / ABS(AAPQ);
                              NOTROT = 0;
*[RTD]      ROTATED  = ROTATED + 1
                              PSKIPPED = 0;
                              ISWROT = ISWROT + 1;

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP;
                                 APOAQ = AAPP / AAQQ;
                                 THETA = -HALF*ABS( AQOAP-APOAQ )/ AAPQ1;
                                 if (AAQQ > AAPP0) THETA = -THETA;

                                 if ( ABS( THETA ) > BIGTHETA ) {
                                    T  = HALF / THETA;
                                    CS = ONE;
                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
                                    }
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) );
                                    MXSINJ = MAX( MXSINJ, ABS( T ) );
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ1 );
                                    if (AAQQ > AAPP0) THSIGN = -THSIGN;
                                    T = ONE / ( THETA+THSIGN* SQRT( ONE+THETA*THETA ) );
                                    CS = SQRT( ONE / ( ONE+T*T ) );
                                    SN = T*CS;
                                    MXSINJ = MAX( MXSINJ, ABS( SN ) );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) );

                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 D(p) = -D(q) * OMPQ;

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                               if ( AAPP > AAQQ ) {
                                    zcopy(M, A( 1, p ), 1, WORK, 1 );
                                    zlascl('G', 0, 0, AAPP, ONE, M, 1, WORK,LDA, IERR );
                                    zlascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                    zaxpy(M, -AAPQ, WORK, 1, A( 1, q ), 1 );
                                    zlascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                    MXSINJ = MAX( MXSINJ, SFMIN );
                               } else {
                                   zcopy(M, A( 1, q ), 1, WORK, 1 );
                                    zlascl('G', 0, 0, AAQQ, ONE, M, 1, WORK,LDA, IERR );
                                    zlascl('G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR );
                                    zaxpy(M, -CONJG(AAPQ), WORK, 1, A( 1, p ), 1 );
                                    zlascl('G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA( p ) = AAPP*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                    MXSINJ = MAX( MXSINJ, SFMIN );
                               }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // .. recompute SVA(q), SVA(p)
                              if ( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) THEN                                  IF( ( AAQQ < ROOTBIG ) && ( AAQQ > ROOTSFMIN ) ) {
                                    SVA( q ) = DZNRM2( M, A( 1, q ), 1);
                                  } else {
                                    T = ZERO;
                                    AAQQ = ONE;
                                    zlassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ );
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2 <= ROOTEPS ) {
                                 if ( ( AAPP < ROOTBIG ) && ( AAPP > ROOTSFMIN ) ) {
                                    AAPP = DZNRM2( M, A( 1, p ), 1 );
                                 } else {
                                    T = ZERO;
                                    AAPP = ONE;
                                    zlassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*SQRT( AAPP );
                                 }
                                 SVA( p ) = AAPP;
                              }
               // end of OK rotation
                           } else {
                              NOTROT = NOTROT + 1;
*[RTD]      SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1;
                              IJBLSK = IJBLSK + 1;
                           }
                        } else {
                           NOTROT = NOTROT + 1;
                           PSKIPPED = PSKIPPED + 1;
                           IJBLSK = IJBLSK + 1;
                        }

                        if ( ( i <= SWBAND ) && ( IJBLSK >= BLSKIP ) ) {
                           SVA( p ) = AAPP;
                           NOTROT = 0;
                           GO TO 2011;
                        }
                        if ( ( i <= SWBAND ) && ( PSKIPPED > ROWSKIP ) ) {
                           AAPP = -AAPP;
                           NOTROT = 0;
                           GO TO 2203;
                        }

                     } // 2200
         // end of the q-loop
                     } // 2203

                     SVA( p ) = AAPP;

                  } else {

                     if (AAPP == ZERO) NOTROT = NOTROT + MIN( jgl+KBL-1, N ) - jgl + 1;
                     if (AAPP < ZERO) NOTROT = 0;

                  }

               } // 2100
      // end of the p-loop
            } // 2010
      // end of the jbc-loop
            } // 2011
*2011 bailed out of the jbc-loop
            DO 2012 p = igl, MIN( igl+KBL-1, N );
               SVA( p ) = ABS( SVA( p ) );
            } // 2012
***
         } // 2000
*2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ) < ROOTBIG ) && ( SVA( N ) > ROOTSFMIN ) ) {
            SVA( N ) = DZNRM2( M, A( 1, N ), 1 );
         } else {
            T = ZERO;
            AAPP = ONE;
            zlassq(M, A( 1, N ), 1, T, AAPP );
            SVA( N ) = T*SQRT( AAPP );
         }

      // Additional steering devices

         IF( ( i < SWBAND ) && ( ( MXAAPQ <= ROOTTOL ) || ( ISWROT <= N ) ) )SWBAND = i;

         if ( ( i > SWBAND+1 ) && ( MXAAPQ < SQRT( DBLE( N ) )* TOL ) && ( DBLE( N )*MXAAPQ*MXSINJ < TOL ) ) {
            GO TO 1994;
         }

         if (NOTROT >= EMPTSW) GO TO 1994;

      } // 1993
      // end i=1:NSWEEP loop

* #:( Reaching this point means that the procedure has not converged.
      INFO = NSWEEP - 1;
      GO TO 1995;

      } // 1994
* #:) Reaching this point means numerical convergence after the i-th
      // sweep.

      INFO = 0;
* #:) INFO = 0 confirms successful iterations.
      } // 1995

      // Sort the vector SVA() of column norms.
      for (p = 1; p <= N - 1; p++) { // 5991
         q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1;
         if ( p != q ) {
            TEMP1 = SVA( p );
            SVA( p ) = SVA( q );
            SVA( q ) = TEMP1;
            AAPQ = D( p );
            D( p ) = D( q );
            D( q ) = AAPQ;
            zswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            if (RSVEC) CALL ZSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 );
         }
      } // 5991


      return;
      // ..
      // .. END OF ZGSVJ1
      // ..
      }
