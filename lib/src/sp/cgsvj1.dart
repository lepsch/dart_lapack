      void cgsvj1(JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               EPS, SFMIN, TOL;
      int                INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK );
      REAL               SVA( N );
      // ..

// =====================================================================

      // .. Local Parameters ..
      REAL               ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0;
      // ..
      // .. Local Scalars ..
      Complex            AAPQ, OMPQ;
      REAL               AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN;
      int                BLSKIP, EMPTSW, i, ibr, igl, IERR, IJBLSK, ISWROT, jbc, jgl, KBL, MVL, NOTROT, nblc, nblr, p, PSKIPPED, q, ROWSKIP, SWBAND;
      bool               APPLV, ROTOK, RSVEC;
      // ..
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, CONJG, REAL, MIN, SIGN, SQRT
      // ..
      // .. External Functions ..
      //- REAL               SCNRM2;
      //- COMPLEX            CDOTC;
      //- int                ISAMAX;
      //- bool               lsame;
      // EXTERNAL ISAMAX, lsame, CDOTC, SCNRM2
      // ..
      // .. External Subroutines ..
      // .. from BLAS
      // EXTERNAL CCOPY, CROT, CSWAP, CAXPY
      // .. from LAPACK
      // EXTERNAL CLASCL, CLASSQ, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      APPLV = lsame( JOBV, 'A' );
      RSVEC = lsame( JOBV, 'V' );
      if ( !( RSVEC || APPLV || lsame( JOBV, 'N' ) ) ) {
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
         xerbla('CGSVJ1', -INFO );
         return;
      }

      if ( RSVEC ) {
         MVL = N;
      } else if ( APPLV ) {
         MVL = MV;
      }
      RSVEC = RSVEC || APPLV;

      ROOTEPS = sqrt( EPS );
      ROOTSFMIN = sqrt( SFMIN );
      SMALL = SFMIN / EPS;
      BIG = ONE / SFMIN;
      ROOTBIG = ONE / ROOTSFMIN;
      // LARGE = BIG / sqrt( REAL( M*N ) )
      BIGTHETA = ONE / ROOTEPS;
      ROOTTOL = sqrt( TOL );

      // .. Initialize the right singular vector matrix ..

      // RSVEC = lsame( JOBV, 'Y' )

      EMPTSW = N1*( N-N1 );
      NOTROT = 0;

      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

      KBL = min( 8, N );
      NBLR = N1 / KBL;
      if( ( NBLR*KBL ) != N1 )NBLR = NBLR + 1;

      // .. the tiling is nblr-by-nblc [tiles]

      NBLC = ( N-N1 ) / KBL;
      if( ( NBLC*KBL ) != ( N-N1 ) )NBLC = NBLC + 1;
      BLSKIP = ( KBL**2 ) + 1;
// [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = min( 5, KBL );
// [TP] ROWSKIP is a tuning parameter.
      SWBAND = 0;
// [TP] SWBAND is a tuning parameter. It is meaningful and effective
      // if CGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm CGEJSV.


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



// ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1;

             // DO 2010 jbc = ibr + 1, NBL
            for (jbc = 1; jbc <= NBLC; jbc++) { // 2010

               jgl = ( jbc-1 )*KBL + N1 + 1;

         // doing the block at ( ibr, jbc )

               IJBLSK = 0;
               for (p = igl; p <= min( igl+KBL-1, N1 ); p++) { // 2100

                  AAPP = SVA( p );
                  if ( AAPP > ZERO ) {

                     PSKIPPED = 0;

                     for (q = jgl; q <= min( jgl+KBL-1, N ); q++) { // 2200

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
                                 AAPQ = ( CDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP;
                              } else {
                                 ccopy(M, A( 1, p ), 1, WORK, 1 );
                                 clascl('G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR );
                                 AAPQ = CDOTC( M, WORK, 1, A( 1, q ), 1 ) / AAQQ;
                              }
                           } else {
                              if ( AAPP >= AAQQ ) {
                                 ROTOK = AAPP <= ( AAQQ / SMALL );
                              } else {
                                 ROTOK = AAQQ <= ( AAPP / SMALL );
                              }
                              if ( AAPP > ( SMALL / AAQQ ) ) {
                                 AAPQ = ( CDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / max(AAQQ,AAPP) ) / min(AAQQ,AAPP);
                              } else {
                                 ccopy(M, A( 1, q ), 1, WORK, 1 );
                                 clascl('G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR );
                                 AAPQ = CDOTC( M, A( 1, p ), 1, WORK, 1 ) / AAPP;
                              }
                           }

                            // AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                           AAPQ1  = -(AAPQ).abs();
                           MXAAPQ = max( MXAAPQ, -AAPQ1 );

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ( AAPQ1 ).abs() > TOL ) {
                              OMPQ = AAPQ / (AAPQ).abs();
                              NOTROT = 0;
// [RTD]      ROTATED  = ROTATED + 1
                              PSKIPPED = 0;
                              ISWROT = ISWROT + 1;

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP;
                                 APOAQ = AAPP / AAQQ;
                                 THETA = -HALF*( AQOAP-APOAQ ).abs()/ AAPQ1;
                                 if (AAQQ > AAPP0) THETA = -THETA;

                                 if ( ( THETA ).abs() > BIGTHETA ) {
                                    T  = HALF / THETA;
                                    CS = ONE;
                                    crot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        crot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
                                    }
                                    SVA[q] = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ1 ) );
                                    MXSINJ = max( MXSINJ, ( T ).abs() );
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ1 );
                                    if (AAQQ > AAPP0) THSIGN = -THSIGN;
                                    T = ONE / ( THETA+THSIGN* sqrt( ONE+THETA*THETA ) );
                                    CS = sqrt( ONE / ( ONE+T*T ) );
                                    SN = T*CS;
                                    MXSINJ = max( MXSINJ, ( SN ).abs() );
                                    SVA[q] = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ1 ) );

                                    crot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        crot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 D[p] = -D(q) * OMPQ;

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                               if ( AAPP > AAQQ ) {
                                    ccopy(M, A( 1, p ), 1, WORK, 1 );
                                    clascl('G', 0, 0, AAPP, ONE, M, 1, WORK,LDA, IERR );
                                    clascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                    caxpy(M, -AAPQ, WORK, 1, A( 1, q ), 1 );
                                    clascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA[q] = AAQQ*sqrt( max( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                               } else {
                                   ccopy(M, A( 1, q ), 1, WORK, 1 );
                                    clascl('G', 0, 0, AAQQ, ONE, M, 1, WORK,LDA, IERR );
                                    clascl('G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR );
                                    caxpy(M, -CONJG(AAPQ), WORK, 1, A( 1, p ), 1 );
                                    clascl('G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA[p] = AAPP*sqrt( max( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                               }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // .. recompute SVA(q), SVA(p)
                              if ( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) {
                                 if( ( AAQQ < ROOTBIG ) && ( AAQQ > ROOTSFMIN ) ) {
                                    SVA[q] = SCNRM2( M, A( 1, q ), 1);
                                  } else {
                                    T = ZERO;
                                    AAQQ = ONE;
                                    classq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA[q] = T*sqrt( AAQQ );
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2 <= ROOTEPS ) {
                                 if ( ( AAPP < ROOTBIG ) && ( AAPP > ROOTSFMIN ) ) {
                                    AAPP = SCNRM2( M, A( 1, p ), 1 );
                                 } else {
                                    T = ZERO;
                                    AAPP = ONE;
                                    classq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*sqrt( AAPP );
                                 }
                                 SVA[p] = AAPP;
                              }
               // end of OK rotation
                           } else {
                              NOTROT = NOTROT + 1;
// [RTD]      SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1;
                              IJBLSK = IJBLSK + 1;
                           }
                        } else {
                           NOTROT = NOTROT + 1;
                           PSKIPPED = PSKIPPED + 1;
                           IJBLSK = IJBLSK + 1;
                        }

                        if ( ( i <= SWBAND ) && ( IJBLSK >= BLSKIP ) ) {
                           SVA[p] = AAPP;
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

                     SVA[p] = AAPP;

                  } else {

                     if (AAPP == ZERO) NOTROT = NOTROT + min( jgl+KBL-1, N ) - jgl + 1;
                     if (AAPP < ZERO) NOTROT = 0;

                  }

               } // 2100
      // end of the p-loop
            } // 2010
      // end of the jbc-loop
            } // 2011
// 2011 bailed out of the jbc-loop
            for (p = igl; p <= min( igl+KBL-1, N ); p++) { // 2012
               SVA[p] = ( SVA( p ) ).abs();
            } // 2012
// **
         } // 2000
// 2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ) < ROOTBIG ) && ( SVA( N ) > ROOTSFMIN ) ) {
            SVA[N] = SCNRM2( M, A( 1, N ), 1 );
         } else {
            T = ZERO;
            AAPP = ONE;
            classq(M, A( 1, N ), 1, T, AAPP );
            SVA[N] = T*sqrt( AAPP );
         }

      // Additional steering devices

         if( ( i < SWBAND ) && ( ( MXAAPQ <= ROOTTOL ) || ( ISWROT <= N ) ) )SWBAND = i;

         if ( ( i > SWBAND+1 ) && ( MXAAPQ < sqrt( REAL( N ) )* TOL ) && ( REAL( N )*MXAAPQ*MXSINJ < TOL ) ) {
            GO TO 1994;
         }

         if (NOTROT >= EMPTSW) GO TO 1994;

      } // 1993
      // end i=1:NSWEEP loop

// #:( Reaching this point means that the procedure has not converged.
      INFO = NSWEEP - 1;
      GO TO 1995;

      } // 1994
// #:) Reaching this point means numerical convergence after the i-th
      // sweep.

      INFO = 0;
// #:) INFO = 0 confirms successful iterations.
      } // 1995

      // Sort the vector SVA() of column norms.
      for (p = 1; p <= N - 1; p++) { // 5991
         q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1;
         if ( p != q ) {
            TEMP1 = SVA( p );
            SVA[p] = SVA( q );
            SVA[q] = TEMP1;
            AAPQ = D( p );
            D[p] = D( q );
            D[q] = AAPQ;
            cswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            if (RSVEC) cswap( MVL, V( 1, p ), 1, V( 1, q ), 1 );
         }
      } // 5991


      return;
      // ..
      // .. END OF CGSVJ1
      // ..
      }
