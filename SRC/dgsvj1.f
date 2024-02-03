      SUBROUTINE DGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             EPS, SFMIN, TOL;
      int                INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( N ), SVA( N ), V( LDV, * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Local Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0 ;
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
         xerbla('DGSVJ1', -INFO );
         return;
      }

      if ( RSVEC ) {
         MVL = N;
      } else if ( APPLV ) {
         MVL = MV;
      }
      RSVEC = RSVEC || APPLV;

      ROOTEPS = DSQRT( EPS );
      ROOTSFMIN = DSQRT( SFMIN );
      SMALL = SFMIN / EPS;
      BIG = ONE / SFMIN;
      ROOTBIG = ONE / ROOTSFMIN;
      LARGE = BIG / DSQRT( DBLE( M*N ) );
      BIGTHETA = ONE / ROOTEPS;
      ROOTTOL = DSQRT( TOL );

      // .. Initialize the right singular vector matrix ..

      // RSVEC = LSAME( JOBV, 'Y' )

      EMPTSW = N1*( N-N1 );
      NOTROT = 0;
      FASTR( 1 ) = ZERO;

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
      // if SGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm SGESVJ.


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

         for (ibr = 1; ibr <= NBLR; ibr++) { // 2000

            igl = ( ibr-1 )*KBL + 1;


// ........................................................
// ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1;

            for (jbc = 1; jbc <= NBLC; jbc++) { // 2010

               jgl = N1 + ( jbc-1 )*KBL + 1;

         // doing the block at ( ibr, jbc )

               IJBLSK = 0;
               DO 2100 p = igl, min( igl+KBL-1, N1 );

                  AAPP = SVA( p );

                  if ( AAPP > ZERO ) {

                     PSKIPPED = 0;

                     DO 2200 q = jgl, min( jgl+KBL-1, N );

                        AAQQ = SVA( q );

                        if ( AAQQ > ZERO ) {
                           AAPP0 = AAPP;

      // .. M x 2 Jacobi SVD ..

         // .. Safe Gram matrix computation ..

                           if ( AAQQ >= ONE ) {
                              if ( AAPP >= AAQQ ) {
                                 ROTOK = ( SMALL*AAPP ) <= AAQQ;
                              } else {
                                 ROTOK = ( SMALL*AAQQ ) <= AAPP;
                              }
                              if ( AAPP < ( BIG / AAQQ ) ) {
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP;
                              } else {
                                 dcopy(M, A( 1, p ), 1, WORK, 1 );
                                 CALL DLASCL( 'G', 0, 0, AAPP, D( p ), M, 1, WORK, LDA, IERR )                                  AAPQ = DDOT( M, WORK, 1, A( 1, q ), 1 )*D( q ) / AAQQ;
                              }
                           } else {
                              if ( AAPP >= AAQQ ) {
                                 ROTOK = AAPP <= ( AAQQ / SMALL );
                              } else {
                                 ROTOK = AAQQ <= ( AAPP / SMALL );
                              }
                              if ( AAPP > ( SMALL / AAQQ ) ) {
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP;
                              } else {
                                 dcopy(M, A( 1, q ), 1, WORK, 1 );
                                 CALL DLASCL( 'G', 0, 0, AAQQ, D( q ), M, 1, WORK, LDA, IERR )                                  AAPQ = DDOT( M, WORK, 1, A( 1, p ), 1 )*D( p ) / AAPP;
                              }
                           }

                           MXAAPQ = max( MXAAPQ, DABS( AAPQ ) );

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( DABS( AAPQ ) > TOL ) {
                              NOTROT = 0;
            // ROTATED  = ROTATED + 1
                              PSKIPPED = 0;
                              ISWROT = ISWROT + 1;

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP;
                                 APOAQ = AAPP / AAQQ;
                                 THETA = -HALF*DABS(AQOAP-APOAQ) / AAPQ;
                                 if (AAQQ > AAPP0) THETA = -THETA;

                                 if ( DABS( THETA ) > BIGTHETA ) {
                                    T = HALF / THETA;
                                    FASTR( 3 ) = T*D( p ) / D( q );
                                    FASTR( 4 ) = -T*D( q ) / D( p );
                                    drotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                    SVA( q ) = AAQQ*DSQRT( max( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( max( ZERO, ONE-T*AQOAP*AAPQ ) );
                                    MXSINJ = max( MXSINJ, DABS( T ) );
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -DSIGN( ONE, AAPQ );
                                    if (AAQQ > AAPP0) THSIGN = -THSIGN;
                                    T = ONE / ( THETA+THSIGN* DSQRT( ONE+THETA*THETA ) );
                                    CS = DSQRT( ONE / ( ONE+T*T ) );
                                    SN = T*CS;
                                    MXSINJ = max( MXSINJ, DABS( SN ) );
                                    SVA( q ) = AAQQ*DSQRT( max( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( max( ZERO, ONE-T*AQOAP*AAPQ ) );

                                    APOAQ = D( p ) / D( q );
                                    AQOAP = D( q ) / D( p );
                                    if ( D( p ) >= ONE ) {

                                       if ( D( q ) >= ONE ) {
                                          FASTR( 3 ) = T*APOAQ;
                                          FASTR( 4 ) = -T*AQOAP;
                                          D( p ) = D( p )*CS;
                                          D( q ) = D( q )*CS;
                                          drotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                       } else {
                                          daxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          daxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          if ( RSVEC ) {
                                             daxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             daxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                          }
                                          D( p ) = D( p )*CS;
                                          D( q ) = D( q ) / CS;
                                       }
                                    } else {
                                       if ( D( q ) >= ONE ) {
                                          daxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          daxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          if ( RSVEC ) {
                                             daxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             daxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                          }
                                          D( p ) = D( p ) / CS;
                                          D( q ) = D( q )*CS;
                                       } else {
                                          if ( D( p ) >= D( q ) ) {
                                             daxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             daxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             D( p ) = D( p )*CS;
                                             D( q ) = D( q ) / CS;
                                             if ( RSVEC ) {
                                                daxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                                daxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             }
                                          } else {
                                             daxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             daxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             D( p ) = D( p ) / CS;
                                             D( q ) = D( q )*CS;
                                             if ( RSVEC ) {
                                                daxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                                daxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             }
                                          }
                                       }
                                    }
                                 }

                              } else {
                                 if ( AAPP > AAQQ ) {
                                    dcopy(M, A( 1, p ), 1, WORK, 1 );
                                    dlascl('G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR );
                                    dlascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                    TEMP1 = -AAPQ*D( p ) / D( q );
                                    daxpy(M, TEMP1, WORK, 1, A( 1, q ), 1 );
                                    dlascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA( q ) = AAQQ*DSQRT( max( ZERO, ONE-AAPQ*AAPQ ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                                 } else {
                                    dcopy(M, A( 1, q ), 1, WORK, 1 );
                                    dlascl('G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR );
                                    dlascl('G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR );
                                    TEMP1 = -AAPQ*D( q ) / D( p );
                                    daxpy(M, TEMP1, WORK, 1, A( 1, p ), 1 );
                                    dlascl('G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA( p ) = AAPP*DSQRT( max( ZERO, ONE-AAPQ*AAPQ ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                                 }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q)
            // .. recompute SVA(q)
                              if( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) {
                                 if( ( AAQQ < ROOTBIG ) && ( AAQQ > ROOTSFMIN ) ) {
                                    SVA( q ) = DNRM2( M, A( 1, q ), 1 )* D( q );
                                 } else {
                                    T = ZERO;
                                    AAQQ = ONE;
                                    dlassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*DSQRT( AAQQ )*D( q );
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2 <= ROOTEPS ) {
                                 if( ( AAPP < ROOTBIG ) && ( AAPP > ROOTSFMIN ) ) {
                                    AAPP = DNRM2( M, A( 1, p ), 1 )* D( p );
                                 } else {
                                    T = ZERO;
                                    AAPP = ONE;
                                    dlassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*DSQRT( AAPP )*D( p );
                                 }
                                 SVA( p ) = AAPP;
                              }
               // end of OK rotation
                           } else {
                              NOTROT = NOTROT + 1;
            // SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1;
                              IJBLSK = IJBLSK + 1;
                           }
                        } else {
                           NOTROT = NOTROT + 1;
                           PSKIPPED = PSKIPPED + 1;
                           IJBLSK = IJBLSK + 1;
                        }

       // IF ( NOTROT >= EMPTSW )  GO TO 2011
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
                     if (AAPP == ZERO) NOTROT = NOTROT + min( jgl+KBL-1, N ) - jgl + 1;
                     if (AAPP < ZERO) NOTROT = 0;
// **      IF ( NOTROT >= EMPTSW )  GO TO 2011
                  }

               } // 2100
      // end of the p-loop
            } // 2010
      // end of the jbc-loop
            } // 2011
// 2011 bailed out of the jbc-loop
            DO 2012 p = igl, min( igl+KBL-1, N );
               SVA( p ) = DABS( SVA( p ) );
            } // 2012
// **   IF ( NOTROT >= EMPTSW ) GO TO 1994
         } // 2000
// 2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ) < ROOTBIG ) && ( SVA( N ) > ROOTSFMIN ) ) {
            SVA( N ) = DNRM2( M, A( 1, N ), 1 )*D( N );
         } else {
            T = ZERO;
            AAPP = ONE;
            dlassq(M, A( 1, N ), 1, T, AAPP );
            SVA( N ) = T*DSQRT( AAPP )*D( N );
         }

      // Additional steering devices

         if( ( i < SWBAND ) && ( ( MXAAPQ <= ROOTTOL ) || ( ISWROT <= N ) ) )SWBAND = i;
          if ( ( i > SWBAND+1 ) && ( MXAAPQ < DBLE( N )*TOL ) && ( DBLE( N )*MXAAPQ*MXSINJ < TOL ) ) {
            GO TO 1994;
         }


         if (NOTROT >= EMPTSW) GO TO 1994;

      } // 1993
      // end i=1:NSWEEP loop
// #:) Reaching this point means that the procedure has completed the given
      // number of sweeps.
      INFO = NSWEEP - 1;
      GO TO 1995;
      } // 1994
// #:) Reaching this point means that during the i-th sweep all pivots were
      // below the given threshold, causing early exit.

      INFO = 0;
// #:) INFO = 0 confirms successful iterations.
      } // 1995

      // Sort the vector D

      for (p = 1; p <= N - 1; p++) { // 5991
         q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1;
         if ( p != q ) {
            TEMP1 = SVA( p );
            SVA( p ) = SVA( q );
            SVA( q ) = TEMP1;
            TEMP1 = D( p );
            D( p ) = D( q );
            D( q ) = TEMP1;
            dswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            if (RSVEC) CALL DSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 );
         }
      } // 5991

      return;
      // ..
      // .. END OF DGSVJ1
      // ..
      }
