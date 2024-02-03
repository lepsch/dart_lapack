      void sgsvj0(JOBV, M, N, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDV, LWORK, M, MV, N, NSWEEP;
      REAL               EPS, SFMIN, TOL;
      String             JOBV;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), SVA( N ), D( N ), V( LDV, * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Local Parameters ..
      REAL               ZERO, HALF, ONE;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0;
      // ..
      // .. Local Scalars ..
      REAL               AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN;
      int                BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND;
      bool               APPLV, ROTOK, RSVEC;
      // ..
      // .. Local Arrays ..
      REAL               FASTR( 5 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, FLOAT, MIN, SIGN, SQRT
      // ..
      // .. External Functions ..
      //- REAL               SDOT, SNRM2;
      //- int                ISAMAX;
      //- bool               LSAME;
      // EXTERNAL ISAMAX, LSAME, SDOT, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SLASCL, SLASSQ, SROTM, SSWAP, XERBLA
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
      } else if ( LDA < M ) {
         INFO = -5;
      } else if ( ( RSVEC || APPLV ) && ( MV < 0 ) ) {
         INFO = -8;
      } else if ( ( RSVEC && ( LDV < N ) ) || ( APPLV && ( LDV < MV ) ) ) {
         INFO = -10;
      } else if ( TOL <= EPS ) {
         INFO = -13;
      } else if ( NSWEEP < 0 ) {
         INFO = -14;
      } else if ( LWORK < M ) {
         INFO = -16;
      } else {
         INFO = 0;
      }

      // #:(
      if ( INFO != 0 ) {
         xerbla('SGSVJ0', -INFO );
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
      BIGTHETA = ONE / ROOTEPS;
      ROOTTOL = sqrt( TOL );

      // .. Row-cyclic Jacobi SVD algorithm with column pivoting ..

      EMPTSW = ( N*( N-1 ) ) / 2;
      NOTROT = 0;
      FASTR( 1 ) = ZERO;

      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..


      SWBAND = 0;
// [TP] SWBAND is a tuning parameter. It is meaningful and effective
      // if SGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
      // ......

      KBL = min( 8, N );
// [TP] KBL is a tuning parameter that defines the tile size in the
      // tiling of the p-q loops of pivot pairs. In general, an optimal
      // value of KBL depends on the matrix dimensions and on the
      // parameters of the computer's memory.

      NBL = N / KBL;
      if( ( NBL*KBL ) != N )NBL = NBL + 1;

      BLSKIP = ( KBL**2 ) + 1;
// [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = min( 5, KBL );
// [TP] ROWSKIP is a tuning parameter.

      LKAHEAD = 1;
// [TP] LKAHEAD is a tuning parameter.
      SWBAND = 0;
      PSKIPPED = 0;

      for (i = 1; i <= NSWEEP; i++) { // 1993
      // .. go go go ...

         MXAAPQ = ZERO;
         MXSINJ = ZERO;
         ISWROT = 0;

         NOTROT = 0;
         PSKIPPED = 0;

         for (ibr = 1; ibr <= NBL; ibr++) { // 2000

            igl = ( ibr-1 )*KBL + 1;

            DO 1002 ir1 = 0, min( LKAHEAD, NBL-ibr );

               igl = igl + ir1*KBL;

               DO 2001 p = igl, min( igl+KBL-1, N-1 );

      // .. de Rijk's pivoting
                  q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1;
                  if ( p != q ) {
                     sswap(M, A( 1, p ), 1, A( 1, q ), 1 );
                     if (RSVEC) sswap( MVL, V( 1, p ), 1, V( 1, q ), 1 );
                     TEMP1 = SVA( p );
                     SVA( p ) = SVA( q );
                     SVA( q ) = TEMP1;
                     TEMP1 = D( p );
                     D( p ) = D( q );
                     D( q ) = TEMP1;
                  }

                  if ( ir1 == 0 ) {

         // Column norms are periodically updated by explicit
         // norm computation.
         // Caveat:
         // Some BLAS implementations compute SNRM2(M,A(1,p),1)
         // as sqrt(SDOT(M,A(1,p),1,A(1,p),1)), which may result in
         // overflow for ||A(:,p)||_2 > sqrt(overflow_threshold), and
         // underflow for ||A(:,p)||_2 < sqrt(underflow_threshold).
         // Hence, SNRM2 cannot be trusted, not even in the case when
         // the true norm is far from the under(over)flow boundaries.
         // If properly implemented SNRM2 is available, the IF-THEN-ELSE
         // below should read "AAPP = SNRM2( M, A(1,p), 1 ) * D(p)".

                     if ( ( SVA( p ) < ROOTBIG ) && ( SVA( p ) > ROOTSFMIN ) ) {
                        SVA( p ) = SNRM2( M, A( 1, p ), 1 )*D( p );
                     } else {
                        TEMP1 = ZERO;
                        AAPP = ONE;
                        slassq(M, A( 1, p ), 1, TEMP1, AAPP );
                        SVA( p ) = TEMP1*sqrt( AAPP )*D( p );
                     }
                     AAPP = SVA( p );
                  } else {
                     AAPP = SVA( p );
                  }


                  if ( AAPP > ZERO ) {

                     PSKIPPED = 0;

                     DO 2002 q = p + 1, min( igl+KBL-1, N );

                        AAQQ = SVA( q );

                        if ( AAQQ > ZERO ) {

                           AAPP0 = AAPP;
                           if ( AAQQ >= ONE ) {
                              ROTOK = ( SMALL*AAPP ) <= AAQQ;
                              if ( AAPP < ( BIG / AAQQ ) ) {
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP;
                              } else {
                                 scopy(M, A( 1, p ), 1, WORK, 1 );
                                 CALL SLASCL( 'G', 0, 0, AAPP, D( p ), M, 1, WORK, LDA, IERR )                                  AAPQ = SDOT( M, WORK, 1, A( 1, q ), 1 )*D( q ) / AAQQ;
                              }
                           } else {
                              ROTOK = AAPP <= ( AAQQ / SMALL );
                              if ( AAPP > ( SMALL / AAQQ ) ) {
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP;
                              } else {
                                 scopy(M, A( 1, q ), 1, WORK, 1 );
                                 CALL SLASCL( 'G', 0, 0, AAQQ, D( q ), M, 1, WORK, LDA, IERR )                                  AAPQ = SDOT( M, WORK, 1, A( 1, p ), 1 )*D( p ) / AAPP;
                              }
                           }

                           MXAAPQ = max( MXAAPQ, ( AAPQ ).abs() );

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ( AAPQ ).abs() > TOL ) {

            // .. rotate
            // ROTATED = ROTATED + ONE

                              if ( ir1 == 0 ) {
                                 NOTROT = 0;
                                 PSKIPPED = 0;
                                 ISWROT = ISWROT + 1;
                              }

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP;
                                 APOAQ = AAPP / AAQQ;
                                 THETA = -HALF*( AQOAP-APOAQ ).abs() / AAPQ;

                                 if ( ( THETA ).abs() > BIGTHETA ) {

                                    T = HALF / THETA;
                                    FASTR( 3 ) = T*D( p ) / D( q );
                                    FASTR( 4 ) = -T*D( q ) / D( p );
                                    srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                    SVA( q ) = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ ) );
                                    MXSINJ = max( MXSINJ, ( T ).abs() );

                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ );
                                    T = ONE / ( THETA+THSIGN* sqrt( ONE+THETA*THETA ) );
                                    CS = sqrt( ONE / ( ONE+T*T ) );
                                    SN = T*CS;

                                    MXSINJ = max( MXSINJ, ( SN ).abs() );
                                    SVA( q ) = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ ) );

                                    APOAQ = D( p ) / D( q );
                                    AQOAP = D( q ) / D( p );
                                    if ( D( p ) >= ONE ) {
                                       if ( D( q ) >= ONE ) {
                                          FASTR( 3 ) = T*APOAQ;
                                          FASTR( 4 ) = -T*AQOAP;
                                          D( p ) = D( p )*CS;
                                          D( q ) = D( q )*CS;
                                          srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                       } else {
                                          saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          D( p ) = D( p )*CS;
                                          D( q ) = D( q ) / CS;
                                          if ( RSVEC ) {
                                             saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                          }
                                       }
                                    } else {
                                       if ( D( q ) >= ONE ) {
                                          saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          D( p ) = D( p ) / CS;
                                          D( q ) = D( q )*CS;
                                          if ( RSVEC ) {
                                             saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             saxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                          }
                                       } else {
                                          if ( D( p ) >= D( q ) ) {
                                             saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             D( p ) = D( p )*CS;
                                             D( q ) = D( q ) / CS;
                                             if ( RSVEC ) {
                                                saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                                saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             }
                                          } else {
                                             saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             D( p ) = D( p ) / CS;
                                             D( q ) = D( q )*CS;
                                             if ( RSVEC ) {
                                                saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                                saxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             }
                                          }
                                       }
                                    }
                                 }

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                                 scopy(M, A( 1, p ), 1, WORK, 1 );
                                 slascl('G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR );
                                 slascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                 TEMP1 = -AAPQ*D( p ) / D( q );
                                 saxpy(M, TEMP1, WORK, 1, A( 1, q ), 1 );
                                 slascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )                                  SVA( q ) = AAQQ*sqrt( max( ZERO, ONE-AAPQ*AAPQ ) );
                                 MXSINJ = max( MXSINJ, SFMIN );
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // recompute SVA(q), SVA(p).
                              if( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) {
                                 if( ( AAQQ < ROOTBIG ) && ( AAQQ > ROOTSFMIN ) ) {
                                    SVA( q ) = SNRM2( M, A( 1, q ), 1 )* D( q );
                                 } else {
                                    T = ZERO;
                                    AAQQ = ONE;
                                    slassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*sqrt( AAQQ )*D( q );
                                 }
                              }
                              if ( ( AAPP / AAPP0 ) <= ROOTEPS ) {
                                 if( ( AAPP < ROOTBIG ) && ( AAPP > ROOTSFMIN ) ) {
                                    AAPP = SNRM2( M, A( 1, p ), 1 )* D( p );
                                 } else {
                                    T = ZERO;
                                    AAPP = ONE;
                                    slassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*sqrt( AAPP )*D( p );
                                 }
                                 SVA( p ) = AAPP;
                              }

                           } else {
         // A(:,p) and A(:,q) already numerically orthogonal
                              if (ir1 == 0) NOTROT = NOTROT + 1;
                              PSKIPPED = PSKIPPED + 1;
                           }
                        } else {
         // A(:,q) is zero column
                           if (ir1 == 0) NOTROT = NOTROT + 1;
                           PSKIPPED = PSKIPPED + 1;
                        }

                        if ( ( i <= SWBAND ) && ( PSKIPPED > ROWSKIP ) ) {
                           if (ir1 == 0) AAPP = -AAPP;
                           NOTROT = 0;
                           GO TO 2103;
                        }

                     } // 2002
      // END q-LOOP

                     } // 2103
      // bailed out of q-loop

                     SVA( p ) = AAPP;

                  } else {
                     SVA( p ) = AAPP;
                     if( ( ir1 == 0 ) && ( AAPP == ZERO ) ) NOTROT = NOTROT + min( igl+KBL-1, N ) - p;
                  }

               } // 2001
      // end of the p-loop
      // end of doing the block ( ibr, ibr )
            } // 1002
      // end of ir1-loop

// ........................................................
// ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1;

            for (jbc = ibr + 1; jbc <= NBL; jbc++) { // 2010

               jgl = ( jbc-1 )*KBL + 1;

         // doing the block at ( ibr, jbc )

               IJBLSK = 0;
               DO 2100 p = igl, min( igl+KBL-1, N );

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
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP;
                              } else {
                                 scopy(M, A( 1, p ), 1, WORK, 1 );
                                 CALL SLASCL( 'G', 0, 0, AAPP, D( p ), M, 1, WORK, LDA, IERR )                                  AAPQ = SDOT( M, WORK, 1, A( 1, q ), 1 )*D( q ) / AAQQ;
                              }
                           } else {
                              if ( AAPP >= AAQQ ) {
                                 ROTOK = AAPP <= ( AAQQ / SMALL );
                              } else {
                                 ROTOK = AAQQ <= ( AAPP / SMALL );
                              }
                              if ( AAPP > ( SMALL / AAQQ ) ) {
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*D( p )*D( q ) / AAQQ ) / AAPP;
                              } else {
                                 scopy(M, A( 1, q ), 1, WORK, 1 );
                                 CALL SLASCL( 'G', 0, 0, AAQQ, D( q ), M, 1, WORK, LDA, IERR )                                  AAPQ = SDOT( M, WORK, 1, A( 1, p ), 1 )*D( p ) / AAPP;
                              }
                           }

                           MXAAPQ = max( MXAAPQ, ( AAPQ ).abs() );

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ( AAPQ ).abs() > TOL ) {
                              NOTROT = 0;
            // ROTATED  = ROTATED + 1
                              PSKIPPED = 0;
                              ISWROT = ISWROT + 1;

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP;
                                 APOAQ = AAPP / AAQQ;
                                 THETA = -HALF*( AQOAP-APOAQ ).abs() / AAPQ;
                                 if (AAQQ > AAPP0) THETA = -THETA;

                                 if ( ( THETA ).abs() > BIGTHETA ) {
                                    T = HALF / THETA;
                                    FASTR( 3 ) = T*D( p ) / D( q );
                                    FASTR( 4 ) = -T*D( q ) / D( p );
                                    srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                    SVA( q ) = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ ) );
                                    MXSINJ = max( MXSINJ, ( T ).abs() );
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ );
                                    if (AAQQ > AAPP0) THSIGN = -THSIGN;
                                    T = ONE / ( THETA+THSIGN* sqrt( ONE+THETA*THETA ) );
                                    CS = sqrt( ONE / ( ONE+T*T ) );
                                    SN = T*CS;
                                    MXSINJ = max( MXSINJ, ( SN ).abs() );
                                    SVA( q ) = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ ) );

                                    APOAQ = D( p ) / D( q );
                                    AQOAP = D( q ) / D( p );
                                    if ( D( p ) >= ONE ) {

                                       if ( D( q ) >= ONE ) {
                                          FASTR( 3 ) = T*APOAQ;
                                          FASTR( 4 ) = -T*AQOAP;
                                          D( p ) = D( p )*CS;
                                          D( q ) = D( q )*CS;
                                          srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                       } else {
                                          saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          if ( RSVEC ) {
                                             saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                          }
                                          D( p ) = D( p )*CS;
                                          D( q ) = D( q ) / CS;
                                       }
                                    } else {
                                       if ( D( q ) >= ONE ) {
                                          saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          if ( RSVEC ) {
                                             saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             saxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                          }
                                          D( p ) = D( p ) / CS;
                                          D( q ) = D( q )*CS;
                                       } else {
                                          if ( D( p ) >= D( q ) ) {
                                             saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             D( p ) = D( p )*CS;
                                             D( q ) = D( q ) / CS;
                                             if ( RSVEC ) {
                                                saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                                saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             }
                                          } else {
                                             saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             D( p ) = D( p ) / CS;
                                             D( q ) = D( q )*CS;
                                             if ( RSVEC ) {
                                                saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                                saxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             }
                                          }
                                       }
                                    }
                                 }

                              } else {
                                 if ( AAPP > AAQQ ) {
                                    scopy(M, A( 1, p ), 1, WORK, 1 );
                                    slascl('G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR );
                                    slascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                    TEMP1 = -AAPQ*D( p ) / D( q );
                                    saxpy(M, TEMP1, WORK, 1, A( 1, q ), 1 );
                                    slascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA( q ) = AAQQ*sqrt( max( ZERO, ONE-AAPQ*AAPQ ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                                 } else {
                                    scopy(M, A( 1, q ), 1, WORK, 1 );
                                    slascl('G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR );
                                    slascl('G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR );
                                    TEMP1 = -AAPQ*D( q ) / D( p );
                                    saxpy(M, TEMP1, WORK, 1, A( 1, p ), 1 );
                                    slascl('G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA( p ) = AAPP*sqrt( max( ZERO, ONE-AAPQ*AAPQ ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                                 }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q)
            // .. recompute SVA(q)
                              if( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) {
                                 if( ( AAQQ < ROOTBIG ) && ( AAQQ > ROOTSFMIN ) ) {
                                    SVA( q ) = SNRM2( M, A( 1, q ), 1 )* D( q );
                                 } else {
                                    T = ZERO;
                                    AAQQ = ONE;
                                    slassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*sqrt( AAQQ )*D( q );
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2 <= ROOTEPS ) {
                                 if( ( AAPP < ROOTBIG ) && ( AAPP > ROOTSFMIN ) ) {
                                    AAPP = SNRM2( M, A( 1, p ), 1 )* D( p );
                                 } else {
                                    T = ZERO;
                                    AAPP = ONE;
                                    slassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*sqrt( AAPP )*D( p );
                                 }
                                 SVA( p ) = AAPP;
                              }
               // end of OK rotation
                           } else {
                              NOTROT = NOTROT + 1;
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
                     if (AAPP == ZERO) NOTROT = NOTROT + min( jgl+KBL-1, N ) - jgl + 1;
                     if (AAPP < ZERO) NOTROT = 0;
                  }

               } // 2100
      // end of the p-loop
            } // 2010
      // end of the jbc-loop
            } // 2011
// 2011 bailed out of the jbc-loop
            DO 2012 p = igl, min( igl+KBL-1, N );
               SVA( p ) = ( SVA( p ) ).abs();
            } // 2012

         } // 2000
// 2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ) < ROOTBIG ) && ( SVA( N ) > ROOTSFMIN ) ) {
            SVA( N ) = SNRM2( M, A( 1, N ), 1 )*D( N );
         } else {
            T = ZERO;
            AAPP = ONE;
            slassq(M, A( 1, N ), 1, T, AAPP );
            SVA( N ) = T*sqrt( AAPP )*D( N );
         }

      // Additional steering devices

         if( ( i < SWBAND ) && ( ( MXAAPQ <= ROOTTOL ) || ( ISWROT <= N ) ) )SWBAND = i;

         if ( ( i > SWBAND+1 ) && ( MXAAPQ < FLOAT( N )*TOL ) && ( FLOAT( N )*MXAAPQ*MXSINJ < TOL ) ) {
            GO TO 1994;
         }

         if (NOTROT >= EMPTSW) GO TO 1994;

      } // 1993
      // end i=1:NSWEEP loop
// #:) Reaching this point means that the procedure has completed the given
      // number of iterations.
      INFO = NSWEEP - 1;
      GO TO 1995;
      } // 1994
// #:) Reaching this point means that during the i-th sweep all pivots were
      // below the given tolerance, causing early exit.

      INFO = 0;
// #:) INFO = 0 confirms successful iterations.
      } // 1995

      // Sort the vector D.
      for (p = 1; p <= N - 1; p++) { // 5991
         q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1;
         if ( p != q ) {
            TEMP1 = SVA( p );
            SVA( p ) = SVA( q );
            SVA( q ) = TEMP1;
            TEMP1 = D( p );
            D( p ) = D( q );
            D( q ) = TEMP1;
            sswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            if (RSVEC) sswap( MVL, V( 1, p ), 1, V( 1, q ), 1 );
         }
      } // 5991

      return;
      // ..
      // .. END OF SGSVJ0
      // ..
      }
