      void zgesvj(final int JOBA, final int JOBU, final int JOBV, final int M, final int N, final Matrix<double> A, final int LDA, final int SVA, final int MV, final Matrix<double> V, final int LDV, final int CWORK, final int LWORK, final Array<int> RWORK, final int LRWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDV, LWORK, LRWORK, M, MV, N;
      String             JOBA, JOBU, JOBV;
      Complex         A( LDA, * ),  V( LDV, * ), CWORK( LWORK );
      double             RWORK( LRWORK ), SVA( N );
      // ..

// =====================================================================

      // .. Local Parameters ..
      double             ZERO,         HALF,         ONE;
      const            ZERO = 0.0, HALF = 0.5, ONE = 1.0;
      Complex         CZERO,                  CONE;
      const            CZERO = (0.0, 0.0), CONE = (1.0, 0.0) ;
      int                NSWEEP;
      const            NSWEEP = 30 ;
      Complex         AAPQ, OMPQ;
      double             AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, CTOL, EPSLN, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SKL, SFMIN, SMALL, SN, T, TEMP1, THETA, THSIGN, TOL;
      int                BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, N2, N34, N4, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND, MINMN, LWMIN, LRWMIN;
      bool               APPLV, GOSCALE, LOWER, LQUERY, LSVEC, NOSCALE, ROTOK, RSVEC, UCTOL, UPPER;
      // ..
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, CONJG, DBLE, SIGN, SQRT
      // ..
      // .. External Functions ..
      // ..
      // from BLAS
      double             DZNRM2;
      Complex         ZDOTC;
      // EXTERNAL ZDOTC, DZNRM2
      int                idamax;
      // EXTERNAL idamax
      // from LAPACK
      double             DLAMCH;
      // EXTERNAL DLAMCH
      bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // ..
      // from BLAS
      // EXTERNAL ZCOPY, ZROT, ZDSCAL, ZSWAP, ZAXPY
      // from LAPACK
      // EXTERNAL DLASCL, ZLASCL, ZLASET, ZLASSQ, XERBLA
      // EXTERNAL ZGSVJ0, ZGSVJ1

      // Test the input arguments

      LSVEC = lsame( JOBU, 'U' ) || lsame( JOBU, 'F' );
      UCTOL = lsame( JOBU, 'C' );
      RSVEC = lsame( JOBV, 'V' ) || lsame( JOBV, 'J' );
      APPLV = lsame( JOBV, 'A' );
      UPPER = lsame( JOBA, 'U' );
      LOWER = lsame( JOBA, 'L' );

      MINMN = min( M, N );
      if ( MINMN == 0 ) {
         LWMIN  = 1;
         LRWMIN = 1;
      } else {
         LWMIN  = M+N;
         LRWMIN = max( 6, N );
      }

      LQUERY = ( LWORK == -1 ) || ( LRWORK == -1 );
      if ( !( UPPER || LOWER || lsame( JOBA, 'G' ) ) ) {
         INFO = -1;
      } else if ( !( LSVEC || UCTOL || lsame( JOBU, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( RSVEC || APPLV || lsame( JOBV, 'N' ) ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( ( N < 0 ) || ( N > M ) ) {
         INFO = -5;
      } else if ( LDA < M ) {
         INFO = -7;
      } else if ( MV < 0 ) {
         INFO = -9;
      } else if ( ( RSVEC && ( LDV < N ) ) || ( APPLV && ( LDV < MV ) ) ) {
         INFO = -11;
      } else if ( UCTOL && ( RWORK( 1 ) <= ONE ) ) {
         INFO = -12;
      } else if ( LWORK < LWMIN && ( !LQUERY ) ) {
         INFO = -13;
      } else if ( LRWORK < LRWMIN && ( !LQUERY ) ) {
         INFO = -15;
      } else {
         INFO = 0;
      }

      // #:(
      if ( INFO != 0 ) {
         xerbla('ZGESVJ', -INFO );
         return;
      } else if ( LQUERY ) {
         CWORK[1] = LWMIN;
         RWORK[1] = LRWMIN;
         return;
      }

// #:) Quick return for void matrix

      if (MINMN == 0) return;

      // Set numerical parameters
      // The stopping criterion for Jacobi rotations is

      // max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS

      // where EPS is the round-off and CTOL is defined as follows:

      if ( UCTOL ) {
         // ... user controlled
         CTOL = RWORK( 1 );
      } else {
         // ... default
         if ( LSVEC || RSVEC || APPLV ) {
            CTOL = sqrt( M.toDouble() );
         } else {
            CTOL = M.toDouble();
         }
      }
      // ... and the machine dependent parameters are
// [!]  (Make sure that SLAMCH() works properly on the target machine.)

      EPSLN = dlamch( 'Epsilon' );
      ROOTEPS = sqrt( EPSLN );
      SFMIN = dlamch( 'SafeMinimum' );
      ROOTSFMIN = sqrt( SFMIN );
      SMALL = SFMIN / EPSLN;
      BIG = dlamch( 'Overflow' );
      // BIG         = ONE    / SFMIN
      ROOTBIG = ONE / ROOTSFMIN;
       // LARGE = BIG / sqrt( (M*N).toDouble() )
      BIGTHETA = ONE / ROOTEPS;

      TOL = CTOL*EPSLN;
      ROOTTOL = sqrt( TOL );

      if ( M.toDouble()*EPSLN >= ONE ) {
         INFO = -4;
         xerbla('ZGESVJ', -INFO );
         return;
      }

      // Initialize the right singular vector matrix.

      if ( RSVEC ) {
         MVL = N;
         zlaset('A', MVL, N, CZERO, CONE, V, LDV );
      } else if ( APPLV ) {
         MVL = MV;
      }
      RSVEC = RSVEC || APPLV;

      // Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N )
// (!)  If necessary, scale A to protect the largest singular value
      // from overflow. It is possible that saving the largest singular
      // value destroys the information about the small ones.
      // This initial scaling is almost minimal in the sense that the
      // goal is to make sure that no column norm overflows, and that
      // sqrt(N)*max_i SVA(i) does not overflow. If INFinite entries
      // in A are detected, the procedure returns with INFO=-6.

      SKL = ONE / sqrt( (M).toDouble()*N.toDouble() );
      NOSCALE = true;
      GOSCALE = true;

      if ( LOWER ) {
         // the input matrix is M-by-N lower triangular (trapezoidal)
         for (p = 1; p <= N; p++) { // 1874
            AAPP = ZERO;
            AAQQ = ONE;
            zlassq(M-p+1, A( p, p ), 1, AAPP, AAQQ );
            if ( AAPP > BIG ) {
               INFO = -6;
               xerbla('ZGESVJ', -INFO );
               return;
            }
            AAQQ = sqrt( AAQQ );
            if ( ( AAPP < ( BIG / AAQQ ) ) && NOSCALE ) {
               SVA[p] = AAPP*AAQQ;
            } else {
               NOSCALE = false;
               SVA[p] = AAPP*( AAQQ*SKL );
               if ( GOSCALE ) {
                  GOSCALE = false;
                  for (q = 1; q <= p - 1; q++) { // 1873
                     SVA[q] = SVA( q )*SKL;
                  } // 1873
               }
            }
         } // 1874
      } else if ( UPPER ) {
         // the input matrix is M-by-N upper triangular (trapezoidal)
         for (p = 1; p <= N; p++) { // 2874
            AAPP = ZERO;
            AAQQ = ONE;
            zlassq(p, A( 1, p ), 1, AAPP, AAQQ );
            if ( AAPP > BIG ) {
               INFO = -6;
               xerbla('ZGESVJ', -INFO );
               return;
            }
            AAQQ = sqrt( AAQQ );
            if ( ( AAPP < ( BIG / AAQQ ) ) && NOSCALE ) {
               SVA[p] = AAPP*AAQQ;
            } else {
               NOSCALE = false;
               SVA[p] = AAPP*( AAQQ*SKL );
               if ( GOSCALE ) {
                  GOSCALE = false;
                  for (q = 1; q <= p - 1; q++) { // 2873
                     SVA[q] = SVA( q )*SKL;
                  } // 2873
               }
            }
         } // 2874
      } else {
         // the input matrix is M-by-N general dense
         for (p = 1; p <= N; p++) { // 3874
            AAPP = ZERO;
            AAQQ = ONE;
            zlassq(M, A( 1, p ), 1, AAPP, AAQQ );
            if ( AAPP > BIG ) {
               INFO = -6;
               xerbla('ZGESVJ', -INFO );
               return;
            }
            AAQQ = sqrt( AAQQ );
            if ( ( AAPP < ( BIG / AAQQ ) ) && NOSCALE ) {
               SVA[p] = AAPP*AAQQ;
            } else {
               NOSCALE = false;
               SVA[p] = AAPP*( AAQQ*SKL );
               if ( GOSCALE ) {
                  GOSCALE = false;
                  for (q = 1; q <= p - 1; q++) { // 3873
                     SVA[q] = SVA( q )*SKL;
                  } // 3873
               }
            }
         } // 3874
      }

      if (NOSCALE) SKL = ONE;

      // Move the smaller part of the spectrum from the underflow threshold
// (!)  Start by determining the position of the nonzero entries of the
      // array SVA() relative to ( SFMIN, BIG ).

      AAPP = ZERO;
      AAQQ = BIG;
      for (p = 1; p <= N; p++) { // 4781
         if( SVA( p ) != ZERO )AAQQ = min( AAQQ, SVA( p ) );
         AAPP = max( AAPP, SVA( p ) );
      } // 4781

// #:) Quick return for zero matrix

      if ( AAPP == ZERO ) {
         if (LSVEC) zlaset( 'G', M, N, CZERO, CONE, A, LDA );
         RWORK[1] = ONE;
         RWORK[2] = ZERO;
         RWORK[3] = ZERO;
         RWORK[4] = ZERO;
         RWORK[5] = ZERO;
         RWORK[6] = ZERO;
         return;
      }

// #:) Quick return for one-column matrix

      if ( N == 1 ) {
         if (LSVEC) zlascl( 'G', 0, 0, SVA( 1 ), SKL, M, 1, A( 1, 1 ), LDA, IERR );
         RWORK[1] = ONE / SKL;
         if ( SVA( 1 ) >= SFMIN ) {
            RWORK[2] = ONE;
         } else {
            RWORK[2] = ZERO;
         }
         RWORK[3] = ZERO;
         RWORK[4] = ZERO;
         RWORK[5] = ZERO;
         RWORK[6] = ZERO;
         return;
      }

      // Protect small singular values from underflow, and try to
      // avoid underflows/overflows in computing Jacobi rotations.

      SN = sqrt( SFMIN / EPSLN );
      TEMP1 = sqrt( BIG / N.toDouble() );
      if ( ( AAPP <= SN ) || ( AAQQ >= TEMP1 ) || ( ( SN <= AAQQ ) && ( AAPP <= TEMP1 ) ) ) {
         TEMP1 = min( BIG, TEMP1 / AAPP );
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ <= SN ) && ( AAPP <= TEMP1 ) ) {
         TEMP1 = min( SN / AAQQ, BIG / (AAPP*sqrt( N.toDouble()) ) );
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ >= SN ) && ( AAPP >= TEMP1 ) ) {
         TEMP1 = max( SN / AAQQ, TEMP1 / AAPP );
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ <= SN ) && ( AAPP >= TEMP1 ) ) {
         TEMP1 = min( SN / AAQQ, BIG / ( sqrt( N.toDouble() )*AAPP ) );
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else {
         TEMP1 = ONE;
      }

      // Scale, if necessary

      if ( TEMP1 != ONE ) {
         dlascl('G', 0, 0, ONE, TEMP1, N, 1, SVA, N, IERR );
      }
      SKL = TEMP1*SKL;
      if ( SKL != ONE ) {
         zlascl(JOBA, 0, 0, ONE, SKL, M, N, A, LDA, IERR );
         SKL = ONE / SKL;
      }

      // Row-cyclic Jacobi SVD algorithm with column pivoting

      EMPTSW = ( N*( N-1 ) ) / 2;
      NOTROT = 0;

      for (q = 1; q <= N; q++) { // 1868
         CWORK[q] = CONE;
      } // 1868



      SWBAND = 3;
// [TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
      // if ZGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm ZGEJSV. For sweeps i=1:SWBAND the procedure
      // works on pivots inside a band-like region around the diagonal.
      // The boundaries are determined dynamically, based on the number of
      // pivots above a threshold.

      KBL = min( 8, N );
// [TP] KBL is a tuning parameter that defines the tile size in the
      // tiling of the p-q loops of pivot pairs. In general, an optimal
      // value of KBL depends on the matrix dimensions and on the
      // parameters of the computer's memory.

      NBL = N / KBL;
      if( ( NBL*KBL ) != N )NBL = NBL + 1;

      BLSKIP = KBL**2;
// [TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = min( 5, KBL );
// [TP] ROWSKIP is a tuning parameter.

      LKAHEAD = 1;
// [TP] LKAHEAD is a tuning parameter.

      // Quasi block transformations, using the lower (upper) triangular
      // structure of the input matrix. The quasi-block-cycling usually
      // invokes cubic convergence. Big part of this cycle is done inside
      // canonical subspaces of dimensions less than M.

      if ( ( LOWER || UPPER ) && ( N > max( 64, 4*KBL ) ) ) {
// [TP] The number of partition levels and the actual partition are
      // tuning parameters.
         N4 = N / 4;
         N2 = N / 2;
         N34 = 3*N4;
         if ( APPLV ) {
            q = 0;
         } else {
            q = 1;
         }

         if ( LOWER ) {

      // This works very well on lower triangular matrices, in particular
      // in the framework of the preconditioned Jacobi SVD (xGEJSV).
      // The idea is simple:
      // [+ 0 0 0]   Note that Jacobi transformations of [0 0]
      // [+ + 0 0]                                       [0 0]
      // [+ + x 0]   actually work on [x 0]              [x 0]
      // [+ + x x]                    [x x].             [x x]

            zgsvj0(JOBV, M-N34, N-N34, A( N34+1, N34+1 ), LDA, CWORK( N34+1 ), SVA( N34+1 ), MVL, V( N34*q+1, N34+1 ), LDV, EPSLN, SFMIN, TOL, 2, CWORK( N+1 ), LWORK-N, IERR );
             zgsvj0(JOBV, M-N2, N34-N2, A( N2+1, N2+1 ), LDA, CWORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 2, CWORK( N+1 ), LWORK-N, IERR );
             zgsvj1(JOBV, M-N2, N-N2, N4, A( N2+1, N2+1 ), LDA, CWORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );
             zgsvj0(JOBV, M-N4, N2-N4, A( N4+1, N4+1 ), LDA, CWORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            zgsvj0(JOBV, M, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            zgsvj1(JOBV, M, N2, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );


         } else if ( UPPER ) {


            zgsvj0(JOBV, N4, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 2, CWORK( N+1 ), LWORK-N, IERR );

            zgsvj0(JOBV, N2, N4, A( 1, N4+1 ), LDA, CWORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            zgsvj1(JOBV, N2, N2, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            zgsvj0(JOBV, N2+N4, N4, A( 1, N2+1 ), LDA, CWORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

         }

      }

      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

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

         for (ibr = 1; ibr <= NBL; ibr++) { // 2000

            igl = ( ibr-1 )*KBL + 1;

            for (ir1 = 0; ir1 <= min( LKAHEAD, NBL-ibr ); ir1++) { // 1002

               igl = igl + ir1*KBL;

               for (p = igl; p <= min( igl+KBL-1, N-1 ); p++) { // 2001

      // .. de Rijk's pivoting

                  q = idamax( N-p+1, SVA( p ), 1 ) + p - 1;
                  if ( p != q ) {
                     zswap(M, A( 1, p ), 1, A( 1, q ), 1 );
                     if (RSVEC) zswap( MVL, V( 1, p ), 1, V( 1, q ), 1 );
                     TEMP1 = SVA( p );
                     SVA[p] = SVA( q );
                     SVA[q] = TEMP1;
                     AAPQ = CWORK(p);
                     CWORK[p] = CWORK(q);
                     CWORK[q] = AAPQ;
                  }

                  if ( ir1 == 0 ) {

         // Column norms are periodically updated by explicit
         // norm computation.
// [!]     Caveat:
         // Unfortunately, some BLAS implementations compute DZNRM2(M,A(1,p),1)
         // as sqrt(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to
         // overflow for ||A(:,p)||_2 > sqrt(overflow_threshold), and to
         // underflow for ||A(:,p)||_2 < sqrt(underflow_threshold).
         // Hence, DZNRM2 cannot be trusted, not even in the case when
         // the true norm is far from the under(over)flow boundaries.
         // If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF
         // below should be replaced with "AAPP = DZNRM2( M, A(1,p), 1 )".

                     if ( ( SVA( p ) < ROOTBIG ) && ( SVA( p ) > ROOTSFMIN ) ) {
                        SVA[p] = DZNRM2( M, A( 1, p ), 1 );
                     } else {
                        TEMP1 = ZERO;
                        AAPP = ONE;
                        zlassq(M, A( 1, p ), 1, TEMP1, AAPP );
                        SVA[p] = TEMP1*sqrt( AAPP );
                     }
                     AAPP = SVA( p );
                  } else {
                     AAPP = SVA( p );
                  }

                  if ( AAPP > ZERO ) {

                     PSKIPPED = 0;

                     for (q = p + 1; q <= min( igl+KBL-1, N ); q++) { // 2002

                        AAQQ = SVA( q );

                        if ( AAQQ > ZERO ) {

                           AAPP0 = AAPP;
                           if ( AAQQ >= ONE ) {
                              ROTOK = ( SMALL*AAPP ) <= AAQQ;
                              if ( AAPP < ( BIG / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP;
                              } else {
                                 zcopy(M, A( 1, p ), 1, CWORK(N+1), 1 );
                                 CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1), LDA, IERR )                                  AAPQ = ZDOTC( M, CWORK(N+1), 1, A( 1, q ), 1 ) / AAQQ;
                              }
                           } else {
                              ROTOK = AAPP <= ( AAQQ / SMALL );
                              if ( AAPP > ( SMALL / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAPP ) / AAQQ;
                              } else {
                                 zcopy(M, A( 1, q ), 1, CWORK(N+1), 1 );
                                 zlascl('G', 0, 0, AAQQ, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 AAPQ = ZDOTC( M, A(1, p ), 1, CWORK(N+1), 1 ) / AAPP;
                              }
                           }


                            // AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q)
                           AAPQ1  = -(AAPQ).abs();
                           MXAAPQ = max( MXAAPQ, -AAPQ1 );

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ( AAPQ1 ).abs() > TOL ) {
                           OMPQ = AAPQ / (AAPQ).abs();

            // .. rotate
// [RTD]      ROTATED = ROTATED + ONE

                              if ( ir1 == 0 ) {
                                 NOTROT = 0;
                                 PSKIPPED = 0;
                                 ISWROT = ISWROT + 1;
                              }

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP;
                                 APOAQ = AAPP / AAQQ;
                                 THETA = -HALF*( AQOAP-APOAQ ).abs()/AAPQ1;

                                 if ( ( THETA ).abs() > BIGTHETA ) {

                                    T  = HALF / THETA;
                                    CS = ONE;
                                     zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
                                    }
                                     SVA[q] = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ1 ) );
                                    MXSINJ = max( MXSINJ, ( T ).abs() );

                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -sign( ONE, AAPQ1 );
                                    T = ONE / ( THETA+THSIGN* sqrt( ONE+THETA*THETA ) );
                                    CS = sqrt( ONE / ( ONE+T*T ) );
                                    SN = T*CS;

                                    MXSINJ = max( MXSINJ, ( SN ).abs() );
                                    SVA[q] = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ1 ) );

                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 CWORK[p] = -CWORK(q) * OMPQ;

                                 } else {
               // .. have to use modified Gram-Schmidt like transformation
                                 zcopy(M, A( 1, p ), 1, CWORK(N+1), 1 );
                                 zlascl('G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 zlascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                 zaxpy(M, -AAPQ, CWORK(N+1), 1, A( 1, q ), 1 );
                                 zlascl['G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )                                  SVA( q] = AAQQ*sqrt( max( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                 MXSINJ = max( MXSINJ, SFMIN );
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // recompute SVA(q), SVA(p).

                              if ( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) {
                                 if( ( AAQQ < ROOTBIG ) && ( AAQQ > ROOTSFMIN ) ) {
                                    SVA[q] = DZNRM2( M, A( 1, q ), 1 );
                                 } else {
                                    T = ZERO;
                                    AAQQ = ONE;
                                    zlassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA[q] = T*sqrt( AAQQ );
                                 }
                              }
                              if ( ( AAPP / AAPP0 ) <= ROOTEPS ) {
                                 if ( ( AAPP < ROOTBIG ) && ( AAPP > ROOTSFMIN ) ) {
                                    AAPP = DZNRM2( M, A( 1, p ), 1 );
                                 } else {
                                    T = ZERO;
                                    AAPP = ONE;
                                    zlassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*sqrt( AAPP );
                                 }
                                 SVA[p] = AAPP;
                              }

                           } else {
                              // A(:,p) and A(:,q) already numerically orthogonal
                              if (ir1 == 0) NOTROT = NOTROT + 1;
// [RTD]      SKIPPED  = SKIPPED + 1
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

                     SVA[p] = AAPP;

                  } else {
                     SVA[p] = AAPP;
                     if( ( ir1 == 0 ) && ( AAPP == ZERO ) ) NOTROT = NOTROT + min( igl+KBL-1, N ) - p;
                  }

               } // 2001
      // end of the p-loop
      // end of doing the block ( ibr, ibr )
            } // 1002
      // end of ir1-loop

// ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1;

            for (jbc = ibr + 1; jbc <= NBL; jbc++) { // 2010

               jgl = ( jbc-1 )*KBL + 1;

         // doing the block at ( ibr, jbc )

               IJBLSK = 0;
               for (p = igl; p <= min( igl+KBL-1, N ); p++) { // 2100

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
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP;
                              } else {
                                 zcopy(M, A( 1, p ), 1, CWORK(N+1), 1 );
                                 zlascl('G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 AAPQ = ZDOTC( M, CWORK(N+1), 1, A( 1, q ), 1 ) / AAQQ;
                              }
                           } else {
                              if ( AAPP >= AAQQ ) {
                                 ROTOK = AAPP <= ( AAQQ / SMALL );
                              } else {
                                 ROTOK = AAQQ <= ( AAPP / SMALL );
                              }
                              if ( AAPP > ( SMALL / AAQQ ) ) {
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / max(AAQQ,AAPP) ) / min(AAQQ,AAPP);
                              } else {
                                 zcopy(M, A( 1, q ), 1, CWORK(N+1), 1 );
                                 zlascl('G', 0, 0, AAQQ, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 AAPQ = ZDOTC( M, A( 1, p ), 1, CWORK(N+1),  1 ) / AAPP;
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
                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
                                    }
                                    SVA[q] = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ1 ) );
                                    MXSINJ = max( MXSINJ, ( T ).abs() );
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -sign( ONE, AAPQ1 );
                                    if (AAQQ > AAPP0) THSIGN = -THSIGN;
                                    T = ONE / ( THETA+THSIGN* sqrt( ONE+THETA*THETA ) );
                                    CS = sqrt( ONE / ( ONE+T*T ) );
                                    SN = T*CS;
                                    MXSINJ = max( MXSINJ, ( SN ).abs() );
                                    SVA[q] = AAQQ*sqrt( max( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*sqrt( max( ZERO, ONE-T*AQOAP*AAPQ1 ) );

                                    zrot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        zrot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 CWORK[p] = -CWORK(q) * OMPQ;

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                               if ( AAPP > AAQQ ) {
                                    zcopy(M, A( 1, p ), 1, CWORK(N+1), 1 );
                                    zlascl('G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1),LDA, IERR );
                                    zlascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                    zaxpy(M, -AAPQ, CWORK(N+1), 1, A( 1, q ), 1 );
                                    zlascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA[q] = AAQQ*sqrt( max( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                               } else {
                                   zcopy(M, A( 1, q ), 1, CWORK(N+1), 1 );
                                    zlascl('G', 0, 0, AAQQ, ONE, M, 1, CWORK(N+1),LDA, IERR );
                                    zlascl('G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR );
                                    zaxpy(M, -CONJG(AAPQ), CWORK(N+1), 1, A( 1, p ), 1 );
                                    zlascl('G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA[p] = AAPP*sqrt( max( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                    MXSINJ = max( MXSINJ, SFMIN );
                               }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // .. recompute SVA(q), SVA(p)
                              if ( ( SVA( q ) / AAQQ )**2 <= ROOTEPS ) {
                                 if( ( AAQQ < ROOTBIG ) && ( AAQQ > ROOTSFMIN ) ) {
                                    SVA[q] = DZNRM2( M, A( 1, q ), 1);
                                  } else {
                                    T = ZERO;
                                    AAQQ = ONE;
                                    zlassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA[q] = T*sqrt( AAQQ );
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2 <= ROOTEPS ) {
                                 if ( ( AAPP < ROOTBIG ) && ( AAPP > ROOTSFMIN ) ) {
                                    AAPP = DZNRM2( M, A( 1, p ), 1 );
                                 } else {
                                    T = ZERO;
                                    AAPP = ONE;
                                    zlassq(M, A( 1, p ), 1, T, AAPP );
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
            SVA[N] = DZNRM2( M, A( 1, N ), 1 );
         } else {
            T = ZERO;
            AAPP = ONE;
            zlassq(M, A( 1, N ), 1, T, AAPP );
            SVA[N] = T*sqrt( AAPP );
         }

      // Additional steering devices

         if( ( i < SWBAND ) && ( ( MXAAPQ <= ROOTTOL ) || ( ISWROT <= N ) ) )SWBAND = i;

         if ( ( i > SWBAND+1 ) && ( MXAAPQ < sqrt( (N).toDouble() )* TOL ) && ( N.toDouble()*MXAAPQ*MXSINJ < TOL ) ) {
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

      // Sort the singular values and find how many are above
      // the underflow threshold.

      N2 = 0;
      N4 = 0;
      for (p = 1; p <= N - 1; p++) { // 5991
         q = idamax( N-p+1, SVA( p ), 1 ) + p - 1;
         if ( p != q ) {
            TEMP1 = SVA( p );
            SVA[p] = SVA( q );
            SVA[q] = TEMP1;
            zswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            if (RSVEC) zswap( MVL, V( 1, p ), 1, V( 1, q ), 1 );
         }
         if ( SVA( p ) != ZERO ) {
            N4 = N4 + 1;
            if( SVA( p )*SKL > SFMIN )N2 = N2 + 1;
         }
      } // 5991
      if ( SVA( N ) != ZERO ) {
         N4 = N4 + 1;
         if( SVA( N )*SKL > SFMIN )N2 = N2 + 1;
      }

      // Normalize the left singular vectors.

      if ( LSVEC || UCTOL ) {
         for (p = 1; p <= N4; p++) { // 1998
             // CALL ZDSCAL( M, ONE / SVA( p ), A( 1, p ), 1 )
            zlascl('G',0,0, SVA(p), ONE, M, 1, A(1,p), M, IERR );
         } // 1998
      }

      // Scale the product of Jacobi rotations.

      if ( RSVEC ) {
            for (p = 1; p <= N; p++) { // 2399
               TEMP1 = ONE / DZNRM2( MVL, V( 1, p ), 1 );
               zdscal(MVL, TEMP1, V( 1, p ), 1 );
            } // 2399
      }

      // Undo scaling, if necessary (and possible).
      if ( ( ( SKL > ONE ) && ( SVA( 1 ) < ( BIG / SKL ) ) ) || ( ( SKL < ONE ) && ( SVA( max( N2, 1 ) ) > ( SFMIN / SKL ) ) ) ) {
         for (p = 1; p <= N; p++) { // 2400
            SVA[p] = SKL*SVA( p );
         } // 2400
         SKL = ONE;
      }

      RWORK[1] = SKL;
      // The singular values of A are SKL*SVA(1:N). If SKL != ONE
      // then some of the singular values may overflow or underflow and
      // the spectrum is given in this factored representation.

      RWORK[2] = N4.toDouble();
      // N4 is the number of computed nonzero singular values of A.

      RWORK[3] = N2.toDouble();
      // N2 is the number of singular values of A greater than SFMIN.
      // If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
      // that may carry some information.

      RWORK[4] = i.toDouble();
      // i is the index of the last sweep before declaring convergence.

      RWORK[5] = MXAAPQ;
      // MXAAPQ is the largest absolute value of scaled pivots in the
      // last sweep

      RWORK[6] = MXSINJ;
      // MXSINJ is the largest absolute value of the sines of Jacobi angles
      // in the last sweep

      return;
      // ..
      // .. END OF ZGESVJ
      // ..
      }
