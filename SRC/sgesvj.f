      SUBROUTINE SGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, LDV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDV, LWORK, M, MV, N;
      String             JOBA, JOBU, JOBV;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), SVA( N ), V( LDV, * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Local Parameters ..
      REAL               ZERO, HALF, ONE
      const              ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0;
      int                NSWEEP;
      const              NSWEEP = 30 ;
      // ..
      // .. Local Scalars ..
      REAL               AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, CTOL, EPSLN, LARGE, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SKL, SFMIN, SMALL, SN, T, TEMP1, THETA, THSIGN, TOL
      int                BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, N2, N34, N4, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND, MINMN, LWMIN;
      bool               APPLV, GOSCALE, LOWER, LQUERY, LSVEC, NOSCALE, ROTOK, RSVEC, UCTOL, UPPER;
      // ..
      // .. Local Arrays ..
      REAL               FASTR( 5 )
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, FLOAT, SIGN, SQRT
      // ..
      // .. External Functions ..
      // ..
      // from BLAS
      REAL               SDOT, SNRM2
      // EXTERNAL SDOT, SNRM2
      int                ISAMAX;
      // EXTERNAL ISAMAX
      // from LAPACK
      REAL               SLAMCH, SROUNDUP_LWORK
      // EXTERNAL SLAMCH, SROUNDUP_LWORK
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // ..
      // from BLAS
      // EXTERNAL SAXPY, SCOPY, SROTM, SSCAL, SSWAP
      // from LAPACK
      // EXTERNAL SLASCL, SLASET, SLASSQ, XERBLA

      // EXTERNAL SGSVJ0, SGSVJ1
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      LSVEC = LSAME( JOBU, 'U' )
      UCTOL = LSAME( JOBU, 'C' )
      RSVEC = LSAME( JOBV, 'V' )
      APPLV = LSAME( JOBV, 'A' )
      UPPER = LSAME( JOBA, 'U' )
      LOWER = LSAME( JOBA, 'L' )

      MINMN = MIN( M, N )
      if ( MINMN == 0 ) {
         LWMIN = 1
      } else {
         LWMIN = MAX( 6, M+N )
      }

      LQUERY = ( LWORK == -1 )
      if ( .NOT.( UPPER .OR. LOWER .OR. LSAME( JOBA, 'G' ) ) ) {
         INFO = -1
      } else if ( .NOT.( LSVEC .OR. UCTOL .OR. LSAME( JOBU, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( RSVEC .OR. APPLV .OR. LSAME( JOBV, 'N' ) ) ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( ( N.LT.0 ) .OR. ( N.GT.M ) ) {
         INFO = -5
      } else if ( LDA.LT.M ) {
         INFO = -7
      } else if ( MV.LT.0 ) {
         INFO = -9
      } else if ( ( RSVEC && ( LDV.LT.N ) ) .OR. ( APPLV && ( LDV.LT.MV ) ) ) {
         INFO = -11
      } else if ( UCTOL && ( WORK( 1 ).LE.ONE ) ) {
         INFO = -12
      } else if ( LWORK.LT.LWMIN && ( .NOT.LQUERY ) ) {
         INFO = -13
      } else {
         INFO = 0
      }

      // #:(
      if ( INFO != 0 ) {
         xerbla('SGESVJ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
         RETURN
      }

* #:) Quick return for void matrix

      if (MINMN == 0) RETURN;

      // Set numerical parameters
      // The stopping criterion for Jacobi rotations is

      // max_{i<>j}|A(:,i)^T * A(:,j)|/(||A(:,i)||*||A(:,j)||) < CTOL*EPS

      // where EPS is the round-off and CTOL is defined as follows:

      if ( UCTOL ) {
         // ... user controlled
         CTOL = WORK( 1 )
      } else {
         // ... default
         if ( LSVEC .OR. RSVEC .OR. APPLV ) {
            CTOL = SQRT( FLOAT( M ) )
         } else {
            CTOL = FLOAT( M )
         }
      }
      // ... and the machine dependent parameters are
*[!]  (Make sure that SLAMCH() works properly on the target machine.)

      EPSLN = SLAMCH( 'Epsilon' )
      ROOTEPS = SQRT( EPSLN )
      SFMIN = SLAMCH( 'SafeMinimum' )
      ROOTSFMIN = SQRT( SFMIN )
      SMALL = SFMIN / EPSLN
      BIG = SLAMCH( 'Overflow' )
      // BIG         = ONE    / SFMIN
      ROOTBIG = ONE / ROOTSFMIN
      LARGE = BIG / SQRT( FLOAT( M*N ) )
      BIGTHETA = ONE / ROOTEPS

      TOL = CTOL*EPSLN
      ROOTTOL = SQRT( TOL )

      if ( FLOAT( M )*EPSLN.GE.ONE ) {
         INFO = -4
         xerbla('SGESVJ', -INFO );
         RETURN
      }

      // Initialize the right singular vector matrix.

      if ( RSVEC ) {
         MVL = N
         slaset('A', MVL, N, ZERO, ONE, V, LDV );
      } else if ( APPLV ) {
         MVL = MV
      }
      RSVEC = RSVEC .OR. APPLV

      // Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N )
*(!)  If necessary, scale A to protect the largest singular value
      // from overflow. It is possible that saving the largest singular
      // value destroys the information about the small ones.
      // This initial scaling is almost minimal in the sense that the
      // goal is to make sure that no column norm overflows, and that
      // SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries
      // in A are detected, the procedure returns with INFO=-6.

      SKL = ONE / SQRT( FLOAT( M )*FLOAT( N ) )
      NOSCALE = true;
      GOSCALE = true;

      if ( LOWER ) {
         // the input matrix is M-by-N lower triangular (trapezoidal)
         for (p = 1; p <= N; p++) { // 1874
            AAPP = ZERO
            AAQQ = ONE
            slassq(M-p+1, A( p, p ), 1, AAPP, AAQQ );
            if ( AAPP.GT.BIG ) {
               INFO = -6
               xerbla('SGESVJ', -INFO );
               RETURN
            }
            AAQQ = SQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) && NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = false;
               SVA( p ) = AAPP*( AAQQ*SKL )
               if ( GOSCALE ) {
                  GOSCALE = false;
                  for (q = 1; q <= p - 1; q++) { // 1873
                     SVA( q ) = SVA( q )*SKL
                  } // 1873
               }
            }
         } // 1874
      } else if ( UPPER ) {
         // the input matrix is M-by-N upper triangular (trapezoidal)
         for (p = 1; p <= N; p++) { // 2874
            AAPP = ZERO
            AAQQ = ONE
            slassq(p, A( 1, p ), 1, AAPP, AAQQ );
            if ( AAPP.GT.BIG ) {
               INFO = -6
               xerbla('SGESVJ', -INFO );
               RETURN
            }
            AAQQ = SQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) && NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = false;
               SVA( p ) = AAPP*( AAQQ*SKL )
               if ( GOSCALE ) {
                  GOSCALE = false;
                  for (q = 1; q <= p - 1; q++) { // 2873
                     SVA( q ) = SVA( q )*SKL
                  } // 2873
               }
            }
         } // 2874
      } else {
         // the input matrix is M-by-N general dense
         for (p = 1; p <= N; p++) { // 3874
            AAPP = ZERO
            AAQQ = ONE
            slassq(M, A( 1, p ), 1, AAPP, AAQQ );
            if ( AAPP.GT.BIG ) {
               INFO = -6
               xerbla('SGESVJ', -INFO );
               RETURN
            }
            AAQQ = SQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) && NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = false;
               SVA( p ) = AAPP*( AAQQ*SKL )
               if ( GOSCALE ) {
                  GOSCALE = false;
                  for (q = 1; q <= p - 1; q++) { // 3873
                     SVA( q ) = SVA( q )*SKL
                  } // 3873
               }
            }
         } // 3874
      }

      if (NOSCALE) SKL = ONE;

      // Move the smaller part of the spectrum from the underflow threshold
*(!)  Start by determining the position of the nonzero entries of the
      // array SVA() relative to ( SFMIN, BIG ).

      AAPP = ZERO
      AAQQ = BIG
      for (p = 1; p <= N; p++) { // 4781
         IF( SVA( p ) != ZERO )AAQQ = MIN( AAQQ, SVA( p ) )
         AAPP = MAX( AAPP, SVA( p ) )
      } // 4781

* #:) Quick return for zero matrix

      if ( AAPP == ZERO ) {
         if (LSVEC) CALL SLASET( 'G', M, N, ZERO, ONE, A, LDA );
         WORK( 1 ) = ONE
         WORK( 2 ) = ZERO
         WORK( 3 ) = ZERO
         WORK( 4 ) = ZERO
         WORK( 5 ) = ZERO
         WORK( 6 ) = ZERO
         RETURN
      }

* #:) Quick return for one-column matrix

      if ( N == 1 ) {
         if (LSVEC) CALL SLASCL( 'G', 0, 0, SVA( 1 ), SKL, M, 1, A( 1, 1 ), LDA, IERR );
         WORK( 1 ) = ONE / SKL
         if ( SVA( 1 ).GE.SFMIN ) {
            WORK( 2 ) = ONE
         } else {
            WORK( 2 ) = ZERO
         }
         WORK( 3 ) = ZERO
         WORK( 4 ) = ZERO
         WORK( 5 ) = ZERO
         WORK( 6 ) = ZERO
         RETURN
      }

      // Protect small singular values from underflow, and try to
      // avoid underflows/overflows in computing Jacobi rotations.

      SN = SQRT( SFMIN / EPSLN )
      TEMP1 = SQRT( BIG / FLOAT( N ) )
      if ( ( AAPP.LE.SN ) .OR. ( AAQQ.GE.TEMP1 ) .OR. ( ( SN.LE.AAQQ ) && ( AAPP.LE.TEMP1 ) ) ) {
         TEMP1 = MIN( BIG, TEMP1 / AAPP )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.LE.SN ) && ( AAPP.LE.TEMP1 ) ) {
         TEMP1 = MIN( SN / AAQQ, BIG / ( AAPP*SQRT( FLOAT( N ) ) ) )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.GE.SN ) && ( AAPP.GE.TEMP1 ) ) {
         TEMP1 = MAX( SN / AAQQ, TEMP1 / AAPP )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.LE.SN ) && ( AAPP.GE.TEMP1 ) ) {
         TEMP1 = MIN( SN / AAQQ, BIG / ( SQRT( FLOAT( N ) )*AAPP ) )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else {
         TEMP1 = ONE
      }

      // Scale, if necessary

      if ( TEMP1 != ONE ) {
         slascl('G', 0, 0, ONE, TEMP1, N, 1, SVA, N, IERR );
      }
      SKL = TEMP1*SKL
      if ( SKL != ONE ) {
         slascl(JOBA, 0, 0, ONE, SKL, M, N, A, LDA, IERR );
         SKL = ONE / SKL
      }

      // Row-cyclic Jacobi SVD algorithm with column pivoting

      EMPTSW = ( N*( N-1 ) ) / 2
      NOTROT = 0
      FASTR( 1 ) = ZERO

      // A is represented in factored form A = A * diag(WORK), where diag(WORK)
      // is initialized to identity. WORK is updated during fast scaled
      // rotations.

      for (q = 1; q <= N; q++) { // 1868
         WORK( q ) = ONE
      } // 1868


      SWBAND = 3
*[TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
      // if SGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
      // works on pivots inside a band-like region around the diagonal.
      // The boundaries are determined dynamically, based on the number of
      // pivots above a threshold.

      KBL = MIN( 8, N )
*[TP] KBL is a tuning parameter that defines the tile size in the
      // tiling of the p-q loops of pivot pairs. In general, an optimal
      // value of KBL depends on the matrix dimensions and on the
      // parameters of the computer's memory.

      NBL = N / KBL
      IF( ( NBL*KBL ) != N )NBL = NBL + 1

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

      if ( ( LOWER .OR. UPPER ) && ( N.GT.MAX( 64, 4*KBL ) ) ) {
*[TP] The number of partition levels and the actual partition are
      // tuning parameters.
         N4 = N / 4
         N2 = N / 2
         N34 = 3*N4
         if ( APPLV ) {
            q = 0
         } else {
            q = 1
         }

         if ( LOWER ) {

      // This works very well on lower triangular matrices, in particular
      // in the framework of the preconditioned Jacobi SVD (xGEJSV).
      // The idea is simple:
      // [+ 0 0 0]   Note that Jacobi transformations of [0 0]
      // [+ + 0 0]                                       [0 0]
      // [+ + x 0]   actually work on [x 0]              [x 0]
      // [+ + x x]                    [x x].             [x x]

            sgsvj0(JOBV, M-N34, N-N34, A( N34+1, N34+1 ), LDA, WORK( N34+1 ), SVA( N34+1 ), MVL, V( N34*q+1, N34+1 ), LDV, EPSLN, SFMIN, TOL, 2, WORK( N+1 ), LWORK-N, IERR );

            sgsvj0(JOBV, M-N2, N34-N2, A( N2+1, N2+1 ), LDA, WORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 2, WORK( N+1 ), LWORK-N, IERR );

            sgsvj1(JOBV, M-N2, N-N2, N4, A( N2+1, N2+1 ), LDA, WORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR );

            sgsvj0(JOBV, M-N4, N2-N4, A( N4+1, N4+1 ), LDA, WORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR );

            sgsvj0(JOBV, M, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR );

            sgsvj1(JOBV, M, N2, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR );


         } else if ( UPPER ) {


            sgsvj0(JOBV, N4, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 2, WORK( N+1 ), LWORK-N, IERR );

            sgsvj0(JOBV, N2, N4, A( 1, N4+1 ), LDA, WORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR );

            sgsvj1(JOBV, N2, N2, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR );

            sgsvj0(JOBV, N2+N4, N4, A( 1, N2+1 ), LDA, WORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR );

         }

      }

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

                  q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
                  if ( p != q ) {
                     sswap(M, A( 1, p ), 1, A( 1, q ), 1 );
                     if (RSVEC) CALL SSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 );
                     TEMP1 = SVA( p )
                     SVA( p ) = SVA( q )
                     SVA( q ) = TEMP1
                     TEMP1 = WORK( p )
                     WORK( p ) = WORK( q )
                     WORK( q ) = TEMP1
                  }

                  if ( ir1 == 0 ) {

         // Column norms are periodically updated by explicit
         // norm computation.
         // Caveat:
         // Unfortunately, some BLAS implementations compute SNRM2(M,A(1,p),1)
         // as SQRT(SDOT(M,A(1,p),1,A(1,p),1)), which may cause the result to
         // overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to
         // underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
         // Hence, SNRM2 cannot be trusted, not even in the case when
         // the true norm is far from the under(over)flow boundaries.
         // If properly implemented SNRM2 is available, the IF-THEN-ELSE
         // below should read "AAPP = SNRM2( M, A(1,p), 1 ) * WORK(p)".

                     if ( ( SVA( p ).LT.ROOTBIG ) && ( SVA( p ).GT.ROOTSFMIN ) ) {
                        SVA( p ) = SNRM2( M, A( 1, p ), 1 )*WORK( p )
                     } else {
                        TEMP1 = ZERO
                        AAPP = ONE
                        slassq(M, A( 1, p ), 1, TEMP1, AAPP );
                        SVA( p ) = TEMP1*SQRT( AAPP )*WORK( p )
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
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 scopy(M, A( 1, p ), 1, WORK( N+1 ), 1 );
                                 slascl('G', 0, 0, AAPP, WORK( p ), M, 1, WORK( N+1 ), LDA, IERR );
                                 AAPQ = SDOT( M, WORK( N+1 ), 1, A( 1, q ), 1 )*WORK( q ) / AAQQ
                              }
                           } else {
                              ROTOK = AAPP.LE.( AAQQ / SMALL )
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 scopy(M, A( 1, q ), 1, WORK( N+1 ), 1 );
                                 slascl('G', 0, 0, AAQQ, WORK( q ), M, 1, WORK( N+1 ), LDA, IERR );
                                 AAPQ = SDOT( M, WORK( N+1 ), 1, A( 1, p ), 1 )*WORK( p ) / AAPP
                              }
                           }

                           MXAAPQ = MAX( MXAAPQ, ABS( AAPQ ) )

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ABS( AAPQ ).GT.TOL ) {

            // .. rotate
*[RTD]      ROTATED = ROTATED + ONE

                              if ( ir1 == 0 ) {
                                 NOTROT = 0
                                 PSKIPPED = 0
                                 ISWROT = ISWROT + 1
                              }

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*ABS( AQOAP-APOAQ ) / AAPQ

                                 if ( ABS( THETA ).GT.BIGTHETA ) {

                                    T = HALF / THETA
                                    FASTR( 3 ) = T*WORK( p ) / WORK( q )
                                    FASTR( 4 ) = -T*WORK( q ) / WORK( p );
                                    srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, ABS( T ) )

                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ )
                                    T = ONE / ( THETA+THSIGN* SQRT( ONE+THETA*THETA ) )
                                    CS = SQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS

                                    MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )

                                    APOAQ = WORK( p ) / WORK( q )
                                    AQOAP = WORK( q ) / WORK( p )
                                    if ( WORK( p ).GE.ONE ) {
                                       if ( WORK( q ).GE.ONE ) {
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q )*CS
                                          srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                       } else {
                                          saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q ) / CS
                                          if ( RSVEC ) {
                                             saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                          }
                                       }
                                    } else {
                                       if ( WORK( q ).GE.ONE ) {
                                          saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          WORK( p ) = WORK( p ) / CS
                                          WORK( q ) = WORK( q )*CS
                                          if ( RSVEC ) {
                                             saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             saxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                          }
                                       } else {
                                          IF( WORK( p ).GE.WORK( q ) ) THEN;
                                             saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             WORK( p ) = WORK( p )*CS
                                             WORK( q ) = WORK( q ) / CS
                                             if ( RSVEC ) {
                                                saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                                saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             }
                                          } else {
                                             saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             WORK( p ) = WORK( p ) / CS
                                             WORK( q ) = WORK( q )*CS
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
                                 scopy(M, A( 1, p ), 1, WORK( N+1 ), 1 );
                                 slascl('G', 0, 0, AAPP, ONE, M, 1, WORK( N+1 ), LDA, IERR );
                                 slascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                 TEMP1 = -AAPQ*WORK( p ) / WORK( q )
                                 saxpy(M, TEMP1, WORK( N+1 ), 1, A( 1, q ), 1 );
                                 slascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )                                  SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) );
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // recompute SVA(q), SVA(p).

                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) && ( AAQQ.GT.ROOTSFMIN ) ) THEN                                     SVA( q ) = SNRM2( M, A( 1, q ), 1 )* WORK( q )
                                 } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    slassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ )*WORK( q )
                                 }
                              }
                              if ( ( AAPP / AAPP0 ).LE.ROOTEPS ) {
                                 IF( ( AAPP.LT.ROOTBIG ) && ( AAPP.GT.ROOTSFMIN ) ) THEN                                     AAPP = SNRM2( M, A( 1, p ), 1 )* WORK( p )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    slassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*SQRT( AAPP )*WORK( p )
                                 }
                                 SVA( p ) = AAPP
                              }

                           } else {
         // A(:,p) and A(:,q) already numerically orthogonal
                              if (ir1 == 0) NOTROT = NOTROT + 1;
*[RTD]      SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1
                           }
                        } else {
         // A(:,q) is zero column
                           if (ir1 == 0) NOTROT = NOTROT + 1;
                           PSKIPPED = PSKIPPED + 1
                        }

                        if ( ( i.LE.SWBAND ) && ( PSKIPPED.GT.ROWSKIP ) ) {
                           if (ir1 == 0) AAPP = -AAPP;
                           NOTROT = 0
                           GO TO 2103
                        }

                     } // 2002
      // END q-LOOP

                     } // 2103
      // bailed out of q-loop

                     SVA( p ) = AAPP

                  } else {
                     SVA( p ) = AAPP
                     IF( ( ir1 == 0 ) && ( AAPP == ZERO ) ) NOTROT = NOTROT + MIN( igl+KBL-1, N ) - p
                  }

               } // 2001
      // end of the p-loop
      // end of doing the block ( ibr, ibr )
            } // 1002
      // end of ir1-loop

* ... go to the off diagonal blocks

            igl = ( ibr-1 )*KBL + 1

            for (jbc = ibr + 1; jbc <= NBL; jbc++) { // 2010

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
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 scopy(M, A( 1, p ), 1, WORK( N+1 ), 1 );
                                 slascl('G', 0, 0, AAPP, WORK( p ), M, 1, WORK( N+1 ), LDA, IERR );
                                 AAPQ = SDOT( M, WORK( N+1 ), 1, A( 1, q ), 1 )*WORK( q ) / AAQQ
                              }
                           } else {
                              if ( AAPP.GE.AAQQ ) {
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              } else {
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              }
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( SDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 scopy(M, A( 1, q ), 1, WORK( N+1 ), 1 );
                                 slascl('G', 0, 0, AAQQ, WORK( q ), M, 1, WORK( N+1 ), LDA, IERR );
                                 AAPQ = SDOT( M, WORK( N+1 ), 1, A( 1, p ), 1 )*WORK( p ) / AAPP
                              }
                           }

                           MXAAPQ = MAX( MXAAPQ, ABS( AAPQ ) )

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( ABS( AAPQ ).GT.TOL ) {
                              NOTROT = 0
*[RTD]      ROTATED  = ROTATED + 1
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*ABS( AQOAP-APOAQ ) / AAPQ
                                 if (AAQQ.GT.AAPP0) THETA = -THETA;

                                 if ( ABS( THETA ).GT.BIGTHETA ) {
                                    T = HALF / THETA
                                    FASTR( 3 ) = T*WORK( p ) / WORK( q )
                                    FASTR( 4 ) = -T*WORK( q ) / WORK( p );
                                    srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, ABS( T ) )
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -SIGN( ONE, AAPQ )
                                    if (AAQQ.GT.AAPP0) THSIGN = -THSIGN;
                                    T = ONE / ( THETA+THSIGN* SQRT( ONE+THETA*THETA ) )
                                    CS = SQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
                                    MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )

                                    APOAQ = WORK( p ) / WORK( q )
                                    AQOAP = WORK( q ) / WORK( p )
                                    if ( WORK( p ).GE.ONE ) {

                                       if ( WORK( q ).GE.ONE ) {
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q )*CS
                                          srotm(M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL SROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR );
                                       } else {
                                          saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          if ( RSVEC ) {
                                             saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                          }
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q ) / CS
                                       }
                                    } else {
                                       if ( WORK( q ).GE.ONE ) {
                                          saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                          saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                          if ( RSVEC ) {
                                             saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             saxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                          }
                                          WORK( p ) = WORK( p ) / CS
                                          WORK( q ) = WORK( q )*CS
                                       } else {
                                          IF( WORK( p ).GE.WORK( q ) ) THEN;
                                             saxpy(M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             saxpy(M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             WORK( p ) = WORK( p )*CS
                                             WORK( q ) = WORK( q ) / CS
                                             if ( RSVEC ) {
                                                saxpy(MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                                saxpy(MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                             }
                                          } else {
                                             saxpy(M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 );
                                             saxpy(M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 );
                                             WORK( p ) = WORK( p ) / CS
                                             WORK( q ) = WORK( q )*CS
                                             if ( RSVEC ) {
                                                saxpy(MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 );
                                                saxpy(MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 );
                                             }
                                          }
                                       }
                                    }
                                 }

                              } else {
                                 if ( AAPP.GT.AAQQ ) {
                                    scopy(M, A( 1, p ), 1, WORK( N+1 ), 1 );
                                    slascl('G', 0, 0, AAPP, ONE, M, 1, WORK( N+1 ), LDA, IERR );
                                    slascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR );
                                    TEMP1 = -AAPQ*WORK( p ) / WORK( q )
                                    saxpy(M, TEMP1, WORK( N+1 ), 1, A( 1, q ), 1 );
                                    slascl('G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 } else {
                                    scopy(M, A( 1, q ), 1, WORK( N+1 ), 1 );
                                    slascl('G', 0, 0, AAQQ, ONE, M, 1, WORK( N+1 ), LDA, IERR );
                                    slascl('G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR );
                                    TEMP1 = -AAPQ*WORK( q ) / WORK( p )
                                    saxpy(M, TEMP1, WORK( N+1 ), 1, A( 1, p ), 1 );
                                    slascl('G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
                                    SVA( p ) = AAPP*SQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q)
            // .. recompute SVA(q)
                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) && ( AAQQ.GT.ROOTSFMIN ) ) THEN                                     SVA( q ) = SNRM2( M, A( 1, q ), 1 )* WORK( q )
                                 } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    slassq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ )*WORK( q )
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) {
                                 IF( ( AAPP.LT.ROOTBIG ) && ( AAPP.GT.ROOTSFMIN ) ) THEN                                     AAPP = SNRM2( M, A( 1, p ), 1 )* WORK( p )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    slassq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*SQRT( AAPP )*WORK( p )
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

                        if ( ( i.LE.SWBAND ) && ( IJBLSK.GE.BLSKIP ) ) {
                           SVA( p ) = AAPP
                           NOTROT = 0
                           GO TO 2011
                        }
                        if ( ( i.LE.SWBAND ) && ( PSKIPPED.GT.ROWSKIP ) ) {
                           AAPP = -AAPP
                           NOTROT = 0
                           GO TO 2203
                        }

                     } // 2200
         // end of the q-loop
                     } // 2203

                     SVA( p ) = AAPP

                  } else {

                     if (AAPP == ZERO) NOTROT = NOTROT + MIN( jgl+KBL-1, N ) - jgl + 1;
                     if (AAPP.LT.ZERO) NOTROT = 0;

                  }

               } // 2100
      // end of the p-loop
            } // 2010
      // end of the jbc-loop
            } // 2011
*2011 bailed out of the jbc-loop
            DO 2012 p = igl, MIN( igl+KBL-1, N )
               SVA( p ) = ABS( SVA( p ) )
            } // 2012
***
         } // 2000
*2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ).LT.ROOTBIG ) && ( SVA( N ).GT.ROOTSFMIN ) ) {
            SVA( N ) = SNRM2( M, A( 1, N ), 1 )*WORK( N )
         } else {
            T = ZERO
            AAPP = ONE
            slassq(M, A( 1, N ), 1, T, AAPP );
            SVA( N ) = T*SQRT( AAPP )*WORK( N )
         }

      // Additional steering devices

         IF( ( i.LT.SWBAND ) && ( ( MXAAPQ.LE.ROOTTOL ) .OR. ( ISWROT.LE.N ) ) )SWBAND = i

         if ( ( i.GT.SWBAND+1 ) && ( MXAAPQ.LT.SQRT( FLOAT( N ) )* TOL ) && ( FLOAT( N )*MXAAPQ*MXSINJ.LT.TOL ) ) {
            GO TO 1994
         }

         if (NOTROT.GE.EMPTSW) GO TO 1994;

      } // 1993
      // end i=1:NSWEEP loop

* #:( Reaching this point means that the procedure has not converged.
      INFO = NSWEEP - 1
      GO TO 1995

      } // 1994
* #:) Reaching this point means numerical convergence after the i-th
      // sweep.

      INFO = 0
* #:) INFO = 0 confirms successful iterations.
      } // 1995

      // Sort the singular values and find how many are above
      // the underflow threshold.

      N2 = 0
      N4 = 0
      for (p = 1; p <= N - 1; p++) { // 5991
         q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         if ( p != q ) {
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            TEMP1 = WORK( p )
            WORK( p ) = WORK( q )
            WORK( q ) = TEMP1
            sswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            if (RSVEC) CALL SSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 );
         }
         if ( SVA( p ) != ZERO ) {
            N4 = N4 + 1
            IF( SVA( p )*SKL.GT.SFMIN )N2 = N2 + 1
         }
      } // 5991
      if ( SVA( N ) != ZERO ) {
         N4 = N4 + 1
         IF( SVA( N )*SKL.GT.SFMIN )N2 = N2 + 1
      }

      // Normalize the left singular vectors.

      if ( LSVEC .OR. UCTOL ) {
         for (p = 1; p <= N2; p++) { // 1998
            sscal(M, WORK( p ) / SVA( p ), A( 1, p ), 1 );
         } // 1998
      }

      // Scale the product of Jacobi rotations (assemble the fast rotations).

      if ( RSVEC ) {
         if ( APPLV ) {
            for (p = 1; p <= N; p++) { // 2398
               sscal(MVL, WORK( p ), V( 1, p ), 1 );
            } // 2398
         } else {
            for (p = 1; p <= N; p++) { // 2399
               TEMP1 = ONE / SNRM2( MVL, V( 1, p ), 1 )
               sscal(MVL, TEMP1, V( 1, p ), 1 );
            } // 2399
         }
      }

      // Undo scaling, if necessary (and possible).
      if ( ( ( SKL.GT.ONE ) && ( SVA( 1 ).LT.( BIG / SKL ) ) ) .OR. ( ( SKL.LT.ONE ) && ( SVA( MAX( N2, 1 ) ) .GT. ( SFMIN / SKL ) ) ) ) {
         for (p = 1; p <= N; p++) { // 2400
            SVA( P ) = SKL*SVA( P )
         } // 2400
         SKL = ONE
      }

      WORK( 1 ) = SKL
      // The singular values of A are SKL*SVA(1:N). If SKL != ONE
      // then some of the singular values may overflow or underflow and
      // the spectrum is given in this factored representation.

      WORK( 2 ) = FLOAT( N4 )
      // N4 is the number of computed nonzero singular values of A.

      WORK( 3 ) = FLOAT( N2 )
      // N2 is the number of singular values of A greater than SFMIN.
      // If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
      // that may carry some information.

      WORK( 4 ) = FLOAT( i )
      // i is the index of the last sweep before declaring convergence.

      WORK( 5 ) = MXAAPQ
      // MXAAPQ is the largest absolute value of scaled pivots in the
      // last sweep

      WORK( 6 ) = MXSINJ
      // MXSINJ is the largest absolute value of the sines of Jacobi angles
      // in the last sweep

      RETURN
      // ..
      // .. END OF SGESVJ
      // ..
      }
