      SUBROUTINE CGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, LDV, CWORK, LWORK, RWORK, LRWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE
      // .. Scalar Arguments ..
      int                INFO, LDA, LDV, LWORK, LRWORK, M, MV, N;
      String             JOBA, JOBU, JOBV;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ),  V( LDV, * ), CWORK( LWORK )
      REAL               RWORK( LRWORK ), SVA( N )
      // ..

*  =====================================================================

      // .. Local Parameters ..
      REAL         ZERO,         HALF,         ONE
      const      ZERO = 0.0E0, HALF = 0.5E0, ONE = 1.0E0;
      COMPLEX      CZERO,                  CONE
      const      CZERO = (0.0E0, 0.0E0), CONE = (1.0E0, 0.0E0) ;
      int          NSWEEP;
      const      NSWEEP = 30 ;
      // ..
      // .. Local Scalars ..
      COMPLEX    AAPQ, OMPQ
      REAL       AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, CTOL, EPSLN, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SKL, SFMIN, SMALL, SN, T, TEMP1, THETA, THSIGN, TOL;
      int        BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, N2, N34, N4, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND, MINMN, LWMIN, LRWMIN;
      bool       APPLV, GOSCALE, LOWER, LQUERY, LSVEC, NOSCALE, ROTOK, RSVEC, UCTOL, UPPER;
      // ..
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, CONJG, REAL, SIGN, SQRT
      // ..
      // .. External Functions ..
      // ..
      // from BLAS
      REAL               SCNRM2
      COMPLEX            CDOTC
      // EXTERNAL CDOTC, SCNRM2
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
      // EXTERNAL CCOPY, CROT, CSSCAL, CSWAP, CAXPY
      // from LAPACK
      // EXTERNAL CLASCL, CLASET, CLASSQ, SLASCL, XERBLA
      // EXTERNAL CGSVJ0, CGSVJ1
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      LSVEC = LSAME( JOBU, 'U' ) .OR. LSAME( JOBU, 'F' )
      UCTOL = LSAME( JOBU, 'C' )
      RSVEC = LSAME( JOBV, 'V' ) .OR. LSAME( JOBV, 'J' )
      APPLV = LSAME( JOBV, 'A' )
      UPPER = LSAME( JOBA, 'U' )
      LOWER = LSAME( JOBA, 'L' )

      MINMN = MIN( M, N )
      if ( MINMN.EQ.0 ) {
         LWMIN  = 1
         LRWMIN = 1
      } else {
         LWMIN  = M + N
         LRWMIN = MAX( 6, N )
      }

      LQUERY = ( LWORK.EQ.-1 ) .OR. ( LRWORK.EQ.-1 )
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
      } else if ( ( RSVEC .AND. ( LDV.LT.N ) ) .OR. ( APPLV .AND. ( LDV.LT.MV ) ) ) {
         INFO = -11
      } else if ( UCTOL .AND. ( RWORK( 1 ).LE.ONE ) ) {
         INFO = -12
      } else if ( LWORK.LT.LWMIN .AND. ( .NOT.LQUERY ) ) {
         INFO = -13
      } else if ( LRWORK.LT.LRWMIN .AND. ( .NOT.LQUERY ) ) {
         INFO = -15
      } else {
         INFO = 0
      }

      // #:(
      if ( INFO.NE.0 ) {
         xerbla('CGESVJ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         CWORK( 1 ) = SROUNDUP_LWORK( LWMIN )
         RWORK( 1 ) = SROUNDUP_LWORK( LRWMIN )
         RETURN
      }

* #:) Quick return for void matrix

      IF( MINMN.EQ.0 ) RETURN

      // Set numerical parameters
      // The stopping criterion for Jacobi rotations is

      // max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS

      // where EPS is the round-off and CTOL is defined as follows:

      if ( UCTOL ) {
         // ... user controlled
         CTOL = RWORK( 1 )
      } else {
         // ... default
         if ( LSVEC .OR. RSVEC .OR. APPLV ) {
            CTOL = SQRT( REAL( M ) )
         } else {
            CTOL = REAL( M )
         }
      }
      // ... and the machine dependent parameters are
*[!]  (Make sure that SLAMCH() works properly on the target machine.)

      EPSLN = SLAMCH( 'Epsilon' )
      ROOTEPS = SQRT( EPSLN )
      SFMIN = SLAMCH( 'SafeMinimum' )
      ROOTSFMIN = SQRT( SFMIN )
      SMALL = SFMIN / EPSLN
       // BIG = SLAMCH( 'Overflow' )
      BIG     = ONE  / SFMIN
      ROOTBIG = ONE / ROOTSFMIN
      // LARGE = BIG / SQRT( REAL( M*N ) )
      BIGTHETA = ONE / ROOTEPS

      TOL = CTOL*EPSLN
      ROOTTOL = SQRT( TOL )

      if ( REAL( M )*EPSLN.GE.ONE ) {
         INFO = -4
         xerbla('CGESVJ', -INFO );
         RETURN
      }

      // Initialize the right singular vector matrix.

      if ( RSVEC ) {
         MVL = N
         claset('A', MVL, N, CZERO, CONE, V, LDV );
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

      SKL = ONE / SQRT( REAL( M )*REAL( N ) )
      NOSCALE = .TRUE.
      GOSCALE = .TRUE.

      if ( LOWER ) {
         // the input matrix is M-by-N lower triangular (trapezoidal)
         DO 1874 p = 1, N
            AAPP = ZERO
            AAQQ = ONE
            classq(M-p+1, A( p, p ), 1, AAPP, AAQQ );
            if ( AAPP.GT.BIG ) {
               INFO = -6
               xerbla('CGESVJ', -INFO );
               RETURN
            }
            AAQQ = SQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) .AND. NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = .FALSE.
               SVA( p ) = AAPP*( AAQQ*SKL )
               if ( GOSCALE ) {
                  GOSCALE = .FALSE.
                  DO 1873 q = 1, p - 1
                     SVA( q ) = SVA( q )*SKL
 1873             CONTINUE
               }
            }
 1874    CONTINUE
      } else if ( UPPER ) {
         // the input matrix is M-by-N upper triangular (trapezoidal)
         DO 2874 p = 1, N
            AAPP = ZERO
            AAQQ = ONE
            classq(p, A( 1, p ), 1, AAPP, AAQQ );
            if ( AAPP.GT.BIG ) {
               INFO = -6
               xerbla('CGESVJ', -INFO );
               RETURN
            }
            AAQQ = SQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) .AND. NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = .FALSE.
               SVA( p ) = AAPP*( AAQQ*SKL )
               if ( GOSCALE ) {
                  GOSCALE = .FALSE.
                  DO 2873 q = 1, p - 1
                     SVA( q ) = SVA( q )*SKL
 2873             CONTINUE
               }
            }
 2874    CONTINUE
      } else {
         // the input matrix is M-by-N general dense
         DO 3874 p = 1, N
            AAPP = ZERO
            AAQQ = ONE
            classq(M, A( 1, p ), 1, AAPP, AAQQ );
            if ( AAPP.GT.BIG ) {
               INFO = -6
               xerbla('CGESVJ', -INFO );
               RETURN
            }
            AAQQ = SQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) .AND. NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = .FALSE.
               SVA( p ) = AAPP*( AAQQ*SKL )
               if ( GOSCALE ) {
                  GOSCALE = .FALSE.
                  DO 3873 q = 1, p - 1
                     SVA( q ) = SVA( q )*SKL
 3873             CONTINUE
               }
            }
 3874    CONTINUE
      }

      IF( NOSCALE )SKL = ONE

      // Move the smaller part of the spectrum from the underflow threshold
*(!)  Start by determining the position of the nonzero entries of the
      // array SVA() relative to ( SFMIN, BIG ).

      AAPP = ZERO
      AAQQ = BIG
      DO 4781 p = 1, N
         IF( SVA( p ).NE.ZERO )AAQQ = MIN( AAQQ, SVA( p ) )
         AAPP = MAX( AAPP, SVA( p ) )
 4781 CONTINUE

* #:) Quick return for zero matrix

      if ( AAPP.EQ.ZERO ) {
         IF( LSVEC )CALL CLASET( 'G', M, N, CZERO, CONE, A, LDA )
         RWORK( 1 ) = ONE
         RWORK( 2 ) = ZERO
         RWORK( 3 ) = ZERO
         RWORK( 4 ) = ZERO
         RWORK( 5 ) = ZERO
         RWORK( 6 ) = ZERO
         RETURN
      }

* #:) Quick return for one-column matrix

      if ( N.EQ.1 ) {
         IF( LSVEC )CALL CLASCL( 'G', 0, 0, SVA( 1 ), SKL, M, 1, A( 1, 1 ), LDA, IERR )
         RWORK( 1 ) = ONE / SKL
         if ( SVA( 1 ).GE.SFMIN ) {
            RWORK( 2 ) = ONE
         } else {
            RWORK( 2 ) = ZERO
         }
         RWORK( 3 ) = ZERO
         RWORK( 4 ) = ZERO
         RWORK( 5 ) = ZERO
         RWORK( 6 ) = ZERO
         RETURN
      }

      // Protect small singular values from underflow, and try to
      // avoid underflows/overflows in computing Jacobi rotations.

      SN = SQRT( SFMIN / EPSLN )
      TEMP1 = SQRT( BIG / REAL( N ) )
      if ( ( AAPP.LE.SN ) .OR. ( AAQQ.GE.TEMP1 ) .OR. ( ( SN.LE.AAQQ ) .AND. ( AAPP.LE.TEMP1 ) ) ) {
         TEMP1 = MIN( BIG, TEMP1 / AAPP )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.LE.SN ) .AND. ( AAPP.LE.TEMP1 ) ) {
         TEMP1 = MIN( SN / AAQQ, BIG / ( AAPP*SQRT( REAL( N ) ) ) )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.GE.SN ) .AND. ( AAPP.GE.TEMP1 ) ) {
         TEMP1 = MAX( SN / AAQQ, TEMP1 / AAPP )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.LE.SN ) .AND. ( AAPP.GE.TEMP1 ) ) {
         TEMP1 = MIN( SN / AAQQ, BIG / ( SQRT( REAL( N ) )*AAPP ) )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else {
         TEMP1 = ONE
      }

      // Scale, if necessary

      if ( TEMP1.NE.ONE ) {
         slascl('G', 0, 0, ONE, TEMP1, N, 1, SVA, N, IERR );
      }
      SKL = TEMP1*SKL
      if ( SKL.NE.ONE ) {
         clascl(JOBA, 0, 0, ONE, SKL, M, N, A, LDA, IERR );
         SKL = ONE / SKL
      }

      // Row-cyclic Jacobi SVD algorithm with column pivoting

      EMPTSW = ( N*( N-1 ) ) / 2
      NOTROT = 0

      DO 1868 q = 1, N
         CWORK( q ) = CONE
 1868 CONTINUE



      SWBAND = 3
*[TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
      // if CGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm CGEJSV. For sweeps i=1:SWBAND the procedure
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

      if ( ( LOWER .OR. UPPER ) .AND. ( N.GT.MAX( 64, 4*KBL ) ) ) {
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

            cgsvj0(JOBV, M-N34, N-N34, A( N34+1, N34+1 ), LDA, CWORK( N34+1 ), SVA( N34+1 ), MVL, V( N34*q+1, N34+1 ), LDV, EPSLN, SFMIN, TOL, 2, CWORK( N+1 ), LWORK-N, IERR )              CALL CGSVJ0( JOBV, M-N2, N34-N2, A( N2+1, N2+1 ), LDA, CWORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 2, CWORK( N+1 ), LWORK-N, IERR )              CALL CGSVJ1( JOBV, M-N2, N-N2, N4, A( N2+1, N2+1 ), LDA, CWORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            cgsvj0(JOBV, M-N4, N2-N4, A( N4+1, N4+1 ), LDA, CWORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            cgsvj0(JOBV, M, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            cgsvj1(JOBV, M, N2, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );


         } else if ( UPPER ) {


            cgsvj0(JOBV, N4, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 2, CWORK( N+1 ), LWORK-N, IERR );

            cgsvj0(JOBV, N2, N4, A( 1, N4+1 ), LDA, CWORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            cgsvj1(JOBV, N2, N2, N4, A, LDA, CWORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

            cgsvj0(JOBV, N2+N4, N4, A( 1, N2+1 ), LDA, CWORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, CWORK( N+1 ), LWORK-N, IERR );

         }

      }

      // .. Row-cyclic pivot strategy with de Rijk's pivoting ..

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

         DO 2000 ibr = 1, NBL

            igl = ( ibr-1 )*KBL + 1

            DO 1002 ir1 = 0, MIN( LKAHEAD, NBL-ibr )

               igl = igl + ir1*KBL

               DO 2001 p = igl, MIN( igl+KBL-1, N-1 )

      // .. de Rijk's pivoting

                  q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
                  if ( p.NE.q ) {
                     cswap(M, A( 1, p ), 1, A( 1, q ), 1 );
                     IF( RSVEC )CALL CSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
                     TEMP1 = SVA( p )
                     SVA( p ) = SVA( q )
                     SVA( q ) = TEMP1
                     AAPQ = CWORK(p)
                     CWORK(p) = CWORK(q)
                     CWORK(q) = AAPQ
                  }

                  if ( ir1.EQ.0 ) {

         // Column norms are periodically updated by explicit
         // norm computation.
*[!]     Caveat:
         // Unfortunately, some BLAS implementations compute SCNRM2(M,A(1,p),1)
         // as SQRT(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to
         // overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to
         // underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
         // Hence, SCNRM2 cannot be trusted, not even in the case when
         // the true norm is far from the under(over)flow boundaries.
         // If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF
         // below should be replaced with "AAPP = SCNRM2( M, A(1,p), 1 )".

                     if ( ( SVA( p ).LT.ROOTBIG ) .AND. ( SVA( p ).GT.ROOTSFMIN ) ) {
                        SVA( p ) = SCNRM2( M, A( 1, p ), 1 )
                     } else {
                        TEMP1 = ZERO
                        AAPP = ONE
                        classq(M, A( 1, p ), 1, TEMP1, AAPP );
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
                                 AAPQ = ( CDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP
                              } else {
                                 CALL CCOPY( M, A( 1, p ), 1, CWORK(N+1), 1 )                                  CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1), LDA, IERR )                                  AAPQ = CDOTC( M, CWORK(N+1), 1, A( 1, q ), 1 ) / AAQQ
                              }
                           } else {
                              ROTOK = AAPP.LE.( AAQQ / SMALL )
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( CDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAPP ) / AAQQ
                              } else {
                                 ccopy(M, A( 1, q ), 1, CWORK(N+1), 1 )                                  CALL CLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 AAPQ = CDOTC( M, A(1, p ), 1, CWORK(N+1), 1 ) / AAPP
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
                                     crot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        crot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
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

                                    crot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        crot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 CWORK(p) = -CWORK(q) * OMPQ

                                 } else {
               // .. have to use modified Gram-Schmidt like transformation
                                 ccopy(M, A( 1, p ), 1, CWORK(N+1), 1 )                                  CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 clascl('G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )                                  CALL CAXPY( M, -AAPQ, CWORK(N+1), 1, A( 1, q ), 1 )                                  CALL CLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )                                  SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) );
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // recompute SVA(q), SVA(p).

                              if ( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) {
                                    SVA( q ) = SCNRM2( M, A( 1, q ), 1 )
                                 } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    classq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ )
                                 }
                              }
                              if ( ( AAPP / AAPP0 ).LE.ROOTEPS ) {
                                 if ( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) {
                                    AAPP = SCNRM2( M, A( 1, p ), 1 )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    classq(M, A( 1, p ), 1, T, AAPP );
                                    AAPP = T*SQRT( AAPP )
                                 }
                                 SVA( p ) = AAPP
                              }

                           } else {
                              // A(:,p) and A(:,q) already numerically orthogonal
                              IF( ir1.EQ.0 )NOTROT = NOTROT + 1
*[RTD]      SKIPPED  = SKIPPED + 1
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
                                 AAPQ = ( CDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP
                              } else {
                                 ccopy(M, A( 1, p ), 1, CWORK(N+1), 1 )                                  CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 AAPQ = CDOTC( M, CWORK(N+1), 1, A( 1, q ), 1 ) / AAQQ
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
                                 ccopy(M, A( 1, q ), 1, CWORK(N+1), 1 )                                  CALL CLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, CWORK(N+1), LDA, IERR );
                                 AAPQ = CDOTC( M, A( 1, p ), 1, CWORK(N+1),  1 ) / AAPP
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
                                    crot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T );
                                    if ( RSVEC ) {
                                        crot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T );
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

                                    crot(M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    if ( RSVEC ) {
                                        crot(MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN );
                                    }
                                 }
                                 CWORK(p) = -CWORK(q) * OMPQ

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                               if ( AAPP.GT.AAQQ ) {
                                    ccopy(M, A( 1, p ), 1, CWORK(N+1), 1 )                                     CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, CWORK(N+1),LDA, IERR )                                     CALL CLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )                                     CALL CAXPY( M, -AAPQ, CWORK(N+1), 1, A( 1, q ), 1 )                                     CALL CLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR );
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                               } else {
                                   ccopy(M, A( 1, q ), 1, CWORK(N+1), 1 )                                     CALL CLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, CWORK(N+1),LDA, IERR )                                     CALL CLASCL( 'G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR )                                     CALL CAXPY( M, -CONJG(AAPQ), CWORK(N+1), 1, A( 1, p ), 1 )                                     CALL CLASCL( 'G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR );
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
                                    classq(M, A( 1, q ), 1, T, AAQQ );
                                    SVA( q ) = T*SQRT( AAQQ )
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) {
                                 if ( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) {
                                    AAPP = SCNRM2( M, A( 1, p ), 1 )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    classq(M, A( 1, p ), 1, T, AAPP );
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
            classq(M, A( 1, N ), 1, T, AAPP );
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

      // Sort the singular values and find how many are above
      // the underflow threshold.

      N2 = 0
      N4 = 0
      DO 5991 p = 1, N - 1
         q = ISAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         if ( p.NE.q ) {
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            cswap(M, A( 1, p ), 1, A( 1, q ), 1 );
            IF( RSVEC )CALL CSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
         }
         if ( SVA( p ).NE.ZERO ) {
            N4 = N4 + 1
            IF( SVA( p )*SKL.GT.SFMIN )N2 = N2 + 1
         }
 5991 CONTINUE
      if ( SVA( N ).NE.ZERO ) {
         N4 = N4 + 1
         IF( SVA( N )*SKL.GT.SFMIN )N2 = N2 + 1
      }

      // Normalize the left singular vectors.

      if ( LSVEC .OR. UCTOL ) {
         DO 1998 p = 1, N4
            // CALL CSSCAL( M, ONE / SVA( p ), A( 1, p ), 1 )
            clascl('G',0,0, SVA(p), ONE, M, 1, A(1,p), M, IERR );
 1998    CONTINUE
      }

      // Scale the product of Jacobi rotations.

      if ( RSVEC ) {
            DO 2399 p = 1, N
               TEMP1 = ONE / SCNRM2( MVL, V( 1, p ), 1 )
               csscal(MVL, TEMP1, V( 1, p ), 1 );
 2399       CONTINUE
      }

      // Undo scaling, if necessary (and possible).
      if ( ( ( SKL.GT.ONE ) .AND. ( SVA( 1 ).LT.( BIG / SKL ) ) ) .OR. ( ( SKL.LT.ONE ) .AND. ( SVA( MAX( N2, 1 ) ) .GT. ( SFMIN / SKL ) ) ) ) {
         DO 2400 p = 1, N
            SVA( P ) = SKL*SVA( P )
 2400    CONTINUE
         SKL = ONE
      }

      RWORK( 1 ) = SKL
      // The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE
      // then some of the singular values may overflow or underflow and
      // the spectrum is given in this factored representation.

      RWORK( 2 ) = REAL( N4 )
      // N4 is the number of computed nonzero singular values of A.

      RWORK( 3 ) = REAL( N2 )
      // N2 is the number of singular values of A greater than SFMIN.
      // If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
      // that may carry some information.

      RWORK( 4 ) = REAL( i )
      // i is the index of the last sweep before declaring convergence.

      RWORK( 5 ) = MXAAPQ
      // MXAAPQ is the largest absolute value of scaled pivots in the
      // last sweep

      RWORK( 6 ) = MXSINJ
      // MXSINJ is the largest absolute value of the sines of Jacobi angles
      // in the last sweep

      RETURN
      // ..
      // .. END OF CGESVJ
      // ..
      }
