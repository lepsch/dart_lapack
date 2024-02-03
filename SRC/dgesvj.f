      SUBROUTINE DGESVJ( JOBA, JOBU, JOBV, M, N, A, LDA, SVA, MV, V, LDV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDV, LWORK, M, MV, N;
      String             JOBA, JOBU, JOBV;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), SVA( N ), V( LDV, * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Local Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0;
      int                NSWEEP;
      const              NSWEEP = 30 ;
      // ..
      // .. Local Scalars ..
      double             AAPP, AAPP0, AAPQ, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, CTOL, EPSLN, LARGE, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SKL, SFMIN, SMALL, SN, T, TEMP1, THETA, THSIGN, TOL;
      int                BLSKIP, EMPTSW, i, ibr, IERR, igl, IJBLSK, ir1, ISWROT, jbc, jgl, KBL, LKAHEAD, MVL, N2, N34, N4, NBL, NOTROT, p, PSKIPPED, q, ROWSKIP, SWBAND, MINMN, LWMIN;
      bool               APPLV, GOSCALE, LOWER, LQUERY, LSVEC, NOSCALE, ROTOK, RSVEC, UCTOL, UPPER;
      // ..
      // .. Local Arrays ..
      double             FASTR( 5 );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DABS, MAX, MIN, DBLE, DSIGN, DSQRT
      // ..
      // .. External Functions ..
      // ..
      // from BLAS
      double             DDOT, DNRM2;
      // EXTERNAL DDOT, DNRM2
      int                IDAMAX;
      // EXTERNAL IDAMAX
      // from LAPACK
      double             DLAMCH;
      // EXTERNAL DLAMCH
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // ..
      // from BLAS
      // EXTERNAL DAXPY, DCOPY, DROTM, DSCAL, DSWAP
      // from LAPACK
      // EXTERNAL DLASCL, DLASET, DLASSQ, XERBLA

      // EXTERNAL DGSVJ0, DGSVJ1
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
      if ( MINMN.EQ.0 ) {
         LWMIN = 1
      } else {
         LWMIN = MAX( 6, M+N )
      }

      LQUERY = ( LWORK.EQ.-1 )
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
      } else if ( UCTOL .AND. ( WORK( 1 ).LE.ONE ) ) {
         INFO = -12
      } else if ( LWORK.LT.LWMIN .AND. ( .NOT.LQUERY ) ) {
         INFO = -13
      } else {
         INFO = 0
      }

      // #:(
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DGESVJ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = LWMIN
         RETURN
      }

* #:) Quick return for void matrix

      IF( MINMN.EQ.0 ) RETURN

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
            CTOL = DSQRT( DBLE( M ) )
         } else {
            CTOL = DBLE( M )
         }
      }
      // ... and the machine dependent parameters are
*[!]  (Make sure that DLAMCH() works properly on the target machine.)

      EPSLN = DLAMCH( 'Epsilon' )
      ROOTEPS = DSQRT( EPSLN )
      SFMIN = DLAMCH( 'SafeMinimum' )
      ROOTSFMIN = DSQRT( SFMIN )
      SMALL = SFMIN / EPSLN
      BIG = DLAMCH( 'Overflow' )
      // BIG         = ONE    / SFMIN
      ROOTBIG = ONE / ROOTSFMIN
      LARGE = BIG / DSQRT( DBLE( M*N ) )
      BIGTHETA = ONE / ROOTEPS

      TOL = CTOL*EPSLN
      ROOTTOL = DSQRT( TOL )

      if ( DBLE( M )*EPSLN.GE.ONE ) {
         INFO = -4
         CALL XERBLA( 'DGESVJ', -INFO )
         RETURN
      }

      // Initialize the right singular vector matrix.

      if ( RSVEC ) {
         MVL = N
         CALL DLASET( 'A', MVL, N, ZERO, ONE, V, LDV )
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
      // DSQRT(N)*max_i SVA(i) does not overflow. If INFinite entries
      // in A are detected, the procedure returns with INFO=-6.

      SKL= ONE / DSQRT( DBLE( M )*DBLE( N ) )
      NOSCALE = .TRUE.
      GOSCALE = .TRUE.

      if ( LOWER ) {
        t // he input matrix is M-by-N lower triangular (trapezoidal)
         DO 1874 p = 1, N
            AAPP = ZERO
            AAQQ = ONE
            CALL DLASSQ( M-p+1, A( p, p ), 1, AAPP, AAQQ )
            if ( AAPP.GT.BIG ) {
               INFO = -6
               CALL XERBLA( 'DGESVJ', -INFO )
               RETURN
            }
            AAQQ = DSQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) .AND. NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = .FALSE.
               SVA( p ) = AAPP*( AAQQ*SKL)
               if ( GOSCALE ) {
                  GOSCALE = .FALSE.
                  DO 1873 q = 1, p - 1
                     SVA( q ) = SVA( q )*SKL
 1873             CONTINUE
               }
            }
 1874    CONTINUE
      } else if ( UPPER ) {
        t // he input matrix is M-by-N upper triangular (trapezoidal)
         DO 2874 p = 1, N
            AAPP = ZERO
            AAQQ = ONE
            CALL DLASSQ( p, A( 1, p ), 1, AAPP, AAQQ )
            if ( AAPP.GT.BIG ) {
               INFO = -6
               CALL XERBLA( 'DGESVJ', -INFO )
               RETURN
            }
            AAQQ = DSQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) .AND. NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = .FALSE.
               SVA( p ) = AAPP*( AAQQ*SKL)
               if ( GOSCALE ) {
                  GOSCALE = .FALSE.
                  DO 2873 q = 1, p - 1
                     SVA( q ) = SVA( q )*SKL
 2873             CONTINUE
               }
            }
 2874    CONTINUE
      } else {
        t // he input matrix is M-by-N general dense
         DO 3874 p = 1, N
            AAPP = ZERO
            AAQQ = ONE
            CALL DLASSQ( M, A( 1, p ), 1, AAPP, AAQQ )
            if ( AAPP.GT.BIG ) {
               INFO = -6
               CALL XERBLA( 'DGESVJ', -INFO )
               RETURN
            }
            AAQQ = DSQRT( AAQQ )
            if ( ( AAPP.LT.( BIG / AAQQ ) ) .AND. NOSCALE ) {
               SVA( p ) = AAPP*AAQQ
            } else {
               NOSCALE = .FALSE.
               SVA( p ) = AAPP*( AAQQ*SKL)
               if ( GOSCALE ) {
                  GOSCALE = .FALSE.
                  DO 3873 q = 1, p - 1
                     SVA( q ) = SVA( q )*SKL
 3873             CONTINUE
               }
            }
 3874    CONTINUE
      }

      IF( NOSCALE )SKL= ONE

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
         IF( LSVEC )CALL DLASET( 'G', M, N, ZERO, ONE, A, LDA )
         WORK( 1 ) = ONE
         WORK( 2 ) = ZERO
         WORK( 3 ) = ZERO
         WORK( 4 ) = ZERO
         WORK( 5 ) = ZERO
         WORK( 6 ) = ZERO
         RETURN
      }

* #:) Quick return for one-column matrix

      if ( N.EQ.1 ) {
         IF( LSVEC )CALL DLASCL( 'G', 0, 0, SVA( 1 ), SKL, M, 1, A( 1, 1 ), LDA, IERR )
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

      SN = DSQRT( SFMIN / EPSLN )
      TEMP1 = DSQRT( BIG / DBLE( N ) )
      if ( ( AAPP.LE.SN ) .OR. ( AAQQ.GE.TEMP1 ) .OR. ( ( SN.LE.AAQQ ) .AND. ( AAPP.LE.TEMP1 ) ) ) {
         TEMP1 = MIN( BIG, TEMP1 / AAPP )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.LE.SN ) .AND. ( AAPP.LE.TEMP1 ) ) {
         TEMP1 = MIN( SN / AAQQ, BIG / ( AAPP*DSQRT( DBLE( N ) ) ) )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.GE.SN ) .AND. ( AAPP.GE.TEMP1 ) ) {
         TEMP1 = MAX( SN / AAQQ, TEMP1 / AAPP )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else if ( ( AAQQ.LE.SN ) .AND. ( AAPP.GE.TEMP1 ) ) {
         TEMP1 = MIN( SN / AAQQ, BIG / ( DSQRT( DBLE( N ) )*AAPP ) )
          // AAQQ  = AAQQ*TEMP1
          // AAPP  = AAPP*TEMP1
      } else {
         TEMP1 = ONE
      }

      // Scale, if necessary

      if ( TEMP1.NE.ONE ) {
         CALL DLASCL( 'G', 0, 0, ONE, TEMP1, N, 1, SVA, N, IERR )
      }
      SKL= TEMP1*SKL
      if ( SKL.NE.ONE ) {
         CALL DLASCL( JOBA, 0, 0, ONE, SKL, M, N, A, LDA, IERR )
         SKL= ONE / SKL
      }

      // Row-cyclic Jacobi SVD algorithm with column pivoting

      EMPTSW = ( N*( N-1 ) ) / 2
      NOTROT = 0
      FASTR( 1 ) = ZERO

      // A is represented in factored form A = A * diag(WORK), where diag(WORK)
      // is initialized to identity. WORK is updated during fast scaled
      // rotations.

      DO 1868 q = 1, N
         WORK( q ) = ONE
 1868 CONTINUE


      SWBAND = 3
*[TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
      // if DGESVJ is used as a computational routine in the preconditioned
      // Jacobi SVD algorithm DGESVJ. For sweeps i=1:SWBAND the procedure
      // works on pivots inside a band-like region around the diagonal.
      // The boundaries are determined dynamically, based on the number of
      // pivots above a threshold.

      KBL = MIN( 8, N )
*[TP] KBL is a tuning parameter that defines the tile size in the
     t // iling of the p-q loops of pivot pairs. In general, an optimal
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
     t // uning parameters.
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

            CALL DGSVJ0( JOBV, M-N34, N-N34, A( N34+1, N34+1 ), LDA, WORK( N34+1 ), SVA( N34+1 ), MVL, V( N34*q+1, N34+1 ), LDV, EPSLN, SFMIN, TOL, 2, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ0( JOBV, M-N2, N34-N2, A( N2+1, N2+1 ), LDA, WORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 2, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ1( JOBV, M-N2, N-N2, N4, A( N2+1, N2+1 ), LDA, WORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ0( JOBV, M-N4, N2-N4, A( N4+1, N4+1 ), LDA, WORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ0( JOBV, M, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ1( JOBV, M, N2, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR )


         } else if ( UPPER ) {


            CALL DGSVJ0( JOBV, N4, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 2, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ0( JOBV, N2, N4, A( 1, N4+1 ), LDA, WORK( N4+1 ), SVA( N4+1 ), MVL, V( N4*q+1, N4+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ1( JOBV, N2, N2, N4, A, LDA, WORK, SVA, MVL, V, LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR )

            CALL DGSVJ0( JOBV, N2+N4, N4, A( 1, N2+1 ), LDA, WORK( N2+1 ), SVA( N2+1 ), MVL, V( N2*q+1, N2+1 ), LDV, EPSLN, SFMIN, TOL, 1, WORK( N+1 ), LWORK-N, IERR )

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

                  q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
                  if ( p.NE.q ) {
                     CALL DSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
                     IF( RSVEC )CALL DSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
                     TEMP1 = SVA( p )
                     SVA( p ) = SVA( q )
                     SVA( q ) = TEMP1
                     TEMP1 = WORK( p )
                     WORK( p ) = WORK( q )
                     WORK( q ) = TEMP1
                  }

                  if ( ir1.EQ.0 ) {

         // Column norms are periodically updated by explicit
         // norm computation.
         // Caveat:
         // Unfortunately, some BLAS implementations compute DNRM2(M,A(1,p),1)
         // as DSQRT(DDOT(M,A(1,p),1,A(1,p),1)), which may cause the result to
         // overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and to
         // underflow for ||A(:,p)||_2 < DSQRT(underflow_threshold).
         // Hence, DNRM2 cannot be trusted, not even in the case when
        t // he true norm is far from the under(over)flow boundaries.
         // If properly implemented DNRM2 is available, the IF-THEN-ELSE
         // below should read "AAPP = DNRM2( M, A(1,p), 1 ) * WORK(p)".

                     if ( ( SVA( p ).LT.ROOTBIG ) .AND. ( SVA( p ).GT.ROOTSFMIN ) ) {
                        SVA( p ) = DNRM2( M, A( 1, p ), 1 )*WORK( p )
                     } else {
                        TEMP1 = ZERO
                        AAPP = ONE
                        CALL DLASSQ( M, A( 1, p ), 1, TEMP1, AAPP )
                        SVA( p ) = TEMP1*DSQRT( AAPP )*WORK( p )
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
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 CALL DCOPY( M, A( 1, p ), 1, WORK( N+1 ), 1 )                                  CALL DLASCL( 'G', 0, 0, AAPP, WORK( p ), M, 1, WORK( N+1 ), LDA, IERR )
                                 AAPQ = DDOT( M, WORK( N+1 ), 1, A( 1, q ), 1 )*WORK( q ) / AAQQ
                              }
                           } else {
                              ROTOK = AAPP.LE.( AAQQ / SMALL )
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 CALL DCOPY( M, A( 1, q ), 1, WORK( N+1 ), 1 )                                  CALL DLASCL( 'G', 0, 0, AAQQ, WORK( q ), M, 1, WORK( N+1 ), LDA, IERR )
                                 AAPQ = DDOT( M, WORK( N+1 ), 1, A( 1, p ), 1 )*WORK( p ) / AAPP
                              }
                           }

                           MXAAPQ = MAX( MXAAPQ, DABS( AAPQ ) )

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( DABS( AAPQ ).GT.TOL ) {

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
                                 THETA = -HALF*DABS(AQOAP-APOAQ)/AAPQ

                                 if ( DABS( THETA ).GT.BIGTHETA ) {

                                    T = HALF / THETA
                                    FASTR( 3 ) = T*WORK( p ) / WORK( q )
                                    FASTR( 4 ) = -T*WORK( q ) / WORK( p )                                     CALL DROTM( M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, DABS( T ) )

                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -DSIGN( ONE, AAPQ )
                                    T = ONE / ( THETA+THSIGN* DSQRT( ONE+THETA*THETA ) )
                                    CS = DSQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS

                                    MXSINJ = MAX( MXSINJ, DABS( SN ) )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )

                                    APOAQ = WORK( p ) / WORK( q )
                                    AQOAP = WORK( q ) / WORK( p )
                                    if ( WORK( p ).GE.ONE ) {
                                       if ( WORK( q ).GE.ONE ) {
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q )*CS
                                          CALL DROTM( M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR )
                                       } else {
                                          CALL DAXPY( M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                           CALL DAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q ) / CS
                                          if ( RSVEC ) {
                                             CALL DAXPY( MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                              CALL DAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )
                                          }
                                       }
                                    } else {
                                       if ( WORK( q ).GE.ONE ) {
                                          CALL DAXPY( M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                           CALL DAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )
                                          WORK( p ) = WORK( p ) / CS
                                          WORK( q ) = WORK( q )*CS
                                          if ( RSVEC ) {
                                             CALL DAXPY( MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                              CALL DAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )
                                          }
                                       } else {
                                          IF( WORK( p ).GE.WORK( q ) ) THEN                                              CALL DAXPY( M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                              CALL DAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )
                                             WORK( p ) = WORK( p )*CS
                                             WORK( q ) = WORK( q ) / CS
                                             if ( RSVEC ) {
                                                CALL DAXPY( MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                                 CALL DAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )
                                             }
                                          } else {
                                             CALL DAXPY( M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                              CALL DAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )
                                             WORK( p ) = WORK( p ) / CS
                                             WORK( q ) = WORK( q )*CS
                                             if ( RSVEC ) {
                                                CALL DAXPY( MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                                 CALL DAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )
                                             }
                                          }
                                       }
                                    }
                                 }

                              } else {
               // .. have to use modified Gram-Schmidt like transformation
                                 CALL DCOPY( M, A( 1, p ), 1, WORK( N+1 ), 1 )                                  CALL DLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK( N+1 ), LDA, IERR )
                                 CALL DLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )
                                 TEMP1 = -AAPQ*WORK( p ) / WORK( q )
                                 CALL DAXPY( M, TEMP1, WORK( N+1 ), 1, A( 1, q ), 1 )                                  CALL DLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )                                  SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                 MXSINJ = MAX( MXSINJ, SFMIN )
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q), SVA(p)
            // recompute SVA(q), SVA(p).

                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) THEN                                     SVA( q ) = DNRM2( M, A( 1, q ), 1 )* WORK( q )
                                 } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    CALL DLASSQ( M, A( 1, q ), 1, T, AAQQ )
                                    SVA( q ) = T*DSQRT( AAQQ )*WORK( q )
                                 }
                              }
                              if ( ( AAPP / AAPP0 ).LE.ROOTEPS ) {
                                 IF( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) THEN                                     AAPP = DNRM2( M, A( 1, p ), 1 )* WORK( p )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    CALL DLASSQ( M, A( 1, p ), 1, T, AAPP )
                                    AAPP = T*DSQRT( AAPP )*WORK( p )
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
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 CALL DCOPY( M, A( 1, p ), 1, WORK( N+1 ), 1 )                                  CALL DLASCL( 'G', 0, 0, AAPP, WORK( p ), M, 1, WORK( N+1 ), LDA, IERR )
                                 AAPQ = DDOT( M, WORK( N+1 ), 1, A( 1, q ), 1 )*WORK( q ) / AAQQ
                              }
                           } else {
                              if ( AAPP.GE.AAQQ ) {
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              } else {
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              }
                              if ( AAPP.GT.( SMALL / AAQQ ) ) {
                                 AAPQ = ( DDOT( M, A( 1, p ), 1, A( 1, q ), 1 )*WORK( p )*WORK( q ) / AAQQ ) / AAPP
                              } else {
                                 CALL DCOPY( M, A( 1, q ), 1, WORK( N+1 ), 1 )                                  CALL DLASCL( 'G', 0, 0, AAQQ, WORK( q ), M, 1, WORK( N+1 ), LDA, IERR )
                                 AAPQ = DDOT( M, WORK( N+1 ), 1, A( 1, p ), 1 )*WORK( p ) / AAPP
                              }
                           }

                           MXAAPQ = MAX( MXAAPQ, DABS( AAPQ ) )

         // TO rotate or NOT to rotate, THAT is the question ...

                           if ( DABS( AAPQ ).GT.TOL ) {
                              NOTROT = 0
*[RTD]      ROTATED  = ROTATED + 1
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1

                              if ( ROTOK ) {

                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*DABS(AQOAP-APOAQ)/AAPQ
                                 IF( AAQQ.GT.AAPP0 )THETA = -THETA

                                 if ( DABS( THETA ).GT.BIGTHETA ) {
                                    T = HALF / THETA
                                    FASTR( 3 ) = T*WORK( p ) / WORK( q )
                                    FASTR( 4 ) = -T*WORK( q ) / WORK( p )                                     CALL DROTM( M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                     IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, DABS( T ) )
                                 } else {

                  // .. choose correct signum for THETA and rotate

                                    THSIGN = -DSIGN( ONE, AAPQ )
                                    IF( AAQQ.GT.AAPP0 )THSIGN = -THSIGN
                                    T = ONE / ( THETA+THSIGN* DSQRT( ONE+THETA*THETA ) )
                                    CS = DSQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
                                    MXSINJ = MAX( MXSINJ, DABS( SN ) )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ ) )                                     AAPP = AAPP*DSQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ ) )

                                    APOAQ = WORK( p ) / WORK( q )
                                    AQOAP = WORK( q ) / WORK( p )
                                    if ( WORK( p ).GE.ONE ) {

                                       if ( WORK( q ).GE.ONE ) {
                                          FASTR( 3 ) = T*APOAQ
                                          FASTR( 4 ) = -T*AQOAP
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q )*CS
                                          CALL DROTM( M, A( 1, p ), 1, A( 1, q ), 1, FASTR )                                           IF( RSVEC )CALL DROTM( MVL, V( 1, p ), 1, V( 1, q ), 1, FASTR )
                                       } else {
                                          CALL DAXPY( M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                           CALL DAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )
                                          if ( RSVEC ) {
                                             CALL DAXPY( MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                              CALL DAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )
                                          }
                                          WORK( p ) = WORK( p )*CS
                                          WORK( q ) = WORK( q ) / CS
                                       }
                                    } else {
                                       if ( WORK( q ).GE.ONE ) {
                                          CALL DAXPY( M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                           CALL DAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )
                                          if ( RSVEC ) {
                                             CALL DAXPY( MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                              CALL DAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )
                                          }
                                          WORK( p ) = WORK( p ) / CS
                                          WORK( q ) = WORK( q )*CS
                                       } else {
                                          IF( WORK( p ).GE.WORK( q ) ) THEN                                              CALL DAXPY( M, -T*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )                                              CALL DAXPY( M, CS*SN*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )
                                             WORK( p ) = WORK( p )*CS
                                             WORK( q ) = WORK( q ) / CS
                                             if ( RSVEC ) {
                                                CALL DAXPY( MVL, -T*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )                                                 CALL DAXPY( MVL, CS*SN*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )
                                             }
                                          } else {
                                             CALL DAXPY( M, T*APOAQ, A( 1, p ), 1, A( 1, q ), 1 )                                              CALL DAXPY( M, -CS*SN*AQOAP, A( 1, q ), 1, A( 1, p ), 1 )
                                             WORK( p ) = WORK( p ) / CS
                                             WORK( q ) = WORK( q )*CS
                                             if ( RSVEC ) {
                                                CALL DAXPY( MVL, T*APOAQ, V( 1, p ), 1, V( 1, q ), 1 )                                                 CALL DAXPY( MVL, -CS*SN*AQOAP, V( 1, q ), 1, V( 1, p ), 1 )
                                             }
                                          }
                                       }
                                    }
                                 }

                              } else {
                                 if ( AAPP.GT.AAQQ ) {
                                    CALL DCOPY( M, A( 1, p ), 1, WORK( N+1 ), 1 )                                     CALL DLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK( N+1 ), LDA, IERR )                                     CALL DLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )
                                    TEMP1 = -AAPQ*WORK( p ) / WORK( q )
                                    CALL DAXPY( M, TEMP1, WORK( N+1 ), 1, A( 1, q ), 1 )                                     CALL DLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )
                                    SVA( q ) = AAQQ*DSQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 } else {
                                    CALL DCOPY( M, A( 1, q ), 1, WORK( N+1 ), 1 )                                     CALL DLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK( N+1 ), LDA, IERR )                                     CALL DLASCL( 'G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR )
                                    TEMP1 = -AAPQ*WORK( q ) / WORK( p )
                                    CALL DAXPY( M, TEMP1, WORK( N+1 ), 1, A( 1, p ), 1 )                                     CALL DLASCL( 'G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR )
                                    SVA( p ) = AAPP*DSQRT( MAX( ZERO, ONE-AAPQ*AAPQ ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                                 }
                              }
            // END IF ROTOK THEN ... ELSE

            // In the case of cancellation in updating SVA(q)
            // .. recompute SVA(q)
                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) THEN                                     SVA( q ) = DNRM2( M, A( 1, q ), 1 )* WORK( q )
                                 } else {
                                    T = ZERO
                                    AAQQ = ONE
                                    CALL DLASSQ( M, A( 1, q ), 1, T, AAQQ )
                                    SVA( q ) = T*DSQRT( AAQQ )*WORK( q )
                                 }
                              }
                              if ( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) {
                                 IF( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) THEN                                     AAPP = DNRM2( M, A( 1, p ), 1 )* WORK( p )
                                 } else {
                                    T = ZERO
                                    AAPP = ONE
                                    CALL DLASSQ( M, A( 1, p ), 1, T, AAPP )
                                    AAPP = T*DSQRT( AAPP )*WORK( p )
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
               SVA( p ) = DABS( SVA( p ) )
 2012       CONTINUE
***
 2000    CONTINUE
*2000 :: end of the ibr-loop

      // .. update SVA(N)
         if ( ( SVA( N ).LT.ROOTBIG ) .AND. ( SVA( N ).GT.ROOTSFMIN ) ) {
            SVA( N ) = DNRM2( M, A( 1, N ), 1 )*WORK( N )
         } else {
            T = ZERO
            AAPP = ONE
            CALL DLASSQ( M, A( 1, N ), 1, T, AAPP )
            SVA( N ) = T*DSQRT( AAPP )*WORK( N )
         }

      // Additional steering devices

         IF( ( i.LT.SWBAND ) .AND. ( ( MXAAPQ.LE.ROOTTOL ) .OR. ( ISWROT.LE.N ) ) )SWBAND = i

         if ( ( i.GT.SWBAND+1 ) .AND. ( MXAAPQ.LT.DSQRT( DBLE( N ) )* TOL ) .AND. ( DBLE( N )*MXAAPQ*MXSINJ.LT.TOL ) ) {
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
     t // he underflow threshold.

      N2 = 0
      N4 = 0
      DO 5991 p = 1, N - 1
         q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         if ( p.NE.q ) {
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            TEMP1 = WORK( p )
            WORK( p ) = WORK( q )
            WORK( q ) = TEMP1
            CALL DSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
            IF( RSVEC )CALL DSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
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
         DO 1998 p = 1, N2
            CALL DSCAL( M, WORK( p ) / SVA( p ), A( 1, p ), 1 )
 1998    CONTINUE
      }

      // Scale the product of Jacobi rotations (assemble the fast rotations).

      if ( RSVEC ) {
         if ( APPLV ) {
            DO 2398 p = 1, N
               CALL DSCAL( MVL, WORK( p ), V( 1, p ), 1 )
 2398       CONTINUE
         } else {
            DO 2399 p = 1, N
               TEMP1 = ONE / DNRM2( MVL, V( 1, p ), 1 )
               CALL DSCAL( MVL, TEMP1, V( 1, p ), 1 )
 2399       CONTINUE
         }
      }

      // Undo scaling, if necessary (and possible).
      if ( ( ( SKL.GT.ONE ) .AND. ( SVA( 1 ).LT.( BIG / SKL) ) ) .OR. ( ( SKL.LT.ONE ) .AND. ( SVA( MAX( N2, 1 ) ) .GT. ( SFMIN / SKL) ) ) ) {
         DO 2400 p = 1, N
            SVA( P ) = SKL*SVA( P )
 2400    CONTINUE
         SKL= ONE
      }

      WORK( 1 ) = SKL
      // The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE
     t // hen some of the singular values may overflow or underflow and
     t // he spectrum is given in this factored representation.

      WORK( 2 ) = DBLE( N4 )
      // N4 is the number of computed nonzero singular values of A.

      WORK( 3 ) = DBLE( N2 )
      // N2 is the number of singular values of A greater than SFMIN.
      // If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
     t // hat may carry some information.

      WORK( 4 ) = DBLE( i )
      // i is the index of the last sweep before declaring convergence.

      WORK( 5 ) = MXAAPQ
      // MXAAPQ is the largest absolute value of scaled pivots in the
      // last sweep

      WORK( 6 ) = MXSINJ
      // MXSINJ is the largest absolute value of the sines of Jacobi angles
      // in the last sweep

      RETURN
      // ..
      // .. END OF DGESVJ
      // ..
      }
