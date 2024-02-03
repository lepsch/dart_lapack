      SUBROUTINE DGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, M, N, A, LDA, SVA, U, LDU, V, LDV, WORK, LWORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      IMPLICIT    NONE
      int         INFO, LDA, LDU, LDV, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double           A( LDA, * ), SVA( N ), U( LDU, * ), V( LDV, * ), WORK( LWORK );
      int         IWORK( * );
      String      JOBA, JOBP, JOBR, JOBT, JOBU, JOBV;
      // ..

*  ===========================================================================

      // .. Local Parameters ..
      double             ZERO,  ONE;
      const     ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      double           AAPP, AAQQ, AATMAX, AATMIN, BIG, BIG1, COND_OK, CONDR1, CONDR2, ENTRA,  ENTRAT, EPSLN,  MAXPRJ, SCALEM, SCONDA, SFMIN,  SMALL,  TEMP1,  USCAL1, USCAL2, XSC;
      int     IERR,   N1,     NR,     NUMRANK,        p, q,   WARNING;
      bool    ALMORT, DEFR,   ERREST, GOSCAL, JRACC,  KILL,   LSVEC, L2ABER, L2KILL, L2PERT, L2RANK, L2TRAN, NOSCAL, ROWPIV, RSVEC,  TRANSP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DABS, DLOG, MAX, MIN, DBLE, IDNINT, DSIGN, DSQRT
      // ..
      // .. External Functions ..
      double            DLAMCH, DNRM2;
      int       IDAMAX;
      bool      LSAME;
      // EXTERNAL IDAMAX, LSAME, DLAMCH, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY,  DGELQF, DGEQP3, DGEQRF, DLACPY, DLASCL, DLASET, DLASSQ, DLASWP, DORGQR, DORMLQ, DORMQR, DPOCON, DSCAL,  DSWAP,  DTRSM,  XERBLA

      // EXTERNAL DGESVJ
      // ..

      // Test the input arguments

      LSVEC  = LSAME( JOBU, 'U' ) .OR. LSAME( JOBU, 'F' )
      JRACC  = LSAME( JOBV, 'J' )
      RSVEC  = LSAME( JOBV, 'V' ) .OR. JRACC
      ROWPIV = LSAME( JOBA, 'F' ) .OR. LSAME( JOBA, 'G' )
      L2RANK = LSAME( JOBA, 'R' )
      L2ABER = LSAME( JOBA, 'A' )
      ERREST = LSAME( JOBA, 'E' ) .OR. LSAME( JOBA, 'G' )
      L2TRAN = LSAME( JOBT, 'T' )
      L2KILL = LSAME( JOBR, 'R' )
      DEFR   = LSAME( JOBR, 'N' )
      L2PERT = LSAME( JOBP, 'P' )

      if ( .NOT.(ROWPIV .OR. L2RANK .OR. L2ABER .OR. ERREST .OR. LSAME( JOBA, 'C' ) )) {
         INFO = - 1
      } else if ( .NOT.( LSVEC  .OR. LSAME( JOBU, 'N' ) .OR. LSAME( JOBU, 'W' )) ) {
         INFO = - 2
      } else if ( .NOT.( RSVEC .OR. LSAME( JOBV, 'N' ) .OR. LSAME( JOBV, 'W' )) .OR. ( JRACC .AND. (.NOT.LSVEC) ) ) {
         INFO = - 3
      } else if ( .NOT. ( L2KILL .OR. DEFR ) ) {
         INFO = - 4
      } else if ( .NOT. ( L2TRAN .OR. LSAME( JOBT, 'N' ) ) ) {
         INFO = - 5
      } else if ( .NOT. ( L2PERT .OR. LSAME( JOBP, 'N' ) ) ) {
         INFO = - 6
      } else if ( M .LT. 0 ) {
         INFO = - 7
      } else if ( ( N .LT. 0 ) .OR. ( N .GT. M ) ) {
         INFO = - 8
      } else if ( LDA .LT. M ) {
         INFO = - 10
      } else if ( LSVEC .AND. ( LDU .LT. M ) ) {
         INFO = - 13
      } else if ( RSVEC .AND. ( LDV .LT. N ) ) {
         INFO = - 15
      } else if ( (.NOT.(LSVEC .OR. RSVEC .OR. ERREST).AND. (LWORK .LT. MAX(7,4*N+1,2*M+N))) .OR. (.NOT.(LSVEC .OR. RSVEC) .AND. ERREST .AND. (LWORK .LT. MAX(7,4*N+N*N,2*M+N))) .OR. (LSVEC .AND. (.NOT.RSVEC) .AND. (LWORK .LT. MAX(7,2*M+N,4*N+1))) .OR. (RSVEC .AND. (.NOT.LSVEC) .AND. (LWORK .LT. MAX(7,2*M+N,4*N+1))) .OR. (LSVEC .AND. RSVEC .AND. (.NOT.JRACC) .AND. (LWORK.LT.MAX(2*M+N,6*N+2*N*N))) .OR. (LSVEC .AND. RSVEC .AND. JRACC .AND. LWORK.LT.MAX(2*M+N,4*N+N*N,2*N+N*N+6))) {
         INFO = - 17
      } else {
         // #:)
         INFO = 0
      }

      if ( INFO .NE. 0 ) {
        // #:(
         CALL XERBLA( 'DGEJSV', - INFO )
         RETURN
      }

      // Quick return for void matrix (Y3K safe)
* #:)
      if ( ( M .EQ. 0 ) .OR. ( N .EQ. 0 ) ) {
         IWORK(1:3) = 0
         WORK(1:7) = 0
         RETURN
      ENDIF

      // Determine whether the matrix U should be M x N or M x M

      if ( LSVEC ) {
         N1 = N
         IF ( LSAME( JOBU, 'F' ) ) N1 = M
      }

      // Set numerical parameters

*!    NOTE: Make sure DLAMCH() does not fail on the target architecture.

      EPSLN = DLAMCH('Epsilon')
      SFMIN = DLAMCH('SafeMinimum')
      SMALL = SFMIN / EPSLN
      BIG   = DLAMCH('O')
      // BIG   = ONE / SFMIN

      // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N

*(!)  If necessary, scale SVA() to protect the largest norm from
      // overflow. It is possible that this scaling pushes the smallest
      // column norm left from the underflow threshold (extreme case).

      SCALEM  = ONE / DSQRT(DBLE(M)*DBLE(N))
      NOSCAL  = .TRUE.
      GOSCAL  = .TRUE.
      DO 1874 p = 1, N
         AAPP = ZERO
         AAQQ = ONE
         CALL DLASSQ( M, A(1,p), 1, AAPP, AAQQ )
         if ( AAPP .GT. BIG ) {
            INFO = - 9
            CALL XERBLA( 'DGEJSV', -INFO )
            RETURN
         }
         AAQQ = DSQRT(AAQQ)
         if ( ( AAPP .LT. (BIG / AAQQ) ) .AND. NOSCAL  ) {
            SVA(p)  = AAPP * AAQQ
         } else {
            NOSCAL  = .FALSE.
            SVA(p)  = AAPP * ( AAQQ * SCALEM )
            if ( GOSCAL ) {
               GOSCAL = .FALSE.
               CALL DSCAL( p-1, SCALEM, SVA, 1 )
            }
         }
 1874 CONTINUE

      IF ( NOSCAL ) SCALEM = ONE

      AAPP = ZERO
      AAQQ = BIG
      DO 4781 p = 1, N
         AAPP = MAX( AAPP, SVA(p) )
         IF ( SVA(p) .NE. ZERO ) AAQQ = MIN( AAQQ, SVA(p) )
 4781 CONTINUE

      // Quick return for zero M x N matrix
* #:)
      if ( AAPP .EQ. ZERO ) {
         IF ( LSVEC ) CALL DLASET( 'G', M, N1, ZERO, ONE, U, LDU )
         IF ( RSVEC ) CALL DLASET( 'G', N, N,  ZERO, ONE, V, LDV )
         WORK(1) = ONE
         WORK(2) = ONE
         IF ( ERREST ) WORK(3) = ONE
         if ( LSVEC .AND. RSVEC ) {
            WORK(4) = ONE
            WORK(5) = ONE
         }
         if ( L2TRAN ) {
            WORK(6) = ZERO
            WORK(7) = ZERO
         }
         IWORK(1) = 0
         IWORK(2) = 0
         IWORK(3) = 0
         RETURN
      }

      // Issue warning if denormalized column norms detected. Override the
      // high relative accuracy request. Issue licence to kill columns
      // (set them to zero) whose norm is less than sigma_max / BIG (roughly).
* #:(
      WARNING = 0
      if ( AAQQ .LE. SFMIN ) {
         L2RANK = .TRUE.
         L2KILL = .TRUE.
         WARNING = 1
      }

      // Quick return for one-column matrix
* #:)
      if ( N .EQ. 1 ) {

         if ( LSVEC ) {
            CALL DLASCL( 'G',0,0,SVA(1),SCALEM, M,1,A(1,1),LDA,IERR )
            CALL DLACPY( 'A', M, 1, A, LDA, U, LDU )
            // computing all M left singular vectors of the M x 1 matrix
            if ( N1 .NE. N  ) {
               CALL DGEQRF( M, N, U,LDU, WORK, WORK(N+1),LWORK-N,IERR )
               CALL DORGQR( M,N1,1, U,LDU,WORK,WORK(N+1),LWORK-N,IERR )
               CALL DCOPY( M, A(1,1), 1, U(1,1), 1 )
            }
         }
         if ( RSVEC ) {
             V(1,1) = ONE
         }
         if ( SVA(1) .LT. (BIG*SCALEM) ) {
            SVA(1)  = SVA(1) / SCALEM
            SCALEM  = ONE
         }
         WORK(1) = ONE / SCALEM
         WORK(2) = ONE
         if ( SVA(1) .NE. ZERO ) {
            IWORK(1) = 1
            if ( ( SVA(1) / SCALEM) .GE. SFMIN ) {
               IWORK(2) = 1
            } else {
               IWORK(2) = 0
            }
         } else {
            IWORK(1) = 0
            IWORK(2) = 0
         }
         IWORK(3) = 0
         IF ( ERREST ) WORK(3) = ONE
         if ( LSVEC .AND. RSVEC ) {
            WORK(4) = ONE
            WORK(5) = ONE
         }
         if ( L2TRAN ) {
            WORK(6) = ZERO
            WORK(7) = ZERO
         }
         RETURN

      }

      TRANSP = .FALSE.
      L2TRAN = L2TRAN .AND. ( M .EQ. N )

      AATMAX = -ONE
      AATMIN =  BIG
      if ( ROWPIV .OR. L2TRAN ) {

      // Compute the row norms, needed to determine row pivoting sequence
      // (in the case of heavily row weighted A, row pivoting is strongly
      // advised) and to collect information needed to compare the
      // structures of A * A^t and A^t * A (in the case L2TRAN.EQ..TRUE.).

         if ( L2TRAN ) {
            DO 1950 p = 1, M
               XSC   = ZERO
               TEMP1 = ONE
               CALL DLASSQ( N, A(p,1), LDA, XSC, TEMP1 )
               // DLASSQ gets both the ell_2 and the ell_infinity norm
               // in one pass through the vector
               WORK(M+N+p)  = XSC * SCALEM
               WORK(N+p)    = XSC * (SCALEM*DSQRT(TEMP1))
               AATMAX = MAX( AATMAX, WORK(N+p) )
               IF (WORK(N+p) .NE. ZERO) AATMIN = MIN(AATMIN,WORK(N+p))
 1950       CONTINUE
         } else {
            DO 1904 p = 1, M
               WORK(M+N+p) = SCALEM*DABS( A(p,IDAMAX(N,A(p,1),LDA)) )
               AATMAX = MAX( AATMAX, WORK(M+N+p) )
               AATMIN = MIN( AATMIN, WORK(M+N+p) )
 1904       CONTINUE
         }

      }

      // For square matrix A try to determine whether A^t  would be  better
      // input for the preconditioned Jacobi SVD, with faster convergence.
      // The decision is based on an O(N) function of the vector of column
      // and row norms of A, based on the Shannon entropy. This should give
     t // he right choice in most cases when the difference actually matters.
      // It may fail and pick the slower converging side.

      ENTRA  = ZERO
      ENTRAT = ZERO
      if ( L2TRAN ) {

         XSC   = ZERO
         TEMP1 = ONE
         CALL DLASSQ( N, SVA, 1, XSC, TEMP1 )
         TEMP1 = ONE / TEMP1

         ENTRA = ZERO
         DO 1113 p = 1, N
            BIG1  = ( ( SVA(p) / XSC )**2 ) * TEMP1
            IF ( BIG1 .NE. ZERO ) ENTRA = ENTRA + BIG1 * DLOG(BIG1)
 1113    CONTINUE
         ENTRA = - ENTRA / DLOG(DBLE(N))

         // Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex.
         // It is derived from the diagonal of  A^t * A.  Do the same with the
         // diagonal of A * A^t, compute the entropy of the corresponding
         // probability distribution. Note that A * A^t and A^t * A have the
         // same trace.

         ENTRAT = ZERO
         DO 1114 p = N+1, N+M
            BIG1 = ( ( WORK(p) / XSC )**2 ) * TEMP1
            IF ( BIG1 .NE. ZERO ) ENTRAT = ENTRAT + BIG1 * DLOG(BIG1)
 1114    CONTINUE
         ENTRAT = - ENTRAT / DLOG(DBLE(M))

         // Analyze the entropies and decide A or A^t. Smaller entropy
         // usually means better input for the algorithm.

         TRANSP = ( ENTRAT .LT. ENTRA )

         // If A^t is better than A, transpose A.

         if ( TRANSP ) {
            // In an optimal implementation, this trivial transpose
            // should be replaced with faster transpose.
            DO 1115 p = 1, N - 1
               DO 1116 q = p + 1, N
                   TEMP1 = A(q,p)
                  A(q,p) = A(p,q)
                  A(p,q) = TEMP1
 1116          CONTINUE
 1115       CONTINUE
            DO 1117 p = 1, N
               WORK(M+N+p) = SVA(p)
               SVA(p)      = WORK(N+p)
 1117       CONTINUE
            TEMP1  = AAPP
            AAPP   = AATMAX
            AATMAX = TEMP1
            TEMP1  = AAQQ
            AAQQ   = AATMIN
            AATMIN = TEMP1
            KILL   = LSVEC
            LSVEC  = RSVEC
            RSVEC  = KILL
            IF ( LSVEC ) N1 = N

            ROWPIV = .TRUE.
         }

      }
      // END IF L2TRAN

      // Scale the matrix so that its maximal singular value remains less
     t // han DSQRT(BIG) -- the matrix is scaled so that its maximal column
      // has Euclidean norm equal to DSQRT(BIG/N). The only reason to keep
      // DSQRT(BIG) instead of BIG is the fact that DGEJSV uses LAPACK and
      // BLAS routines that, in some implementations, are not capable of
      // working in the full interval [SFMIN,BIG] and that they may provoke
      // overflows in the intermediate results. If the singular values spread
      // from SFMIN to BIG, then DGESVJ will compute them. So, in that case,
      // one should use DGESVJ instead of DGEJSV.

      BIG1   = DSQRT( BIG )
      TEMP1  = DSQRT( BIG / DBLE(N) )

      CALL DLASCL( 'G', 0, 0, AAPP, TEMP1, N, 1, SVA, N, IERR )
      if ( AAQQ .GT. (AAPP * SFMIN) ) {
          AAQQ = ( AAQQ / AAPP ) * TEMP1
      } else {
          AAQQ = ( AAQQ * TEMP1 ) / AAPP
      }
      TEMP1 = TEMP1 * SCALEM
      CALL DLASCL( 'G', 0, 0, AAPP, TEMP1, M, N, A, LDA, IERR )

      // To undo scaling at the end of this procedure, multiply the
      // computed singular values with USCAL2 / USCAL1.

      USCAL1 = TEMP1
      USCAL2 = AAPP

      if ( L2KILL ) {
         // L2KILL enforces computation of nonzero singular values in
        t // he restricted range of condition number of the initial A,
         // sigma_max(A) / sigma_min(A) approx. DSQRT(BIG)/DSQRT(SFMIN).
         XSC = DSQRT( SFMIN )
      } else {
         XSC = SMALL

         // Now, if the condition number of A is too big,
         // sigma_max(A) / sigma_min(A) .GT. DSQRT(BIG/N) * EPSLN / SFMIN,
         // as a precaution measure, the full SVD is computed using DGESVJ
         // with accumulated Jacobi rotations. This provides numerically
         // more robust computation, at the cost of slightly increased run
        t // ime. Depending on the concrete implementation of BLAS and LAPACK
         // (i.e. how they behave in presence of extreme ill-conditioning) the
         // implementor may decide to remove this switch.
         if ( ( AAQQ.LT.DSQRT(SFMIN) ) .AND. LSVEC .AND. RSVEC ) {
            JRACC = .TRUE.
         }

      }
      if ( AAQQ .LT. XSC ) {
         DO 700 p = 1, N
            if ( SVA(p) .LT. XSC ) {
               CALL DLASET( 'A', M, 1, ZERO, ZERO, A(1,p), LDA )
               SVA(p) = ZERO
            }
 700     CONTINUE
      }

      // Preconditioning using QR factorization with pivoting

      if ( ROWPIV ) {
         // Optional row permutation (Bjoerck row pivoting):
         // A result by Cox and Higham shows that the Bjoerck's
         // row pivoting combined with standard column pivoting
         // has similar effect as Powell-Reid complete pivoting.
         // The ell-infinity norms of A are made nonincreasing.
         DO 1952 p = 1, M - 1
            q = IDAMAX( M-p+1, WORK(M+N+p), 1 ) + p - 1
            IWORK(2*N+p) = q
            if ( p .NE. q ) {
               TEMP1       = WORK(M+N+p)
               WORK(M+N+p) = WORK(M+N+q)
               WORK(M+N+q) = TEMP1
            }
 1952    CONTINUE
         CALL DLASWP( N, A, LDA, 1, M-1, IWORK(2*N+1), 1 )
      }

      // End of the preparation phase (scaling, optional sorting and
     t // ransposing, optional flushing of small columns).

      // Preconditioning

      // If the full SVD is needed, the right singular vectors are computed
      // from a matrix equation, and for that we need theoretical analysis
      // of the Businger-Golub pivoting. So we use DGEQP3 as the first RR QRF.
      // In all other cases the first RR QRF can be chosen by other criteria
      // (eg speed by replacing global with restricted window pivoting, such
      // as in SGEQPX from TOMS # 782). Good results will be obtained using
      // SGEQPX with properly (!) chosen numerical parameters.
      // Any improvement of DGEQP3 improves overall performance of DGEJSV.

      // A * P1 = Q1 * [ R1^t 0]^t:
      DO 1963 p = 1, N
         // .. all columns are free columns
         IWORK(p) = 0
 1963 CONTINUE
      CALL DGEQP3( M,N,A,LDA, IWORK,WORK, WORK(N+1),LWORK-N, IERR )

      // The upper triangular matrix R1 from the first QRF is inspected for
      // rank deficiency and possibilities for deflation, or possible
      // ill-conditioning. Depending on the user specified flag L2RANK,
     t // he procedure explores possibilities to reduce the numerical
      // rank by inspecting the computed upper triangular factor. If
      // L2RANK or L2ABER are up, then DGEJSV will compute the SVD of
      // A + dA, where ||dA|| <= f(M,N)*EPSLN.

      NR = 1
      if ( L2ABER ) {
         // Standard absolute error bound suffices. All sigma_i with
         // sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
         // aggressive enforcement of lower numerical rank by introducing a
         // backward error of the order of N*EPSLN*||A||.
         TEMP1 = DSQRT(DBLE(N))*EPSLN
         DO 3001 p = 2, N
            if ( DABS(A(p,p)) .GE. (TEMP1*DABS(A(1,1))) ) {
               NR = NR + 1
            } else {
               GO TO 3002
            }
 3001    CONTINUE
 3002    CONTINUE
      } else if ( L2RANK ) {
         // .. similarly as above, only slightly more gentle (less aggressive).
         // Sudden drop on the diagonal of R1 is used as the criterion for
         // close-to-rank-deficient.
         TEMP1 = DSQRT(SFMIN)
         DO 3401 p = 2, N
            IF ( ( DABS(A(p,p)) .LT. (EPSLN*DABS(A(p-1,p-1))) ) .OR. ( DABS(A(p,p)) .LT. SMALL ) .OR. ( L2KILL .AND. (DABS(A(p,p)) .LT. TEMP1) ) ) GO TO 3402
            NR = NR + 1
 3401    CONTINUE
 3402    CONTINUE

      } else {
         // The goal is high relative accuracy. However, if the matrix
         // has high scaled condition number the relative accuracy is in
         // general not feasible. Later on, a condition number estimator
         // will be deployed to estimate the scaled condition number.
         // Here we just remove the underflowed part of the triangular
         // factor. This prevents the situation in which the code is
         // working hard to get the accuracy not warranted by the data.
         TEMP1  = DSQRT(SFMIN)
         DO 3301 p = 2, N
            IF ( ( DABS(A(p,p)) .LT. SMALL ) .OR. ( L2KILL .AND. (DABS(A(p,p)) .LT. TEMP1) ) ) GO TO 3302
            NR = NR + 1
 3301    CONTINUE
 3302    CONTINUE

      }

      ALMORT = .FALSE.
      if ( NR .EQ. N ) {
         MAXPRJ = ONE
         DO 3051 p = 2, N
            TEMP1  = DABS(A(p,p)) / SVA(IWORK(p))
            MAXPRJ = MIN( MAXPRJ, TEMP1 )
 3051    CONTINUE
         IF ( MAXPRJ**2 .GE. ONE - DBLE(N)*EPSLN ) ALMORT = .TRUE.
      }


      SCONDA = - ONE
      CONDR1 = - ONE
      CONDR2 = - ONE

      if ( ERREST ) {
         if ( N .EQ. NR ) {
            if ( RSVEC ) {
               // .. V is available as workspace
               CALL DLACPY( 'U', N, N, A, LDA, V, LDV )
               DO 3053 p = 1, N
                  TEMP1 = SVA(IWORK(p))
                  CALL DSCAL( p, ONE/TEMP1, V(1,p), 1 )
 3053          CONTINUE
               CALL DPOCON( 'U', N, V, LDV, ONE, TEMP1, WORK(N+1), IWORK(2*N+M+1), IERR )
            } else if ( LSVEC ) {
               // .. U is available as workspace
               CALL DLACPY( 'U', N, N, A, LDA, U, LDU )
               DO 3054 p = 1, N
                  TEMP1 = SVA(IWORK(p))
                  CALL DSCAL( p, ONE/TEMP1, U(1,p), 1 )
 3054          CONTINUE
               CALL DPOCON( 'U', N, U, LDU, ONE, TEMP1, WORK(N+1), IWORK(2*N+M+1), IERR )
            } else {
               CALL DLACPY( 'U', N, N, A, LDA, WORK(N+1), N )
               DO 3052 p = 1, N
                  TEMP1 = SVA(IWORK(p))
                  CALL DSCAL( p, ONE/TEMP1, WORK(N+(p-1)*N+1), 1 )
 3052          CONTINUE
            // .. the columns of R are scaled to have unit Euclidean lengths.
               CALL DPOCON( 'U', N, WORK(N+1), N, ONE, TEMP1, WORK(N+N*N+1), IWORK(2*N+M+1), IERR )
            }
            SCONDA = ONE / DSQRT(TEMP1)
            // SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1).
            // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
         } else {
            SCONDA = - ONE
         }
      }

      L2PERT = L2PERT .AND. ( DABS( A(1,1)/A(NR,NR) ) .GT. DSQRT(BIG1) )
      // If there is no violent scaling, artificial perturbation is not needed.

      // Phase 3:

      if ( .NOT. ( RSVEC .OR. LSVEC ) ) {

          // Singular Values only

          // .. transpose A(1:NR,1:N)
         DO 1946 p = 1, MIN( N-1, NR )
            CALL DCOPY( N-p, A(p,p+1), LDA, A(p+1,p), 1 )
 1946    CONTINUE

         // The following two DO-loops introduce small relative perturbation
         // into the strict upper triangle of the lower triangular matrix.
         // Small entries below the main diagonal are also changed.
         // This modification is useful if the computing environment does not
         // provide/allow FLUSH TO ZERO underflow, for it prevents many
         // annoying denormalized numbers in case of strongly scaled matrices.
         // The perturbation is structured so that it does not introduce any
         // new perturbation of the singular values, and it does not destroy
        t // he job done by the preconditioner.
         // The licence for this perturbation is in the variable L2PERT, which
         // should be .FALSE. if FLUSH TO ZERO underflow is active.

         if ( .NOT. ALMORT ) {

            if ( L2PERT ) {
               // XSC = DSQRT(SMALL)
               XSC = EPSLN / DBLE(N)
               DO 4947 q = 1, NR
                  TEMP1 = XSC*DABS(A(q,q))
                  DO 4949 p = 1, N
                     IF ( ( (p.GT.q) .AND. (DABS(A(p,q)).LE.TEMP1) ) .OR. ( p .LT. q ) ) A(p,q) = DSIGN( TEMP1, A(p,q) )
 4949             CONTINUE
 4947          CONTINUE
            } else {
               CALL DLASET( 'U', NR-1,NR-1, ZERO,ZERO, A(1,2),LDA )
            }

             // .. second preconditioning using the QR factorization

            CALL DGEQRF( N,NR, A,LDA, WORK, WORK(N+1),LWORK-N, IERR )

            // .. and transpose upper to lower triangular
            DO 1948 p = 1, NR - 1
               CALL DCOPY( NR-p, A(p,p+1), LDA, A(p+1,p), 1 )
 1948       CONTINUE

         }

            // Row-cyclic Jacobi SVD algorithm with column pivoting

            // .. again some perturbation (a "background noise") is added
           t // o drown denormals
            if ( L2PERT ) {
               // XSC = DSQRT(SMALL)
               XSC = EPSLN / DBLE(N)
               DO 1947 q = 1, NR
                  TEMP1 = XSC*DABS(A(q,q))
                  DO 1949 p = 1, NR
                     IF ( ( (p.GT.q) .AND. (DABS(A(p,q)).LE.TEMP1) ) .OR. ( p .LT. q ) ) A(p,q) = DSIGN( TEMP1, A(p,q) )
 1949             CONTINUE
 1947          CONTINUE
            } else {
               CALL DLASET( 'U', NR-1, NR-1, ZERO, ZERO, A(1,2), LDA )
            }

            // .. and one-sided Jacobi rotations are started on a lower
           t // riangular matrix (plus perturbation which is ignored in
           t // he part which destroys triangular form (confusing?!))

            CALL DGESVJ( 'L', 'NoU', 'NoV', NR, NR, A, LDA, SVA, N, V, LDV, WORK, LWORK, INFO )

            SCALEM  = WORK(1)
            NUMRANK = IDNINT(WORK(2))


      } else if ( RSVEC .AND. ( .NOT. LSVEC ) ) {

         // -> Singular Values and Right Singular Vectors <-

         if ( ALMORT ) {

            // .. in this case NR equals N
            DO 1998 p = 1, NR
               CALL DCOPY( N-p+1, A(p,p), LDA, V(p,p), 1 )
 1998       CONTINUE
            CALL DLASET( 'Upper', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV )

            CALL DGESVJ( 'L','U','N', N, NR, V,LDV, SVA, NR, A,LDA, WORK, LWORK, INFO )
            SCALEM  = WORK(1)
            NUMRANK = IDNINT(WORK(2))

         } else {

         // .. two more QR factorizations ( one QRF is not enough, two require
         // accumulated product of Jacobi rotations, three are perfect )

            CALL DLASET( 'Lower', NR-1, NR-1, ZERO, ZERO, A(2,1), LDA )
            CALL DGELQF( NR, N, A, LDA, WORK, WORK(N+1), LWORK-N, IERR)
            CALL DLACPY( 'Lower', NR, NR, A, LDA, V, LDV )
            CALL DLASET( 'Upper', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV )
            CALL DGEQRF( NR, NR, V, LDV, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR )
            DO 8998 p = 1, NR
               CALL DCOPY( NR-p+1, V(p,p), LDV, V(p,p), 1 )
 8998       CONTINUE
            CALL DLASET( 'Upper', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV )

            CALL DGESVJ( 'Lower', 'U','N', NR, NR, V,LDV, SVA, NR, U, LDU, WORK(N+1), LWORK, INFO )
            SCALEM  = WORK(N+1)
            NUMRANK = IDNINT(WORK(N+2))
            if ( NR .LT. N ) {
               CALL DLASET( 'A',N-NR, NR, ZERO,ZERO, V(NR+1,1),   LDV )
               CALL DLASET( 'A',NR, N-NR, ZERO,ZERO, V(1,NR+1),   LDV )
               CALL DLASET( 'A',N-NR,N-NR,ZERO,ONE, V(NR+1,NR+1), LDV )
            }

         CALL DORMLQ( 'Left', 'Transpose', N, N, NR, A, LDA, WORK, V, LDV, WORK(N+1), LWORK-N, IERR )

         }

         DO 8991 p = 1, N
            CALL DCOPY( N, V(p,1), LDV, A(IWORK(p),1), LDA )
 8991    CONTINUE
         CALL DLACPY( 'All', N, N, A, LDA, V, LDV )

         if ( TRANSP ) {
            CALL DLACPY( 'All', N, N, V, LDV, U, LDU )
         }

      } else if ( LSVEC .AND. ( .NOT. RSVEC ) ) {

         // .. Singular Values and Left Singular Vectors                 ..

         // .. second preconditioning step to avoid need to accumulate
         // Jacobi rotations in the Jacobi iterations.
         DO 1965 p = 1, NR
            CALL DCOPY( N-p+1, A(p,p), LDA, U(p,p), 1 )
 1965    CONTINUE
         CALL DLASET( 'Upper', NR-1, NR-1, ZERO, ZERO, U(1,2), LDU )

         CALL DGEQRF( N, NR, U, LDU, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR )

         DO 1967 p = 1, NR - 1
            CALL DCOPY( NR-p, U(p,p+1), LDU, U(p+1,p), 1 )
 1967    CONTINUE
         CALL DLASET( 'Upper', NR-1, NR-1, ZERO, ZERO, U(1,2), LDU )

         CALL DGESVJ( 'Lower', 'U', 'N', NR,NR, U, LDU, SVA, NR, A, LDA, WORK(N+1), LWORK-N, INFO )
         SCALEM  = WORK(N+1)
         NUMRANK = IDNINT(WORK(N+2))

         if ( NR .LT. M ) {
            CALL DLASET( 'A',  M-NR, NR,ZERO, ZERO, U(NR+1,1), LDU )
            if ( NR .LT. N1 ) {
               CALL DLASET( 'A',NR, N1-NR, ZERO, ZERO, U(1,NR+1), LDU )
               CALL DLASET( 'A',M-NR,N1-NR,ZERO,ONE,U(NR+1,NR+1), LDU )
            }
         }

         CALL DORMQR( 'Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR )

         IF ( ROWPIV ) CALL DLASWP( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 )

         DO 1974 p = 1, N1
            XSC = ONE / DNRM2( M, U(1,p), 1 )
            CALL DSCAL( M, XSC, U(1,p), 1 )
 1974    CONTINUE

         if ( TRANSP ) {
            CALL DLACPY( 'All', N, N, U, LDU, V, LDV )
         }

      } else {

         // .. Full SVD ..

         if ( .NOT. JRACC ) {

         if ( .NOT. ALMORT ) {

            // Second Preconditioning Step (QRF [with pivoting])
            // Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
            // equivalent to an LQF CALL. Since in many libraries the QRF
            // seems to be better optimized than the LQF, we do explicit
           t // ranspose and use the QRF. This is subject to changes in an
            // optimized implementation of DGEJSV.

            DO 1968 p = 1, NR
               CALL DCOPY( N-p+1, A(p,p), LDA, V(p,p), 1 )
 1968       CONTINUE

            // .. the following two loops perturb small entries to avoid
            // denormals in the second QR factorization, where they are
            // as good as zeros. This is done to avoid painfully slow
            // computation with denormals. The relative size of the perturbation
            // is a parameter that can be changed by the implementer.
            // This perturbation device will be obsolete on machines with
            // properly implemented arithmetic.
            // To switch it off, set L2PERT=.FALSE. To remove it from  the
            // code, remove the action under L2PERT=.TRUE., leave the ELSE part.
            // The following two loops should be blocked and fused with the
           t // ransposed copy above.

            if ( L2PERT ) {
               XSC = DSQRT(SMALL)
               DO 2969 q = 1, NR
                  TEMP1 = XSC*DABS( V(q,q) )
                  DO 2968 p = 1, N
                     IF ( ( p .GT. q ) .AND. ( DABS(V(p,q)) .LE. TEMP1 ) .OR. ( p .LT. q ) ) V(p,q) = DSIGN( TEMP1, V(p,q) )
                     IF ( p .LT. q ) V(p,q) = - V(p,q)
 2968             CONTINUE
 2969          CONTINUE
            } else {
               CALL DLASET( 'U', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV )
            }

            // Estimate the row scaled condition number of R1
            // (If R1 is rectangular, N > NR, then the condition number
            // of the leading NR x NR submatrix is estimated.)

            CALL DLACPY( 'L', NR, NR, V, LDV, WORK(2*N+1), NR )
            DO 3950 p = 1, NR
               TEMP1 = DNRM2(NR-p+1,WORK(2*N+(p-1)*NR+p),1)
               CALL DSCAL(NR-p+1,ONE/TEMP1,WORK(2*N+(p-1)*NR+p),1)
 3950       CONTINUE
            CALL DPOCON('Lower',NR,WORK(2*N+1),NR,ONE,TEMP1, WORK(2*N+NR*NR+1),IWORK(M+2*N+1),IERR)
            CONDR1 = ONE / DSQRT(TEMP1)
            // .. here need a second opinion on the condition number
            // .. then assume worst case scenario
            // R1 is OK for inverse <=> CONDR1 .LT. DBLE(N)
            // more conservative    <=> CONDR1 .LT. DSQRT(DBLE(N))

            COND_OK = DSQRT(DBLE(NR))
*[TP]       COND_OK is a tuning parameter.

            if ( CONDR1 .LT. COND_OK ) {
               // .. the second QRF without pivoting. Note: in an optimized
               // implementation, this QRF should be implemented as the QRF
               // of a lower triangular matrix.
               // R1^t = Q2 * R2
               CALL DGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR )

               if ( L2PERT ) {
                  XSC = DSQRT(SMALL)/EPSLN
                  DO 3959 p = 2, NR
                     DO 3958 q = 1, p - 1
                        TEMP1 = XSC * MIN(DABS(V(p,p)),DABS(V(q,q)))
                        IF ( DABS(V(q,p)) .LE. TEMP1 ) V(q,p) = DSIGN( TEMP1, V(q,p) )
 3958                CONTINUE
 3959             CONTINUE
               }

               IF ( NR .NE. N ) CALL DLACPY( 'A', N, NR, V, LDV, WORK(2*N+1), N )
               // .. save ...

            // .. this transposed copy should be better than naive
               DO 1969 p = 1, NR - 1
                  CALL DCOPY( NR-p, V(p,p+1), LDV, V(p+1,p), 1 )
 1969          CONTINUE

               CONDR2 = CONDR1

            } else {

               // .. ill-conditioned case: second QRF with pivoting
               // Note that windowed pivoting would be equally good
               // numerically, and more run-time efficient. So, in
               // an optimal implementation, the next call to DGEQP3
               // should be replaced with eg. CALL SGEQPX (ACM TOMS #782)
               // with properly (carefully) chosen parameters.

               // R1^t * P2 = Q2 * R2
               DO 3003 p = 1, NR
                  IWORK(N+p) = 0
 3003          CONTINUE
               CALL DGEQP3( N, NR, V, LDV, IWORK(N+1), WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR )
**               CALL DGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1),
**     $              LWORK-2*N, IERR )
               if ( L2PERT ) {
                  XSC = DSQRT(SMALL)
                  DO 3969 p = 2, NR
                     DO 3968 q = 1, p - 1
                        TEMP1 = XSC * MIN(DABS(V(p,p)),DABS(V(q,q)))
                        IF ( DABS(V(q,p)) .LE. TEMP1 ) V(q,p) = DSIGN( TEMP1, V(q,p) )
 3968                CONTINUE
 3969             CONTINUE
               }

               CALL DLACPY( 'A', N, NR, V, LDV, WORK(2*N+1), N )

               if ( L2PERT ) {
                  XSC = DSQRT(SMALL)
                  DO 8970 p = 2, NR
                     DO 8971 q = 1, p - 1
                        TEMP1 = XSC * MIN(DABS(V(p,p)),DABS(V(q,q)))
                        V(p,q) = - DSIGN( TEMP1, V(q,p) )
 8971                CONTINUE
 8970             CONTINUE
               } else {
                  CALL DLASET( 'L',NR-1,NR-1,ZERO,ZERO,V(2,1),LDV )
               }
               // Now, compute R2 = L3 * Q3, the LQ factorization.
               CALL DGELQF( NR, NR, V, LDV, WORK(2*N+N*NR+1), WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR )
               // .. and estimate the condition number
               CALL DLACPY( 'L',NR,NR,V,LDV,WORK(2*N+N*NR+NR+1),NR )
               DO 4950 p = 1, NR
                  TEMP1 = DNRM2( p, WORK(2*N+N*NR+NR+p), NR )
                  CALL DSCAL( p, ONE/TEMP1, WORK(2*N+N*NR+NR+p), NR )
 4950          CONTINUE
               CALL DPOCON( 'L',NR,WORK(2*N+N*NR+NR+1),NR,ONE,TEMP1, WORK(2*N+N*NR+NR+NR*NR+1),IWORK(M+2*N+1),IERR )
               CONDR2 = ONE / DSQRT(TEMP1)

               if ( CONDR2 .GE. COND_OK ) {
                  // .. save the Householder vectors used for Q3
                  // (this overwrites the copy of R2, as it will not be
                  // needed in this branch, but it does not overwrite the
                  // Huseholder vectors of Q2.).
                  CALL DLACPY( 'U', NR, NR, V, LDV, WORK(2*N+1), N )
                  // .. and the rest of the information on Q3 is in
                  // WORK(2*N+N*NR+1:2*N+N*NR+N)
               }

            }

            if ( L2PERT ) {
               XSC = DSQRT(SMALL)
               DO 4968 q = 2, NR
                  TEMP1 = XSC * V(q,q)
                  DO 4969 p = 1, q - 1
                     // V(p,q) = - DSIGN( TEMP1, V(q,p) )
                     V(p,q) = - DSIGN( TEMP1, V(p,q) )
 4969             CONTINUE
 4968          CONTINUE
            } else {
               CALL DLASET( 'U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV )
            }

         // Second preconditioning finished; continue with Jacobi SVD
         // The input matrix is lower triangular.

         // Recover the right singular vectors as solution of a well
         // conditioned triangular matrix equation.

            if ( CONDR1 .LT. COND_OK ) {

               CALL DGESVJ( 'L','U','N',NR,NR,V,LDV,SVA,NR,U, LDU,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,INFO )
               SCALEM  = WORK(2*N+N*NR+NR+1)
               NUMRANK = IDNINT(WORK(2*N+N*NR+NR+2))
               DO 3970 p = 1, NR
                  CALL DCOPY( NR, V(1,p), 1, U(1,p), 1 )
                  CALL DSCAL( NR, SVA(p),    V(1,p), 1 )
 3970          CONTINUE

         // .. pick the right matrix equation and solve it

               if ( NR .EQ. N ) {
* :))             .. best case, R1 is inverted. The solution of this matrix
                  // equation is Q2*V2 = the product of the Jacobi rotations
                  // used in DGESVJ, premultiplied with the orthogonal matrix
                  // from the second QR factorization.
                  CALL DTRSM( 'L','U','N','N', NR,NR,ONE, A,LDA, V,LDV )
               } else {
                  // .. R1 is well conditioned, but non-square. Transpose(R2)
                  // is inverted to get the product of the Jacobi rotations
                  // used in DGESVJ. The Q-factor from the second QR
                  // factorization is then built in explicitly.
                  CALL DTRSM('L','U','T','N',NR,NR,ONE,WORK(2*N+1), N,V,LDV)
                  if ( NR .LT. N ) {
                    CALL DLASET('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV)
                    CALL DLASET('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV)
                    CALL DLASET('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV)
                  }
                  CALL DORMQR('L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR)
               }

            } else if ( CONDR2 .LT. COND_OK ) {

* :)           .. the input matrix A is very likely a relative of
              t // he Kahan matrix :)
               // The matrix R2 is inverted. The solution of the matrix equation
               // is Q3^T*V3 = the product of the Jacobi rotations (applied to
              t // he lower triangular L3 from the LQ factorization of
               // R2=L3*Q3), pre-multiplied with the transposed Q3.
               CALL DGESVJ( 'L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU, WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, INFO )
               SCALEM  = WORK(2*N+N*NR+NR+1)
               NUMRANK = IDNINT(WORK(2*N+N*NR+NR+2))
               DO 3870 p = 1, NR
                  CALL DCOPY( NR, V(1,p), 1, U(1,p), 1 )
                  CALL DSCAL( NR, SVA(p),    U(1,p), 1 )
 3870          CONTINUE
               CALL DTRSM('L','U','N','N',NR,NR,ONE,WORK(2*N+1),N,U,LDU)
               // .. apply the permutation from the second QR factorization
               DO 873 q = 1, NR
                  DO 872 p = 1, NR
                     WORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q)
 872              CONTINUE
                  DO 874 p = 1, NR
                     U(p,q) = WORK(2*N+N*NR+NR+p)
 874              CONTINUE
 873           CONTINUE
               if ( NR .LT. N ) {
                  CALL DLASET( 'A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV )
                  CALL DLASET( 'A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV )
                  CALL DLASET( 'A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV )
               }
               CALL DORMQR( 'L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
            } else {
               // Last line of defense.
* #:(          This is a rather pathological case: no scaled condition
               // improvement after two pivoted QR factorizations. Other
               // possibility is that the rank revealing QR factorization
               // or the condition estimator has failed, or the COND_OK
               // is set very close to ONE (which is unnecessary). Normally,
              t // his branch should never be executed, but in rare cases of
               // failure of the RRQR or condition estimator, the last line of
               // defense ensures that DGEJSV completes the task.
               // Compute the full SVD of L3 using DGESVJ with explicit
               // accumulation of Jacobi rotations.
               CALL DGESVJ( 'L', 'U', 'V', NR, NR, V, LDV, SVA, NR, U, LDU, WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, INFO )
               SCALEM  = WORK(2*N+N*NR+NR+1)
               NUMRANK = IDNINT(WORK(2*N+N*NR+NR+2))
               if ( NR .LT. N ) {
                  CALL DLASET( 'A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV )
                  CALL DLASET( 'A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV )
                  CALL DLASET( 'A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV )
               }
               CALL DORMQR( 'L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )

               CALL DORMLQ( 'L', 'T', NR, NR, NR, WORK(2*N+1), N, WORK(2*N+N*NR+1), U, LDU, WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR )
               DO 773 q = 1, NR
                  DO 772 p = 1, NR
                     WORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q)
 772              CONTINUE
                  DO 774 p = 1, NR
                     U(p,q) = WORK(2*N+N*NR+NR+p)
 774              CONTINUE
 773           CONTINUE

            }

            // Permute the rows of V using the (column) permutation from the
            // first QRF. Also, scale the columns to make them unit in
            // Euclidean norm. This applies to all cases.

            TEMP1 = DSQRT(DBLE(N)) * EPSLN
            DO 1972 q = 1, N
               DO 972 p = 1, N
                  WORK(2*N+N*NR+NR+IWORK(p)) = V(p,q)
  972          CONTINUE
               DO 973 p = 1, N
                  V(p,q) = WORK(2*N+N*NR+NR+p)
  973          CONTINUE
               XSC = ONE / DNRM2( N, V(1,q), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL DSCAL( N, XSC, V(1,q), 1 )
 1972       CONTINUE
            // At this moment, V contains the right singular vectors of A.
            // Next, assemble the left singular vector matrix U (M x N).
            if ( NR .LT. M ) {
               CALL DLASET( 'A', M-NR, NR, ZERO, ZERO, U(NR+1,1), LDU )
               if ( NR .LT. N1 ) {
                  CALL DLASET('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU)
                  CALL DLASET('A',M-NR,N1-NR,ZERO,ONE,U(NR+1,NR+1),LDU)
               }
            }

            // The Q matrix from the first QRF is built into the left singular
            // matrix U. This applies to all cases.

            CALL DORMQR( 'Left', 'No_Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR )

            // The columns of U are normalized. The cost is O(M*N) flops.
            TEMP1 = DSQRT(DBLE(M)) * EPSLN
            DO 1973 p = 1, NR
               XSC = ONE / DNRM2( M, U(1,p), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL DSCAL( M, XSC, U(1,p), 1 )
 1973       CONTINUE

            // If the initial QRF is computed with row pivoting, the left
            // singular vectors must be adjusted.

            IF ( ROWPIV ) CALL DLASWP( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 )

         } else {

         // .. the initial matrix A has almost orthogonal columns and
        t // he second QRF is not needed

            CALL DLACPY( 'Upper', N, N, A, LDA, WORK(N+1), N )
            if ( L2PERT ) {
               XSC = DSQRT(SMALL)
               DO 5970 p = 2, N
                  TEMP1 = XSC * WORK( N + (p-1)*N + p )
                  DO 5971 q = 1, p - 1
                     WORK(N+(q-1)*N+p)=-DSIGN(TEMP1,WORK(N+(p-1)*N+q))
 5971             CONTINUE
 5970          CONTINUE
            } else {
               CALL DLASET( 'Lower',N-1,N-1,ZERO,ZERO,WORK(N+2),N )
            }

            CALL DGESVJ( 'Upper', 'U', 'N', N, N, WORK(N+1), N, SVA, N, U, LDU, WORK(N+N*N+1), LWORK-N-N*N, INFO )

            SCALEM  = WORK(N+N*N+1)
            NUMRANK = IDNINT(WORK(N+N*N+2))
            DO 6970 p = 1, N
               CALL DCOPY( N, WORK(N+(p-1)*N+1), 1, U(1,p), 1 )
               CALL DSCAL( N, SVA(p), WORK(N+(p-1)*N+1), 1 )
 6970       CONTINUE

            CALL DTRSM( 'Left', 'Upper', 'NoTrans', 'No UD', N, N, ONE, A, LDA, WORK(N+1), N )
            DO 6972 p = 1, N
               CALL DCOPY( N, WORK(N+p), N, V(IWORK(p),1), LDV )
 6972       CONTINUE
            TEMP1 = DSQRT(DBLE(N))*EPSLN
            DO 6971 p = 1, N
               XSC = ONE / DNRM2( N, V(1,p), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL DSCAL( N, XSC, V(1,p), 1 )
 6971       CONTINUE

            // Assemble the left singular vector matrix U (M x N).

            if ( N .LT. M ) {
               CALL DLASET( 'A',  M-N, N, ZERO, ZERO, U(N+1,1), LDU )
               if ( N .LT. N1 ) {
                  CALL DLASET( 'A',N,  N1-N, ZERO, ZERO,  U(1,N+1),LDU )
                  CALL DLASET( 'A',M-N,N1-N, ZERO, ONE,U(N+1,N+1),LDU )
               }
            }
            CALL DORMQR( 'Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR )
            TEMP1 = DSQRT(DBLE(M))*EPSLN
            DO 6973 p = 1, N1
               XSC = ONE / DNRM2( M, U(1,p), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL DSCAL( M, XSC, U(1,p), 1 )
 6973       CONTINUE

            IF ( ROWPIV ) CALL DLASWP( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 )

         }

         // end of the  >> almost orthogonal case <<  in the full SVD

         } else {

         // This branch deploys a preconditioned Jacobi SVD with explicitly
         // accumulated rotations. It is included as optional, mainly for
         // experimental purposes. It does perform well, and can also be used.
         // In this implementation, this branch will be automatically activated
         // if the  condition number sigma_max(A) / sigma_min(A) is predicted
        t // o be greater than the overflow threshold. This is because the
         // a posteriori computation of the singular vectors assumes robust
         // implementation of BLAS and some LAPACK procedures, capable of working
         // in presence of extreme values. Since that is not always the case, ...

         DO 7968 p = 1, NR
            CALL DCOPY( N-p+1, A(p,p), LDA, V(p,p), 1 )
 7968    CONTINUE

         if ( L2PERT ) {
            XSC = DSQRT(SMALL/EPSLN)
            DO 5969 q = 1, NR
               TEMP1 = XSC*DABS( V(q,q) )
               DO 5968 p = 1, N
                  IF ( ( p .GT. q ) .AND. ( DABS(V(p,q)) .LE. TEMP1 ) .OR. ( p .LT. q ) ) V(p,q) = DSIGN( TEMP1, V(p,q) )
                  IF ( p .LT. q ) V(p,q) = - V(p,q)
 5968          CONTINUE
 5969       CONTINUE
         } else {
            CALL DLASET( 'U', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV )
         }
          CALL DGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR )
         CALL DLACPY( 'L', N, NR, V, LDV, WORK(2*N+1), N )

         DO 7969 p = 1, NR
            CALL DCOPY( NR-p+1, V(p,p), LDV, U(p,p), 1 )
 7969    CONTINUE

         if ( L2PERT ) {
            XSC = DSQRT(SMALL/EPSLN)
            DO 9970 q = 2, NR
               DO 9971 p = 1, q - 1
                  TEMP1 = XSC * MIN(DABS(U(p,p)),DABS(U(q,q)))
                  U(p,q) = - DSIGN( TEMP1, U(q,p) )
 9971          CONTINUE
 9970       CONTINUE
         } else {
            CALL DLASET('U', NR-1, NR-1, ZERO, ZERO, U(1,2), LDU )
         }
          CALL DGESVJ( 'G', 'U', 'V', NR, NR, U, LDU, SVA, N, V, LDV, WORK(2*N+N*NR+1), LWORK-2*N-N*NR, INFO )
         SCALEM  = WORK(2*N+N*NR+1)
         NUMRANK = IDNINT(WORK(2*N+N*NR+2))

         if ( NR .LT. N ) {
            CALL DLASET( 'A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV )
            CALL DLASET( 'A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV )
            CALL DLASET( 'A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV )
         }
          CALL DORMQR( 'L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )

            // Permute the rows of V using the (column) permutation from the
            // first QRF. Also, scale the columns to make them unit in
            // Euclidean norm. This applies to all cases.

            TEMP1 = DSQRT(DBLE(N)) * EPSLN
            DO 7972 q = 1, N
               DO 8972 p = 1, N
                  WORK(2*N+N*NR+NR+IWORK(p)) = V(p,q)
 8972          CONTINUE
               DO 8973 p = 1, N
                  V(p,q) = WORK(2*N+N*NR+NR+p)
 8973          CONTINUE
               XSC = ONE / DNRM2( N, V(1,q), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL DSCAL( N, XSC, V(1,q), 1 )
 7972       CONTINUE

            // At this moment, V contains the right singular vectors of A.
            // Next, assemble the left singular vector matrix U (M x N).

         if ( NR .LT. M ) {
            CALL DLASET( 'A',  M-NR, NR, ZERO, ZERO, U(NR+1,1), LDU )
            if ( NR .LT. N1 ) {
               CALL DLASET( 'A',NR,  N1-NR, ZERO, ZERO,  U(1,NR+1),LDU )
               CALL DLASET( 'A',M-NR,N1-NR, ZERO, ONE,U(NR+1,NR+1),LDU )
            }
         }

         CALL DORMQR( 'Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR )

            IF ( ROWPIV ) CALL DLASWP( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 )


         }
         if ( TRANSP ) {
            // .. swap U and V because the procedure worked on A^t
            DO 6974 p = 1, N
               CALL DSWAP( N, U(1,p), 1, V(1,p), 1 )
 6974       CONTINUE
         }

      }
      // end of the full SVD

      // Undo scaling, if necessary (and possible)

      if ( USCAL2 .LE. (BIG/SVA(1))*USCAL1 ) {
         CALL DLASCL( 'G', 0, 0, USCAL1, USCAL2, NR, 1, SVA, N, IERR )
         USCAL1 = ONE
         USCAL2 = ONE
      }

      if ( NR .LT. N ) {
         DO 3004 p = NR+1, N
            SVA(p) = ZERO
 3004    CONTINUE
      }

      WORK(1) = USCAL2 * SCALEM
      WORK(2) = USCAL1
      IF ( ERREST ) WORK(3) = SCONDA
      if ( LSVEC .AND. RSVEC ) {
         WORK(4) = CONDR1
         WORK(5) = CONDR2
      }
      if ( L2TRAN ) {
         WORK(6) = ENTRA
         WORK(7) = ENTRAT
      }

      IWORK(1) = NR
      IWORK(2) = NUMRANK
      IWORK(3) = WARNING

      RETURN
      // ..
      // .. END OF DGEJSV
      // ..
      }
