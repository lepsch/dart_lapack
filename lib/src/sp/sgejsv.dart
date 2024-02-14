      void sgejsv(final int JOBA, final int JOBU, final int JOBV, final int JOBR, final int JOBT, final int JOBP, final int M, final int N, final Matrix<double> A_, final int LDA, final int SVA, final Matrix<double> U_, final int LDU, final Matrix<double> V_, final int LDV, final Array<double> WORK_, final int LWORK, final Array<int> IWORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final U = U_.dim();
  final V = V_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int         INFO, LDA, LDU, LDV, LWORK, M, N;
      double        A( LDA, * ), SVA( N ), U( LDU, * ), V( LDV, * ), WORK( LWORK );
      int         IWORK( * );
      String      JOBA, JOBP, JOBR, JOBT, JOBU, JOBV;
      // ..

// ===========================================================================

      // .. Local Parameters ..
      double        ZERO,         ONE;
      const     ZERO = 0.0, ONE = 1.0 ;
      double    AAPP,   AAQQ,   AATMAX, AATMIN, BIG,    BIG1,   COND_OK, CONDR1, CONDR2, ENTRA,  ENTRAT, EPSLN,  MAXPRJ, SCALEM, SCONDA, SFMIN,  SMALL,  TEMP1,  USCAL1, USCAL2, XSC;
      int     IERR,   N1,     NR,     NUMRANK,        p, q,   WARNING;
      bool    ALMORT, DEFR,   ERREST, GOSCAL, JRACC,  KILL,   LSVEC, L2ABER, L2KILL, L2PERT, L2RANK, L2TRAN, NOSCAL, ROWPIV, RSVEC,  TRANSP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, ALOG, MAX, MIN, FLOAT, NINT, SIGN, SQRT
      // ..
      // .. External Functions ..
      //- REAL      SLAMCH, SNRM2;
      //- int       ISAMAX;
      //- bool      lsame;
      // EXTERNAL ISAMAX, lsame, SLAMCH, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY,  SGELQF, SGEQP3, SGEQRF, SLACPY, SLASCL, SLASET, SLASSQ, SLASWP, SORGQR, SORMLQ, SORMQR, SPOCON, SSCAL,  SSWAP,  STRSM,  XERBLA

      // EXTERNAL SGESVJ
      // ..

      // Test the input arguments

      LSVEC  = lsame( JOBU, 'U' ) || lsame( JOBU, 'F' );
      JRACC  = lsame( JOBV, 'J' );
      RSVEC  = lsame( JOBV, 'V' ) || JRACC;
      ROWPIV = lsame( JOBA, 'F' ) || lsame( JOBA, 'G' );
      L2RANK = lsame( JOBA, 'R' );
      L2ABER = lsame( JOBA, 'A' );
      ERREST = lsame( JOBA, 'E' ) || lsame( JOBA, 'G' );
      L2TRAN = lsame( JOBT, 'T' );
      L2KILL = lsame( JOBR, 'R' );
      DEFR   = lsame( JOBR, 'N' );
      L2PERT = lsame( JOBP, 'P' );

      if ( !(ROWPIV || L2RANK || L2ABER || ERREST || lsame( JOBA, 'C' ) )) {
         INFO = - 1;
      } else if ( !( LSVEC || lsame( JOBU, 'N' ) || lsame( JOBU, 'W' )) ) {
         INFO = - 2;
      } else if ( !( RSVEC || lsame( JOBV, 'N' ) || lsame( JOBV, 'W' )) || ( JRACC && ( !LSVEC) ) ) {
         INFO = - 3;
      } else if ( !( L2KILL || DEFR ) ) {
         INFO = - 4;
      } else if ( !( L2TRAN || lsame( JOBT, 'N' ) ) ) {
         INFO = - 5;
      } else if ( !( L2PERT || lsame( JOBP, 'N' ) ) ) {
         INFO = - 6;
      } else if ( M < 0 ) {
         INFO = - 7;
      } else if ( ( N < 0 ) || ( N > M ) ) {
         INFO = - 8;
      } else if ( LDA < M ) {
         INFO = - 10;
      } else if ( LSVEC && ( LDU < M ) ) {
         INFO = - 13;
      } else if ( RSVEC && ( LDV < N ) ) {
         INFO = - 15;
      } else if ( ( !(LSVEC || RSVEC || ERREST) && (LWORK < max(7,4*N+1,2*M+N))) || ( !(LSVEC || RSVEC) && ERREST && (LWORK < max(7,4*N+N*N,2*M+N))) || (LSVEC && ( !RSVEC) && (LWORK < max(7,2*M+N,4*N+1))) || (RSVEC && ( !LSVEC) && (LWORK < max(7,2*M+N,4*N+1))) || (LSVEC && RSVEC && ( !JRACC) && (LWORK < max(2*M+N,6*N+2*N*N))) || (LSVEC && RSVEC && JRACC && LWORK < max(2*M+N,4*N+N*N,2*N+N*N+6))) {
         INFO = - 17;
      } else {
         // #:)
         INFO = 0;
      }

      if ( INFO != 0 ) {
        // #:(
         xerbla('SGEJSV', - INFO );
         return;
      }

      // Quick return for void matrix (Y3K safe)
// #:)
      if ( ( M == 0 ) || ( N == 0 ) ) {
         IWORK[1:3] = 0;
         WORK[1:7] = 0;
         return;
      }

      // Determine whether the matrix U should be M x N or M x M

      if ( LSVEC ) {
         N1 = N;
         if ( lsame( JOBU, 'F' ) ) N1 = M;
      }

      // Set numerical parameters

// !    NOTE: Make sure SLAMCH() does not fail on the target architecture.

      EPSLN = SLAMCH('Epsilon');
      SFMIN = SLAMCH('SafeMinimum');
      SMALL = SFMIN / EPSLN;
      BIG   = SLAMCH('O');
      // BIG   = ONE / SFMIN

      // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N

// (!)  If necessary, scale SVA() to protect the largest norm from
      // overflow. It is possible that this scaling pushes the smallest
      // column norm left from the underflow threshold (extreme case).

      SCALEM  = ONE / sqrt(FLOAT(M)*FLOAT(N));
      NOSCAL  = true;
      GOSCAL  = true;
      for (p = 1; p <= N; p++) { // 1874
         AAPP = ZERO;
         AAQQ = ONE;
         slassq(M, A(1,p), 1, AAPP, AAQQ );
         if ( AAPP > BIG ) {
            INFO = - 9;
            xerbla('SGEJSV', -INFO );
            return;
         }
         AAQQ = sqrt(AAQQ);
         if ( ( AAPP < (BIG / AAQQ) ) && NOSCAL  ) {
            SVA[p] = AAPP * AAQQ;
         } else {
            NOSCAL  = false;
            SVA[p] = AAPP * ( AAQQ * SCALEM );
            if ( GOSCAL ) {
               GOSCAL = false;
               sscal(p-1, SCALEM, SVA, 1 );
            }
         }
      } // 1874

      if (NOSCAL) SCALEM = ONE;

      AAPP = ZERO;
      AAQQ = BIG;
      for (p = 1; p <= N; p++) { // 4781
         AAPP = max( AAPP, SVA(p) );
         if ( SVA(p) != ZERO ) AAQQ = min( AAQQ, SVA(p) );
      } // 4781

      // Quick return for zero M x N matrix
// #:)
      if ( AAPP == ZERO ) {
         if (LSVEC) slaset( 'G', M, N1, ZERO, ONE, U, LDU );
         if (RSVEC) slaset( 'G', N, N,  ZERO, ONE, V, LDV );
         WORK[1] = ONE;
         WORK[2] = ONE;
         if (ERREST) WORK(3) = ONE;
         if ( LSVEC && RSVEC ) {
            WORK[4] = ONE;
            WORK[5] = ONE;
         }
         if ( L2TRAN ) {
            WORK[6] = ZERO;
            WORK[7] = ZERO;
         }
         IWORK[1] = 0;
         IWORK[2] = 0;
         IWORK[3] = 0;
         return;
      }

      // Issue warning if denormalized column norms detected. Override the
      // high relative accuracy request. Issue licence to kill columns
      // (set them to zero) whose norm is less than sigma_max / BIG (roughly).
// #:(
      WARNING = 0;
      if ( AAQQ <= SFMIN ) {
         L2RANK = true;
         L2KILL = true;
         WARNING = 1;
      }

      // Quick return for one-column matrix
// #:)
      if ( N == 1 ) {

         if ( LSVEC ) {
            slascl('G',0,0,SVA(1),SCALEM, M,1,A(1,1),LDA,IERR );
            slacpy('A', M, 1, A, LDA, U, LDU );
            // computing all M left singular vectors of the M x 1 matrix
            if ( N1 != N  ) {
               sgeqrf(M, N, U,LDU, WORK, WORK(N+1),LWORK-N,IERR );
               sorgqr(M,N1,1, U,LDU,WORK,WORK(N+1),LWORK-N,IERR );
               scopy(M, A(1,1), 1, U(1,1), 1 );
            }
         }
         if ( RSVEC ) {
             V[1][1] = ONE;
         }
         if ( SVA(1) < (BIG*SCALEM) ) {
            SVA[1] = SVA(1) / SCALEM;
            SCALEM  = ONE;
         }
         WORK[1] = ONE / SCALEM;
         WORK[2] = ONE;
         if ( SVA(1) != ZERO ) {
            IWORK[1] = 1;
            if ( ( SVA(1) / SCALEM) >= SFMIN ) {
               IWORK[2] = 1;
            } else {
               IWORK[2] = 0;
            }
         } else {
            IWORK[1] = 0;
            IWORK[2] = 0;
         }
         IWORK[3] = 0;
         if (ERREST) WORK(3) = ONE;
         if ( LSVEC && RSVEC ) {
            WORK[4] = ONE;
            WORK[5] = ONE;
         }
         if ( L2TRAN ) {
            WORK[6] = ZERO;
            WORK[7] = ZERO;
         }
         return;

      }

      TRANSP = false;
      L2TRAN = L2TRAN && ( M == N );

      AATMAX = -ONE;
      AATMIN =  BIG;
      if ( ROWPIV || L2TRAN ) {

      // Compute the row norms, needed to determine row pivoting sequence
      // (in the case of heavily row weighted A, row pivoting is strongly
      // advised) and to collect information needed to compare the
      // structures of A * A^t and A^t * A (in the case L2TRAN == true ).

         if ( L2TRAN ) {
            for (p = 1; p <= M; p++) { // 1950
               XSC   = ZERO;
               TEMP1 = ONE;
               slassq(N, A(p,1), LDA, XSC, TEMP1 );
               // SLASSQ gets both the ell_2 and the ell_infinity norm
               // in one pass through the vector
               WORK[M+N+p] = XSC * SCALEM;
               WORK[N+p] = XSC * (SCALEM*sqrt(TEMP1));
               AATMAX = max( AATMAX, WORK(N+p) );
               if (WORK(N+p) != ZERO) AATMIN = min(AATMIN,WORK(N+p));
            } // 1950
         } else {
            for (p = 1; p <= M; p++) { // 1904
               WORK[M+N+p] = SCALEM*ABS( A(p,ISAMAX(N,A(p,1),LDA)) );
               AATMAX = max( AATMAX, WORK(M+N+p) );
               AATMIN = min( AATMIN, WORK(M+N+p) );
            } // 1904
         }

      }

      // For square matrix A try to determine whether A^t  would be  better
      // input for the preconditioned Jacobi SVD, with faster convergence.
      // The decision is based on an O(N) function of the vector of column
      // and row norms of A, based on the Shannon entropy. This should give
      // the right choice in most cases when the difference actually matters.
      // It may fail and pick the slower converging side.

      ENTRA  = ZERO;
      ENTRAT = ZERO;
      if ( L2TRAN ) {

         XSC   = ZERO;
         TEMP1 = ONE;
         slassq(N, SVA, 1, XSC, TEMP1 );
         TEMP1 = ONE / TEMP1;

         ENTRA = ZERO;
         for (p = 1; p <= N; p++) { // 1113
            BIG1  = ( ( SVA(p) / XSC )**2 ) * TEMP1;
            if (BIG1 != ZERO) ENTRA = ENTRA + BIG1 * ALOG(BIG1);
         } // 1113
         ENTRA = - ENTRA / ALOG(FLOAT(N));

         // Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex.
         // It is derived from the diagonal of  A^t * A.  Do the same with the
         // diagonal of A * A^t, compute the entropy of the corresponding
         // probability distribution. Note that A * A^t and A^t * A have the
         // same trace.

         ENTRAT = ZERO;
         for (p = N+1; p <= N+M; p++) { // 1114
            BIG1 = ( ( WORK(p) / XSC )**2 ) * TEMP1;
            if (BIG1 != ZERO) ENTRAT = ENTRAT + BIG1 * ALOG(BIG1);
         } // 1114
         ENTRAT = - ENTRAT / ALOG(FLOAT(M));

         // Analyze the entropies and decide A or A^t. Smaller entropy
         // usually means better input for the algorithm.

         TRANSP = ( ENTRAT < ENTRA );

         // If A^t is better than A, transpose A.

         if ( TRANSP ) {
            // In an optimal implementation, this trivial transpose
            // should be replaced with faster transpose.
            for (p = 1; p <= N - 1; p++) { // 1115
               for (q = p + 1; q <= N; q++) { // 1116
                   TEMP1 = A(q,p);
                  A[q][p] = A(p,q);
                  A[p][q] = TEMP1;
               } // 1116
            } // 1115
            for (p = 1; p <= N; p++) { // 1117
               WORK[M+N+p] = SVA(p);
               SVA[p] = WORK(N+p);
            } // 1117
            TEMP1  = AAPP;
            AAPP   = AATMAX;
            AATMAX = TEMP1;
            TEMP1  = AAQQ;
            AAQQ   = AATMIN;
            AATMIN = TEMP1;
            KILL   = LSVEC;
            LSVEC  = RSVEC;
            RSVEC  = KILL;
            if (LSVEC) N1 = N;

            ROWPIV = true;
         }

      }
      // END IF L2TRAN

      // Scale the matrix so that its maximal singular value remains less
      // than sqrt(BIG) -- the matrix is scaled so that its maximal column
      // has Euclidean norm equal to sqrt(BIG/N). The only reason to keep
      // sqrt(BIG) instead of BIG is the fact that SGEJSV uses LAPACK and
      // BLAS routines that, in some implementations, are not capable of
      // working in the full interval [SFMIN,BIG] and that they may provoke
      // overflows in the intermediate results. If the singular values spread
      // from SFMIN to BIG, then SGESVJ will compute them. So, in that case,
      // one should use SGESVJ instead of SGEJSV.

      BIG1   = sqrt( BIG );
      TEMP1  = sqrt( BIG / FLOAT(N) );

      slascl('G', 0, 0, AAPP, TEMP1, N, 1, SVA, N, IERR );
      if ( AAQQ > (AAPP * SFMIN) ) {
          AAQQ = ( AAQQ / AAPP ) * TEMP1;
      } else {
          AAQQ = ( AAQQ * TEMP1 ) / AAPP;
      }
      TEMP1 = TEMP1 * SCALEM;
      slascl('G', 0, 0, AAPP, TEMP1, M, N, A, LDA, IERR );

      // To undo scaling at the end of this procedure, multiply the
      // computed singular values with USCAL2 / USCAL1.

      USCAL1 = TEMP1;
      USCAL2 = AAPP;

      if ( L2KILL ) {
         // L2KILL enforces computation of nonzero singular values in
         // the restricted range of condition number of the initial A,
         // sigma_max(A) / sigma_min(A) approx. sqrt(BIG)/sqrt(SFMIN).
         XSC = sqrt( SFMIN );
      } else {
         XSC = SMALL;

         // Now, if the condition number of A is too big,
         // sigma_max(A) / sigma_min(A) > sqrt(BIG/N) * EPSLN / SFMIN,
         // as a precaution measure, the full SVD is computed using SGESVJ
         // with accumulated Jacobi rotations. This provides numerically
         // more robust computation, at the cost of slightly increased run
         // time. Depending on the concrete implementation of BLAS and LAPACK
         // (i.e. how they behave in presence of extreme ill-conditioning) the
         // implementor may decide to remove this switch.
         if ( ( AAQQ < sqrt(SFMIN) ) && LSVEC && RSVEC ) {
            JRACC = true;
         }

      }
      if ( AAQQ < XSC ) {
         for (p = 1; p <= N; p++) { // 700
            if ( SVA(p) < XSC ) {
               slaset('A', M, 1, ZERO, ZERO, A(1,p), LDA );
               SVA[p] = ZERO;
            }
         } // 700
      }

      // Preconditioning using QR factorization with pivoting

      if ( ROWPIV ) {
         // Optional row permutation (Bjoerck row pivoting):
         // A result by Cox and Higham shows that the Bjoerck's
         // row pivoting combined with standard column pivoting
         // has similar effect as Powell-Reid complete pivoting.
         // The ell-infinity norms of A are made nonincreasing.
         for (p = 1; p <= M - 1; p++) { // 1952
            q = ISAMAX( M-p+1, WORK(M+N+p), 1 ) + p - 1;
            IWORK[2*N+p] = q;
            if ( p != q ) {
               TEMP1       = WORK(M+N+p);
               WORK[M+N+p] = WORK(M+N+q);
               WORK[M+N+q] = TEMP1;
            }
         } // 1952
         slaswp(N, A, LDA, 1, M-1, IWORK(2*N+1), 1 );
      }

      // End of the preparation phase (scaling, optional sorting and
      // transposing, optional flushing of small columns).

      // Preconditioning

      // If the full SVD is needed, the right singular vectors are computed
      // from a matrix equation, and for that we need theoretical analysis
      // of the Businger-Golub pivoting. So we use SGEQP3 as the first RR QRF.
      // In all other cases the first RR QRF can be chosen by other criteria
      // (eg speed by replacing global with restricted window pivoting, such
      // as in SGEQPX from TOMS # 782). Good results will be obtained using
      // SGEQPX with properly (!) chosen numerical parameters.
      // Any improvement of SGEQP3 improves overall performance of SGEJSV.

      // A * P1 = Q1 * [ R1^t 0]^t:
      for (p = 1; p <= N; p++) { // 1963
         // .. all columns are free columns
         IWORK[p] = 0;
      } // 1963
      sgeqp3(M,N,A,LDA, IWORK,WORK, WORK(N+1),LWORK-N, IERR );

      // The upper triangular matrix R1 from the first QRF is inspected for
      // rank deficiency and possibilities for deflation, or possible
      // ill-conditioning. Depending on the user specified flag L2RANK,
      // the procedure explores possibilities to reduce the numerical
      // rank by inspecting the computed upper triangular factor. If
      // L2RANK or L2ABER are up, then SGEJSV will compute the SVD of
      // A + dA, where ||dA|| <= f(M,N)*EPSLN.

      NR = 1;
      if ( L2ABER ) {
         // Standard absolute error bound suffices. All sigma_i with
         // sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
         // aggressive enforcement of lower numerical rank by introducing a
         // backward error of the order of N*EPSLN*||A||.
         TEMP1 = sqrt(FLOAT(N))*EPSLN;
         for (p = 2; p <= N; p++) { // 3001
            if ( (A(p,p)).abs() >= (TEMP1*(A(1,1) ).abs() ) ) {
               NR = NR + 1;
            } else {
               GO TO 3002;
            }
         } // 3001
         } // 3002
      } else if ( L2RANK ) {
         // .. similarly as above, only slightly more gentle (less aggressive).
         // Sudden drop on the diagonal of R1 is used as the criterion for
         // close-to-rank-deficient.
         TEMP1 = sqrt(SFMIN);
         for (p = 2; p <= N; p++) { // 3401
            if ( ( (A(p,p)).abs() < (EPSLN*(A(p-1,p-1) ).abs() ) ) || ( (A(p,p)).abs() < SMALL ) || ( L2KILL && ((A(p,p)).abs() < TEMP1) ) ) GO TO 3402;
            NR = NR + 1;
         } // 3401
         } // 3402

      } else {
         // The goal is high relative accuracy. However, if the matrix
         // has high scaled condition number the relative accuracy is in
         // general not feasible. Later on, a condition number estimator
         // will be deployed to estimate the scaled condition number.
         // Here we just remove the underflowed part of the triangular
         // factor. This prevents the situation in which the code is
         // working hard to get the accuracy not warranted by the data.
         TEMP1  = sqrt(SFMIN);
         for (p = 2; p <= N; p++) { // 3301
            if ( ( (A(p,p)).abs() < SMALL ) || ( L2KILL && ((A(p,p)).abs() < TEMP1) ) ) GO TO 3302;
            NR = NR + 1;
         } // 3301
         } // 3302

      }

      ALMORT = false;
      if ( NR == N ) {
         MAXPRJ = ONE;
         for (p = 2; p <= N; p++) { // 3051
            TEMP1  = (A(p,p)).abs() / SVA(IWORK(p));
            MAXPRJ = min( MAXPRJ, TEMP1 );
         } // 3051
         if ( MAXPRJ**2 >= ONE - FLOAT(N)*EPSLN ) ALMORT = true;
      }


      SCONDA = - ONE;
      CONDR1 = - ONE;
      CONDR2 = - ONE;

      if ( ERREST ) {
         if ( N == NR ) {
            if ( RSVEC ) {
               // .. V is available as workspace
               slacpy('U', N, N, A, LDA, V, LDV );
               for (p = 1; p <= N; p++) { // 3053
                  TEMP1 = SVA(IWORK(p));
                  sscal(p, ONE/TEMP1, V(1,p), 1 );
               } // 3053
               spocon('U', N, V, LDV, ONE, TEMP1, WORK(N+1), IWORK(2*N+M+1), IERR );
            } else if ( LSVEC ) {
               // .. U is available as workspace
               slacpy('U', N, N, A, LDA, U, LDU );
               for (p = 1; p <= N; p++) { // 3054
                  TEMP1 = SVA(IWORK(p));
                  sscal(p, ONE/TEMP1, U(1,p), 1 );
               } // 3054
               spocon('U', N, U, LDU, ONE, TEMP1, WORK(N+1), IWORK(2*N+M+1), IERR );
            } else {
               slacpy('U', N, N, A, LDA, WORK(N+1), N );
               for (p = 1; p <= N; p++) { // 3052
                  TEMP1 = SVA(IWORK(p));
                  sscal(p, ONE/TEMP1, WORK(N+(p-1)*N+1), 1 );
               } // 3052
            // .. the columns of R are scaled to have unit Euclidean lengths.
               spocon('U', N, WORK(N+1), N, ONE, TEMP1, WORK(N+N*N+1), IWORK(2*N+M+1), IERR );
            }
            SCONDA = ONE / sqrt(TEMP1);
            // SCONDA is an estimate of sqrt(||(R^t * R)^(-1)||_1).
            // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
         } else {
            SCONDA = - ONE;
         }
      }

      L2PERT = L2PERT && ( ABS( A(1,1)/A(NR,NR) ) > sqrt(BIG1) );
      // If there is no violent scaling, artificial perturbation is not needed.

      // Phase 3:

      if ( !( RSVEC || LSVEC ) ) {

          // Singular Values only

          // .. transpose A(1:NR,1:N)
         for (p = 1; p <= min( N-1, NR ); p++) { // 1946
            scopy(N-p, A(p,p+1), LDA, A(p+1,p), 1 );
         } // 1946

         // The following two DO-loops introduce small relative perturbation
         // into the strict upper triangle of the lower triangular matrix.
         // Small entries below the main diagonal are also changed.
         // This modification is useful if the computing environment does not
         // provide/allow FLUSH TO ZERO underflow, for it prevents many
         // annoying denormalized numbers in case of strongly scaled matrices.
         // The perturbation is structured so that it does not introduce any
         // new perturbation of the singular values, and it does not destroy
         // the job done by the preconditioner.
         // The licence for this perturbation is in the variable L2PERT, which
         // should be false if FLUSH TO ZERO underflow is active.

         if ( !ALMORT ) {

            if ( L2PERT ) {
               // XSC = sqrt(SMALL)
               XSC = EPSLN / FLOAT(N);
               for (q = 1; q <= NR; q++) { // 4947
                  TEMP1 = XSC*(A(q,q)).abs();
                  for (p = 1; p <= N; p++) { // 4949
                     if ( ( (p > q) && ((A(p,q)).abs() <= TEMP1) ) || ( p < q ) ) A(p,q) = sign( TEMP1, A(p,q) );
                  } // 4949
               } // 4947
            } else {
               slaset('U', NR-1,NR-1, ZERO,ZERO, A(1,2),LDA );
            }

             // .. second preconditioning using the QR factorization

            sgeqrf(N,NR, A,LDA, WORK, WORK(N+1),LWORK-N, IERR );

            // .. and transpose upper to lower triangular
            for (p = 1; p <= NR - 1; p++) { // 1948
               scopy(NR-p, A(p,p+1), LDA, A(p+1,p), 1 );
            } // 1948

         }

            // Row-cyclic Jacobi SVD algorithm with column pivoting

            // .. again some perturbation (a "background noise") is added
            // to drown denormals
            if ( L2PERT ) {
               // XSC = sqrt(SMALL)
               XSC = EPSLN / FLOAT(N);
               for (q = 1; q <= NR; q++) { // 1947
                  TEMP1 = XSC*(A(q,q)).abs();
                  for (p = 1; p <= NR; p++) { // 1949
                     if ( ( (p > q) && ((A(p,q)).abs() <= TEMP1) ) || ( p < q ) ) A(p,q) = sign( TEMP1, A(p,q) );
                  } // 1949
               } // 1947
            } else {
               slaset('U', NR-1, NR-1, ZERO, ZERO, A(1,2), LDA );
            }

            // .. and one-sided Jacobi rotations are started on a lower
            // triangular matrix (plus perturbation which is ignored in
            // the part which destroys triangular form (confusing?!))

            sgesvj('L', 'NoU', 'NoV', NR, NR, A, LDA, SVA, N, V, LDV, WORK, LWORK, INFO );

            SCALEM  = WORK(1);
            NUMRANK = NINT(WORK(2));


      } else if ( RSVEC && ( !LSVEC ) ) {

         // -> Singular Values and Right Singular Vectors <-

         if ( ALMORT ) {

            // .. in this case NR equals N
            for (p = 1; p <= NR; p++) { // 1998
               scopy(N-p+1, A(p,p), LDA, V(p,p), 1 );
            } // 1998
            slaset('Upper', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV );

            sgesvj('L','U','N', N, NR, V,LDV, SVA, NR, A,LDA, WORK, LWORK, INFO );
            SCALEM  = WORK(1);
            NUMRANK = NINT(WORK(2));

         } else {

         // .. two more QR factorizations ( one QRF is not enough, two require
         // accumulated product of Jacobi rotations, three are perfect )

            slaset('Lower', NR-1, NR-1, ZERO, ZERO, A(2,1), LDA );
            sgelqf(NR, N, A, LDA, WORK, WORK(N+1), LWORK-N, IERR);
            slacpy('Lower', NR, NR, A, LDA, V, LDV );
            slaset('Upper', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV );
            sgeqrf(NR, NR, V, LDV, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR );
            for (p = 1; p <= NR; p++) { // 8998
               scopy(NR-p+1, V(p,p), LDV, V(p,p), 1 );
            } // 8998
            slaset('Upper', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV );

            sgesvj('Lower', 'U','N', NR, NR, V,LDV, SVA, NR, U, LDU, WORK(N+1), LWORK-N, INFO );
            SCALEM  = WORK(N+1);
            NUMRANK = NINT(WORK(N+2));
            if ( NR < N ) {
               slaset('A',N-NR, NR, ZERO,ZERO, V(NR+1,1),   LDV );
               slaset('A',NR, N-NR, ZERO,ZERO, V(1,NR+1),   LDV );
               slaset('A',N-NR,N-NR,ZERO,ONE, V(NR+1,NR+1), LDV );
            }

         sormlq('Left', 'Transpose', N, N, NR, A, LDA, WORK, V, LDV, WORK(N+1), LWORK-N, IERR );

         }

         for (p = 1; p <= N; p++) { // 8991
            scopy(N, V(p,1), LDV, A(IWORK(p),1), LDA );
         } // 8991
         slacpy('All', N, N, A, LDA, V, LDV );

         if ( TRANSP ) {
            slacpy('All', N, N, V, LDV, U, LDU );
         }

      } else if ( LSVEC && ( !RSVEC ) ) {

         // .. Singular Values and Left Singular Vectors                 ..

         // .. second preconditioning step to avoid need to accumulate
         // Jacobi rotations in the Jacobi iterations.
         for (p = 1; p <= NR; p++) { // 1965
            scopy(N-p+1, A(p,p), LDA, U(p,p), 1 );
         } // 1965
         slaset('Upper', NR-1, NR-1, ZERO, ZERO, U(1,2), LDU );

         sgeqrf(N, NR, U, LDU, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR );

         for (p = 1; p <= NR - 1; p++) { // 1967
            scopy(NR-p, U(p,p+1), LDU, U(p+1,p), 1 );
         } // 1967
         slaset('Upper', NR-1, NR-1, ZERO, ZERO, U(1,2), LDU );

         sgesvj('Lower', 'U', 'N', NR,NR, U, LDU, SVA, NR, A, LDA, WORK(N+1), LWORK-N, INFO );
         SCALEM  = WORK(N+1);
         NUMRANK = NINT(WORK(N+2));

         if ( NR < M ) {
            slaset('A',  M-NR, NR,ZERO, ZERO, U(NR+1,1), LDU );
            if ( NR < N1 ) {
               slaset('A',NR, N1-NR, ZERO, ZERO, U(1,NR+1), LDU );
               slaset('A',M-NR,N1-NR,ZERO,ONE,U(NR+1,NR+1), LDU );
            }
         }

         sormqr('Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );

         if (ROWPIV) slaswp( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 );

         for (p = 1; p <= N1; p++) { // 1974
            XSC = ONE / SNRM2( M, U(1,p), 1 );
            sscal(M, XSC, U(1,p), 1 );
         } // 1974

         if ( TRANSP ) {
            slacpy('All', N, N, U, LDU, V, LDV );
         }

      } else {

         // .. Full SVD ..

         if ( !JRACC ) {

         if ( !ALMORT ) {

            // Second Preconditioning Step (QRF [with pivoting])
            // Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
            // equivalent to an LQF CALL. Since in many libraries the QRF
            // seems to be better optimized than the LQF, we do explicit
            // transpose and use the QRF. This is subject to changes in an
            // optimized implementation of SGEJSV.

            for (p = 1; p <= NR; p++) { // 1968
               scopy(N-p+1, A(p,p), LDA, V(p,p), 1 );
            } // 1968

            // .. the following two loops perturb small entries to avoid
            // denormals in the second QR factorization, where they are
            // as good as zeros. This is done to avoid painfully slow
            // computation with denormals. The relative size of the perturbation
            // is a parameter that can be changed by the implementer.
            // This perturbation device will be obsolete on machines with
            // properly implemented arithmetic.
            // To switch it off, set L2PERT= false To remove it from  the
            // code, remove the action under L2PERT= true , leave the ELSE part.
            // The following two loops should be blocked and fused with the
            // transposed copy above.

            if ( L2PERT ) {
               XSC = sqrt(SMALL);
               for (q = 1; q <= NR; q++) { // 2969
                  TEMP1 = XSC*( V(q,q) ).abs();
                  for (p = 1; p <= N; p++) { // 2968
                     if ( ( p > q ) && ( (V(p,q)).abs() <= TEMP1 ) || ( p < q ) ) V(p,q) = sign( TEMP1, V(p,q) );
                     if (p < q) V(p,q) = - V(p,q);
                  } // 2968
               } // 2969
            } else {
               slaset('U', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV );
            }

            // Estimate the row scaled condition number of R1
            // (If R1 is rectangular, N > NR, then the condition number
            // of the leading NR x NR submatrix is estimated.)

            slacpy('L', NR, NR, V, LDV, WORK(2*N+1), NR );
            for (p = 1; p <= NR; p++) { // 3950
               TEMP1 = SNRM2(NR-p+1,WORK(2*N+(p-1)*NR+p),1);
               sscal(NR-p+1,ONE/TEMP1,WORK(2*N+(p-1)*NR+p),1);
            } // 3950
            spocon('Lower',NR,WORK(2*N+1),NR,ONE,TEMP1, WORK(2*N+NR*NR+1),IWORK(M+2*N+1),IERR);
            CONDR1 = ONE / sqrt(TEMP1);
            // .. here need a second opinion on the condition number
            // .. then assume worst case scenario
            // R1 is OK for inverse <=> CONDR1 < FLOAT(N)
            // more conservative    <=> CONDR1 < sqrt(FLOAT(N))

            COND_OK = sqrt(FLOAT(NR));
// [TP]       COND_OK is a tuning parameter.

            if ( CONDR1 < COND_OK ) {
               // .. the second QRF without pivoting. Note: in an optimized
               // implementation, this QRF should be implemented as the QRF
               // of a lower triangular matrix.
               // R1^t = Q2 * R2
               sgeqrf(N, NR, V, LDV, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR );

               if ( L2PERT ) {
                  XSC = sqrt(SMALL)/EPSLN;
                  for (p = 2; p <= NR; p++) { // 3959
                     for (q = 1; q <= p - 1; q++) { // 3958
                        TEMP1 = XSC * min((V(p,p)).abs(),(V(q,q) ).abs() );
                        if ( (V(q,p)).abs() <= TEMP1 ) V(q,p) = sign( TEMP1, V(q,p) );
                     } // 3958
                  } // 3959
               }

               if (NR != N) slacpy( 'A', N, NR, V, LDV, WORK(2*N+1), N );
               // .. save ...

            // .. this transposed copy should be better than naive
               for (p = 1; p <= NR - 1; p++) { // 1969
                  scopy(NR-p, V(p,p+1), LDV, V(p+1,p), 1 );
               } // 1969

               CONDR2 = CONDR1;

            } else {

               // .. ill-conditioned case: second QRF with pivoting
               // Note that windowed pivoting would be equally good
               // numerically, and more run-time efficient. So, in
               // an optimal implementation, the next call to SGEQP3
               // should be replaced with eg. CALL SGEQPX (ACM TOMS #782)
               // with properly (carefully) chosen parameters.

               // R1^t * P2 = Q2 * R2
               for (p = 1; p <= NR; p++) { // 3003
                  IWORK[N+p] = 0;
               } // 3003
               sgeqp3(N, NR, V, LDV, IWORK(N+1), WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR );
// *               CALL SGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1),
// *     $              LWORK-2*N, IERR )
               if ( L2PERT ) {
                  XSC = sqrt(SMALL);
                  for (p = 2; p <= NR; p++) { // 3969
                     for (q = 1; q <= p - 1; q++) { // 3968
                        TEMP1 = XSC * min((V(p,p)).abs(),(V(q,q) ).abs() );
                        if ( (V(q,p)).abs() <= TEMP1 ) V(q,p) = sign( TEMP1, V(q,p) );
                     } // 3968
                  } // 3969
               }

               slacpy('A', N, NR, V, LDV, WORK(2*N+1), N );

               if ( L2PERT ) {
                  XSC = sqrt(SMALL);
                  for (p = 2; p <= NR; p++) { // 8970
                     for (q = 1; q <= p - 1; q++) { // 8971
                        TEMP1 = XSC * min((V(p,p)).abs(),(V(q,q) ).abs() );
                        V[p][q] = - sign( TEMP1, V(q,p) );
                     } // 8971
                  } // 8970
               } else {
                  slaset('L',NR-1,NR-1,ZERO,ZERO,V(2,1),LDV );
               }
               // Now, compute R2 = L3 * Q3, the LQ factorization.
               sgelqf(NR, NR, V, LDV, WORK(2*N+N*NR+1), WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR );
               // .. and estimate the condition number
               slacpy('L',NR,NR,V,LDV,WORK(2*N+N*NR+NR+1),NR );
               for (p = 1; p <= NR; p++) { // 4950
                  TEMP1 = SNRM2( p, WORK(2*N+N*NR+NR+p), NR );
                  sscal(p, ONE/TEMP1, WORK(2*N+N*NR+NR+p), NR );
               } // 4950
               spocon('L',NR,WORK(2*N+N*NR+NR+1),NR,ONE,TEMP1, WORK(2*N+N*NR+NR+NR*NR+1),IWORK(M+2*N+1),IERR );
               CONDR2 = ONE / sqrt(TEMP1);

               if ( CONDR2 >= COND_OK ) {
                  // .. save the Householder vectors used for Q3
                  // (this overwrites the copy of R2, as it will not be
                  // needed in this branch, but it does not overwrite the
                  // Huseholder vectors of Q2.).
                  slacpy('U', NR, NR, V, LDV, WORK(2*N+1), N );
                  // .. and the rest of the information on Q3 is in
                  // WORK(2*N+N*NR+1:2*N+N*NR+N)
               }

            }

            if ( L2PERT ) {
               XSC = sqrt(SMALL);
               for (q = 2; q <= NR; q++) { // 4968
                  TEMP1 = XSC * V(q,q);
                  for (p = 1; p <= q - 1; p++) { // 4969
                     // V(p,q) = - sign( TEMP1, V(q,p) )
                     V[p][q] = - sign( TEMP1, V(p,q) );
                  } // 4969
               } // 4968
            } else {
               slaset('U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV );
            }

         // Second preconditioning finished; continue with Jacobi SVD
         // The input matrix is lower triangular.

         // Recover the right singular vectors as solution of a well
         // conditioned triangular matrix equation.

            if ( CONDR1 < COND_OK ) {

               sgesvj('L','U','N',NR,NR,V,LDV,SVA,NR,U, LDU,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,INFO );
               SCALEM  = WORK(2*N+N*NR+NR+1);
               NUMRANK = NINT(WORK(2*N+N*NR+NR+2));
               for (p = 1; p <= NR; p++) { // 3970
                  scopy(NR, V(1,p), 1, U(1,p), 1 );
                  sscal(NR, SVA(p),    V(1,p), 1 );
               } // 3970

         // .. pick the right matrix equation and solve it

               if ( NR == N ) {
// :))             .. best case, R1 is inverted. The solution of this matrix
                  // equation is Q2*V2 = the product of the Jacobi rotations
                  // used in SGESVJ, premultiplied with the orthogonal matrix
                  // from the second QR factorization.
                  strsm('L','U','N','N', NR,NR,ONE, A,LDA, V,LDV );
               } else {
                  // .. R1 is well conditioned, but non-square. Transpose(R2)
                  // is inverted to get the product of the Jacobi rotations
                  // used in SGESVJ. The Q-factor from the second QR
                  // factorization is then built in explicitly.
                  strsm('L','U','T','N',NR,NR,ONE,WORK(2*N+1), N,V,LDV);
                  if ( NR < N ) {
                    slaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV);
                    slaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV);
                    slaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV);
                  }
                  sormqr('L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR);
               }

            } else if ( CONDR2 < COND_OK ) {

// :)           .. the input matrix A is very likely a relative of
               // the Kahan matrix :)
               // The matrix R2 is inverted. The solution of the matrix equation
               // is Q3^T*V3 = the product of the Jacobi rotations (applied to
               // the lower triangular L3 from the LQ factorization of
               // R2=L3*Q3), pre-multiplied with the transposed Q3.
               sgesvj('L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU, WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, INFO );
               SCALEM  = WORK(2*N+N*NR+NR+1);
               NUMRANK = NINT(WORK(2*N+N*NR+NR+2));
               for (p = 1; p <= NR; p++) { // 3870
                  scopy(NR, V(1,p), 1, U(1,p), 1 );
                  sscal(NR, SVA(p),    U(1,p), 1 );
               } // 3870
               strsm('L','U','N','N',NR,NR,ONE,WORK(2*N+1),N,U,LDU);
               // .. apply the permutation from the second QR factorization
               for (q = 1; q <= NR; q++) { // 873
                  for (p = 1; p <= NR; p++) { // 872
                     WORK[2*N+N*NR+NR+IWORK(N+p)] = U(p,q);
                  } // 872
                  for (p = 1; p <= NR; p++) { // 874
                     U[p][q] = WORK(2*N+N*NR+NR+p);
                  } // 874
               } // 873
               if ( NR < N ) {
                  slaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV );
                  slaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV );
                  slaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV );
               }
               sormqr('L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR );
            } else {
               // Last line of defense.
// #:(          This is a rather pathological case: no scaled condition
               // improvement after two pivoted QR factorizations. Other
               // possibility is that the rank revealing QR factorization
               // or the condition estimator has failed, or the COND_OK
               // is set very close to ONE (which is unnecessary). Normally,
               // this branch should never be executed, but in rare cases of
               // failure of the RRQR or condition estimator, the last line of
               // defense ensures that SGEJSV completes the task.
               // Compute the full SVD of L3 using SGESVJ with explicit
               // accumulation of Jacobi rotations.
               sgesvj('L', 'U', 'V', NR, NR, V, LDV, SVA, NR, U, LDU, WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, INFO );
               SCALEM  = WORK(2*N+N*NR+NR+1);
               NUMRANK = NINT(WORK(2*N+N*NR+NR+2));
               if ( NR < N ) {
                  slaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV );
                  slaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV );
                  slaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV );
               }
               sormqr('L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR );

               sormlq('L', 'T', NR, NR, NR, WORK(2*N+1), N, WORK(2*N+N*NR+1), U, LDU, WORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR );
               for (q = 1; q <= NR; q++) { // 773
                  for (p = 1; p <= NR; p++) { // 772
                     WORK[2*N+N*NR+NR+IWORK(N+p)] = U(p,q);
                  } // 772
                  for (p = 1; p <= NR; p++) { // 774
                     U[p][q] = WORK(2*N+N*NR+NR+p);
                  } // 774
               } // 773

            }

            // Permute the rows of V using the (column) permutation from the
            // first QRF. Also, scale the columns to make them unit in
            // Euclidean norm. This applies to all cases.

            TEMP1 = sqrt(FLOAT(N)) * EPSLN;
            for (q = 1; q <= N; q++) { // 1972
               for (p = 1; p <= N; p++) { // 972
                  WORK[2*N+N*NR+NR+IWORK(p)] = V(p,q);
               } // 972
               for (p = 1; p <= N; p++) { // 973
                  V[p][q] = WORK(2*N+N*NR+NR+p);
               } // 973
               XSC = ONE / SNRM2( N, V(1,q), 1 );
               if ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) sscal( N, XSC, V(1,q), 1 );
            } // 1972
            // At this moment, V contains the right singular vectors of A.
            // Next, assemble the left singular vector matrix U (M x N).
            if ( NR < M ) {
               slaset('A', M-NR, NR, ZERO, ZERO, U(NR+1,1), LDU );
               if ( NR < N1 ) {
                  slaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                  slaset('A',M-NR,N1-NR,ZERO,ONE,U(NR+1,NR+1),LDU);
               }
            }

            // The Q matrix from the first QRF is built into the left singular
            // matrix U. This applies to all cases.

            sormqr('Left', 'No_Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );

            // The columns of U are normalized. The cost is O(M*N) flops.
            TEMP1 = sqrt(FLOAT(M)) * EPSLN;
            for (p = 1; p <= NR; p++) { // 1973
               XSC = ONE / SNRM2( M, U(1,p), 1 );
               if ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) sscal( M, XSC, U(1,p), 1 );
            } // 1973

            // If the initial QRF is computed with row pivoting, the left
            // singular vectors must be adjusted.

            if (ROWPIV) slaswp( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 );

         } else {

         // .. the initial matrix A has almost orthogonal columns and
         // the second QRF is not needed

            slacpy('Upper', N, N, A, LDA, WORK(N+1), N );
            if ( L2PERT ) {
               XSC = sqrt(SMALL);
               for (p = 2; p <= N; p++) { // 5970
                  TEMP1 = XSC * WORK( N + (p-1)*N + p );
                  for (q = 1; q <= p - 1; q++) { // 5971
                     WORK(N+(q-1)*N+p)=-sign(TEMP1,WORK(N+(p-1)*N+q));
                  } // 5971
               } // 5970
            } else {
               slaset('Lower',N-1,N-1,ZERO,ZERO,WORK(N+2),N );
            }

            sgesvj('Upper', 'U', 'N', N, N, WORK(N+1), N, SVA, N, U, LDU, WORK(N+N*N+1), LWORK-N-N*N, INFO );

            SCALEM  = WORK(N+N*N+1);
            NUMRANK = NINT(WORK(N+N*N+2));
            for (p = 1; p <= N; p++) { // 6970
               scopy(N, WORK(N+(p-1)*N+1), 1, U(1,p), 1 );
               sscal(N, SVA(p), WORK(N+(p-1)*N+1), 1 );
            } // 6970

            strsm('Left', 'Upper', 'NoTrans', 'No UD', N, N, ONE, A, LDA, WORK(N+1), N );
            for (p = 1; p <= N; p++) { // 6972
               scopy(N, WORK(N+p), N, V(IWORK(p),1), LDV );
            } // 6972
            TEMP1 = sqrt(FLOAT(N))*EPSLN;
            for (p = 1; p <= N; p++) { // 6971
               XSC = ONE / SNRM2( N, V(1,p), 1 );
               if ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) sscal( N, XSC, V(1,p), 1 );
            } // 6971

            // Assemble the left singular vector matrix U (M x N).

            if ( N < M ) {
               slaset('A',  M-N, N, ZERO, ZERO, U(N+1,1), LDU );
               if ( N < N1 ) {
                  slaset('A',N,  N1-N, ZERO, ZERO,  U(1,N+1),LDU );
                  slaset('A',M-N,N1-N, ZERO, ONE,U(N+1,N+1),LDU );
               }
            }
            sormqr('Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );
            TEMP1 = sqrt(FLOAT(M))*EPSLN;
            for (p = 1; p <= N1; p++) { // 6973
               XSC = ONE / SNRM2( M, U(1,p), 1 );
               if ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) sscal( M, XSC, U(1,p), 1 );
            } // 6973

            if (ROWPIV) slaswp( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 );

         }

         // end of the  >> almost orthogonal case <<  in the full SVD

         } else {

         // This branch deploys a preconditioned Jacobi SVD with explicitly
         // accumulated rotations. It is included as optional, mainly for
         // experimental purposes. It does perform well, and can also be used.
         // In this implementation, this branch will be automatically activated
         // if the  condition number sigma_max(A) / sigma_min(A) is predicted
         // to be greater than the overflow threshold. This is because the
         // a posteriori computation of the singular vectors assumes robust
         // implementation of BLAS and some LAPACK procedures, capable of working
         // in presence of extreme values. Since that is not always the case, ...

         for (p = 1; p <= NR; p++) { // 7968
            scopy(N-p+1, A(p,p), LDA, V(p,p), 1 );
         } // 7968

         if ( L2PERT ) {
            XSC = sqrt(SMALL/EPSLN);
            for (q = 1; q <= NR; q++) { // 5969
               TEMP1 = XSC*( V(q,q) ).abs();
               for (p = 1; p <= N; p++) { // 5968
                  if ( ( p > q ) && ( (V(p,q)).abs() <= TEMP1 ) || ( p < q ) ) V(p,q) = sign( TEMP1, V(p,q) );
                  if (p < q) V(p,q) = - V(p,q);
               } // 5968
            } // 5969
         } else {
            slaset('U', NR-1, NR-1, ZERO, ZERO, V(1,2), LDV );
         }
          sgeqrf(N, NR, V, LDV, WORK(N+1), WORK(2*N+1), LWORK-2*N, IERR );
         slacpy('L', N, NR, V, LDV, WORK(2*N+1), N );

         for (p = 1; p <= NR; p++) { // 7969
            scopy(NR-p+1, V(p,p), LDV, U(p,p), 1 );
         } // 7969

         if ( L2PERT ) {
            XSC = sqrt(SMALL/EPSLN);
            for (q = 2; q <= NR; q++) { // 9970
               for (p = 1; p <= q - 1; p++) { // 9971
                  TEMP1 = XSC * min((U(p,p)).abs(),(U(q,q) ).abs() );
                  U[p][q] = - sign( TEMP1, U(q,p) );
               } // 9971
            } // 9970
         } else {
            slaset('U', NR-1, NR-1, ZERO, ZERO, U(1,2), LDU );
         }
          sgesvj('L', 'U', 'V', NR, NR, U, LDU, SVA, N, V, LDV, WORK(2*N+N*NR+1), LWORK-2*N-N*NR, INFO );
         SCALEM  = WORK(2*N+N*NR+1);
         NUMRANK = NINT(WORK(2*N+N*NR+2));

         if ( NR < N ) {
            slaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV );
            slaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV );
            slaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV );
         }
          sormqr('L','N',N,N,NR,WORK(2*N+1),N,WORK(N+1), V,LDV,WORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR );

            // Permute the rows of V using the (column) permutation from the
            // first QRF. Also, scale the columns to make them unit in
            // Euclidean norm. This applies to all cases.

            TEMP1 = sqrt(FLOAT(N)) * EPSLN;
            for (q = 1; q <= N; q++) { // 7972
               for (p = 1; p <= N; p++) { // 8972
                  WORK[2*N+N*NR+NR+IWORK(p)] = V(p,q);
               } // 8972
               for (p = 1; p <= N; p++) { // 8973
                  V[p][q] = WORK(2*N+N*NR+NR+p);
               } // 8973
               XSC = ONE / SNRM2( N, V(1,q), 1 );
               if ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) sscal( N, XSC, V(1,q), 1 );
            } // 7972

            // At this moment, V contains the right singular vectors of A.
            // Next, assemble the left singular vector matrix U (M x N).

         if ( NR < M ) {
            slaset('A',  M-NR, NR, ZERO, ZERO, U(NR+1,1), LDU );
            if ( NR < N1 ) {
               slaset('A',NR,  N1-NR, ZERO, ZERO,  U(1,NR+1),LDU );
               slaset('A',M-NR,N1-NR, ZERO, ONE,U(NR+1,NR+1),LDU );
            }
         }

         sormqr('Left', 'No Tr', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );

            if (ROWPIV) slaswp( N1, U, LDU, 1, M-1, IWORK(2*N+1), -1 );


         }
         if ( TRANSP ) {
            // .. swap U and V because the procedure worked on A^t
            for (p = 1; p <= N; p++) { // 6974
               sswap(N, U(1,p), 1, V(1,p), 1 );
            } // 6974
         }

      }
      // end of the full SVD

      // Undo scaling, if necessary (and possible)

      if ( USCAL2 <= (BIG/SVA(1))*USCAL1 ) {
         slascl('G', 0, 0, USCAL1, USCAL2, NR, 1, SVA, N, IERR );
         USCAL1 = ONE;
         USCAL2 = ONE;
      }

      if ( NR < N ) {
         for (p = NR+1; p <= N; p++) { // 3004
            SVA[p] = ZERO;
         } // 3004
      }

      WORK[1] = USCAL2 * SCALEM;
      WORK[2] = USCAL1;
      if (ERREST) WORK(3) = SCONDA;
      if ( LSVEC && RSVEC ) {
         WORK[4] = CONDR1;
         WORK[5] = CONDR2;
      }
      if ( L2TRAN ) {
         WORK[6] = ENTRA;
         WORK[7] = ENTRAT;
      }

      IWORK[1] = NR;
      IWORK[2] = NUMRANK;
      IWORK[3] = WARNING;

      return;
      // ..
      // .. END OF SGEJSV
      // ..
      }
