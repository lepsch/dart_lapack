      SUBROUTINE DGESVDQ( JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, INFO )
      // .. Scalar Arguments ..
      IMPLICIT    NONE
      String      JOBA, JOBP, JOBR, JOBU, JOBV;
      int         M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LWORK, LRWORK, INFO;
      // ..
      // .. Array Arguments ..
      double           A( LDA, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      double           S( * ), RWORK( * );
      int              IWORK( * );

*  =====================================================================

      // .. Parameters ..
      double           ZERO,         ONE;
      const          ZERO = 0.0D0, ONE = 1.0D0 ;
      // .. Local Scalars ..
      int         IERR, IWOFF, NR, N1, OPTRATIO, p, q;
      int         LWCON, LWQP3, LWRK_DGELQF, LWRK_DGESVD, LWRK_DGESVD2, LWRK_DGEQP3,  LWRK_DGEQRF, LWRK_DORMLQ, LWRK_DORMQR, LWRK_DORMQR2, LWLQF, LWQRF, LWSVD, LWSVD2, LWORQ, LWORQ2, LWORLQ, MINWRK, MINWRK2, OPTWRK, OPTWRK2, IMINWRK, RMINWRK;
      bool        ACCLA,  ACCLM, ACCLH, ASCALED, CONDA, DNTWU,  DNTWV, LQUERY, LSVC0, LSVEC, ROWPRM,  RSVEC, RTRANS, WNTUA, WNTUF,  WNTUR, WNTUS, WNTVA,   WNTVR;
      double           BIG, EPSLN, RTMP, SCONDA, SFMIN;
      // .. Local Arrays
      double           RDUMMY(1);
      // ..
      // .. External Subroutines (BLAS, LAPACK)
      // EXTERNAL DGELQF, DGEQP3, DGEQRF, DGESVD, DLACPY, DLAPMT, DLASCL, DLASET, DLASWP, DSCAL,  DPOCON, DORMLQ, DORMQR, XERBLA
      // ..
      // .. External Functions (BLAS, LAPACK)
      bool       LSAME;
      int        IDAMAX;
      double            DLANGE, DNRM2, DLAMCH;
      // EXTERNAL DLANGE, LSAME, IDAMAX, DNRM2, DLAMCH
      // ..
      // .. Intrinsic Functions ..

      // INTRINSIC ABS, MAX, MIN, DBLE, SQRT

      // Test the input arguments

      WNTUS  = LSAME( JOBU, 'S' ) .OR. LSAME( JOBU, 'U' )
      WNTUR  = LSAME( JOBU, 'R' )
      WNTUA  = LSAME( JOBU, 'A' )
      WNTUF  = LSAME( JOBU, 'F' )
      LSVC0  = WNTUS .OR. WNTUR .OR. WNTUA
      LSVEC  = LSVC0 .OR. WNTUF
      DNTWU  = LSAME( JOBU, 'N' )

      WNTVR  = LSAME( JOBV, 'R' )
      WNTVA  = LSAME( JOBV, 'A' ) .OR. LSAME( JOBV, 'V' )
      RSVEC  = WNTVR .OR. WNTVA
      DNTWV  = LSAME( JOBV, 'N' )

      ACCLA  = LSAME( JOBA, 'A' )
      ACCLM  = LSAME( JOBA, 'M' )
      CONDA  = LSAME( JOBA, 'E' )
      ACCLH  = LSAME( JOBA, 'H' ) .OR. CONDA

      ROWPRM = LSAME( JOBP, 'P' )
      RTRANS = LSAME( JOBR, 'T' )

      if ( ROWPRM ) {
         if ( CONDA ) {
            IMINWRK = MAX( 1, N + M - 1 + N )
         } else {
            IMINWRK = MAX( 1, N + M - 1 )
         }
         RMINWRK = MAX( 2, M )
      } else {
         if ( CONDA ) {
            IMINWRK = MAX( 1, N + N )
         } else {
            IMINWRK = MAX( 1, N )
         }
         RMINWRK = 2
      }
      LQUERY = (LIWORK .EQ. -1 .OR. LWORK .EQ. -1 .OR. LRWORK .EQ. -1)
      INFO  = 0
      if ( .NOT. ( ACCLA .OR. ACCLM .OR. ACCLH ) ) {
         INFO = -1
      } else if ( .NOT.( ROWPRM .OR. LSAME( JOBP, 'N' ) ) ) {
          INFO = -2
      } else if ( .NOT.( RTRANS .OR. LSAME( JOBR, 'N' ) ) ) {
          INFO = -3
      } else if ( .NOT.( LSVEC .OR. DNTWU ) ) {
         INFO = -4
      } else if ( WNTUR .AND. WNTVA ) {
         INFO = -5
      } else if ( .NOT.( RSVEC .OR. DNTWV )) {
         INFO = -5
      } else if ( M.LT.0 ) {
         INFO = -6
      } else if ( ( N.LT.0 ) .OR. ( N.GT.M ) ) {
         INFO = -7
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -9
      } else if ( LDU.LT.1 .OR. ( LSVC0 .AND. LDU.LT.M ) .OR. ( WNTUF .AND. LDU.LT.N ) ) {
         INFO = -12
      } else if ( LDV.LT.1 .OR. ( RSVEC .AND. LDV.LT.N ) .OR. ( CONDA .AND. LDV.LT.N ) ) {
         INFO = -14
      } else if ( LIWORK .LT. IMINWRK .AND. .NOT. LQUERY ) {
         INFO = -17
      }


      if ( INFO .EQ. 0 ) {
         // .. compute the minimal and the optimal workspace lengths
         // [[The expressions for computing the minimal and the optimal
         // values of LWORK are written with a lot of redundancy and
         // can be simplified. However, this detailed form is easier for
         // maintenance and modifications of the code.]]

         // .. minimal workspace length for DGEQP3 of an M x N matrix
         LWQP3 = 3 * N + 1
         // .. minimal workspace length for DORMQR to build left singular vectors
         if ( WNTUS .OR. WNTUR ) {
             LWORQ  = MAX( N  , 1 )
         } else if ( WNTUA ) {
             LWORQ = MAX( M , 1 )
         }
         // .. minimal workspace length for DPOCON of an N x N matrix
         LWCON = 3 * N
         // .. DGESVD of an N x N matrix
         LWSVD = MAX( 5 * N, 1 )
         if ( LQUERY ) {
             dgeqp3(M, N, A, LDA, IWORK, RDUMMY, RDUMMY, -1, IERR );
             LWRK_DGEQP3 = INT( RDUMMY(1) )
             if ( WNTUS .OR. WNTUR ) {
                 dormqr('L', 'N', M, N, N, A, LDA, RDUMMY, U, LDU, RDUMMY, -1, IERR );
                 LWRK_DORMQR = INT( RDUMMY(1) )
             } else if ( WNTUA ) {
                 dormqr('L', 'N', M, M, N, A, LDA, RDUMMY, U, LDU, RDUMMY, -1, IERR );
                 LWRK_DORMQR = INT( RDUMMY(1) )
             } else {
                 LWRK_DORMQR = 0
             }
         }
         MINWRK = 2
         OPTWRK = 2
         if ( .NOT. (LSVEC .OR. RSVEC )) {
             // .. minimal and optimal sizes of the workspace if
             // only the singular values are requested
             if ( CONDA ) {
                MINWRK = MAX( N+LWQP3, LWCON, LWSVD )
             } else {
                MINWRK = MAX( N+LWQP3, LWSVD )
             }
             if ( LQUERY ) {
                 dgesvd('N', 'N', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                 LWRK_DGESVD = INT( RDUMMY(1) )
                 if ( CONDA ) {
                    OPTWRK = MAX( N+LWRK_DGEQP3, N+LWCON, LWRK_DGESVD )
                 } else {
                    OPTWRK = MAX( N+LWRK_DGEQP3, LWRK_DGESVD )
                 }
             }
         } else if ( LSVEC .AND. (.NOT.RSVEC) ) {
             // .. minimal and optimal sizes of the workspace if the
             // singular values and the left singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD, LWORQ )
             } else {
                 MINWRK = N + MAX( LWQP3, LWSVD, LWORQ )
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   dgesvd('N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                } else {
                   dgesvd('O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                }
                LWRK_DGESVD = INT( RDUMMY(1) )
                if ( CONDA ) {
                    OPTWRK = N + MAX( LWRK_DGEQP3, LWCON, LWRK_DGESVD, LWRK_DORMQR )
                } else {
                    OPTWRK = N + MAX( LWRK_DGEQP3, LWRK_DGESVD, LWRK_DORMQR )
                }
             }
         } else if ( RSVEC .AND. (.NOT.LSVEC) ) {
             // .. minimal and optimal sizes of the workspace if the
             // singular values and the right singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD )
             } else {
                 MINWRK = N + MAX( LWQP3, LWSVD )
             }
             if ( LQUERY ) {
                 if ( RTRANS ) {
                     dgesvd('O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                 } else {
                     dgesvd('N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                 }
                 LWRK_DGESVD = INT( RDUMMY(1) )
                 if ( CONDA ) {
                     OPTWRK = N + MAX( LWRK_DGEQP3, LWCON, LWRK_DGESVD )
                 } else {
                     OPTWRK = N + MAX( LWRK_DGEQP3, LWRK_DGESVD )
                 }
             }
         } else {
             // .. minimal and optimal sizes of the workspace if the
             // full SVD is requested
             if ( RTRANS ) {
                 MINWRK = MAX( LWQP3, LWSVD, LWORQ )
                 if (CONDA) MINWRK = MAX( MINWRK, LWCON );
                 MINWRK = MINWRK + N
                 if ( WNTVA ) {
                    // .. minimal workspace length for N x N/2 DGEQRF
                    LWQRF  = MAX( N/2, 1 )
                    // .. minimal workspace length for N/2 x N/2 DGESVD
                    LWSVD2 = MAX( 5 * (N/2), 1 )
                    LWORQ2 = MAX( N, 1 )
                    MINWRK2 = MAX( LWQP3, N/2+LWQRF, N/2+LWSVD2, N/2+LWORQ2, LWORQ )
                    if (CONDA) MINWRK2 = MAX( MINWRK2, LWCON );
                    MINWRK2 = N + MINWRK2
                    MINWRK = MAX( MINWRK, MINWRK2 )
                 }
             } else {
                 MINWRK = MAX( LWQP3, LWSVD, LWORQ )
                 if (CONDA) MINWRK = MAX( MINWRK, LWCON );
                 MINWRK = MINWRK + N
                 if ( WNTVA ) {
                    // .. minimal workspace length for N/2 x N DGELQF
                    LWLQF  = MAX( N/2, 1 )
                    LWSVD2 = MAX( 5 * (N/2), 1 )
                    LWORLQ = MAX( N , 1 )
                    MINWRK2 = MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWORLQ, LWORQ )
                    if (CONDA) MINWRK2 = MAX( MINWRK2, LWCON );
                    MINWRK2 = N + MINWRK2
                    MINWRK = MAX( MINWRK, MINWRK2 )
                 }
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   dgesvd('O', 'A', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                   LWRK_DGESVD = INT( RDUMMY(1) )
                   OPTWRK = MAX(LWRK_DGEQP3,LWRK_DGESVD,LWRK_DORMQR)
                   if (CONDA) OPTWRK = MAX( OPTWRK, LWCON );
                   OPTWRK = N + OPTWRK
                   if ( WNTVA ) {
                       dgeqrf(N,N/2,U,LDU,RDUMMY,RDUMMY,-1,IERR);
                       LWRK_DGEQRF = INT( RDUMMY(1) )
                       dgesvd('S', 'O', N/2,N/2, V,LDV, S, U,LDU, V, LDV, RDUMMY, -1, IERR );
                       LWRK_DGESVD2 = INT( RDUMMY(1) )
                       dormqr('R', 'C', N, N, N/2, U, LDU, RDUMMY, V, LDV, RDUMMY, -1, IERR );
                       LWRK_DORMQR2 = INT( RDUMMY(1) )
                       OPTWRK2 = MAX( LWRK_DGEQP3, N/2+LWRK_DGEQRF, N/2+LWRK_DGESVD2, N/2+LWRK_DORMQR2 )
                       if (CONDA) OPTWRK2 = MAX( OPTWRK2, LWCON );
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   }
                } else {
                   dgesvd('S', 'O', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                   LWRK_DGESVD = INT( RDUMMY(1) )
                   OPTWRK = MAX(LWRK_DGEQP3,LWRK_DGESVD,LWRK_DORMQR)
                   if (CONDA) OPTWRK = MAX( OPTWRK, LWCON );
                   OPTWRK = N + OPTWRK
                   if ( WNTVA ) {
                      dgelqf(N/2,N,U,LDU,RDUMMY,RDUMMY,-1,IERR);
                      LWRK_DGELQF = INT( RDUMMY(1) )
                      dgesvd('S','O', N/2,N/2, V, LDV, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                      LWRK_DGESVD2 = INT( RDUMMY(1) )
                      dormlq('R', 'N', N, N, N/2, U, LDU, RDUMMY, V, LDV, RDUMMY,-1,IERR );
                      LWRK_DORMLQ = INT( RDUMMY(1) )
                      OPTWRK2 = MAX( LWRK_DGEQP3, N/2+LWRK_DGELQF, N/2+LWRK_DGESVD2, N/2+LWRK_DORMLQ )
                       if (CONDA) OPTWRK2 = MAX( OPTWRK2, LWCON );
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   }
                }
             }
         }

         MINWRK = MAX( 2, MINWRK )
         OPTWRK = MAX( 2, OPTWRK )
         IF ( LWORK .LT. MINWRK .AND. (.NOT.LQUERY) ) INFO = -19

      }

      if (INFO .EQ. 0 .AND. LRWORK .LT. RMINWRK .AND. .NOT. LQUERY) {
         INFO = -21
      }
      if ( INFO.NE.0 ) {
         xerbla('DGESVDQ', -INFO );
         RETURN
      } else if ( LQUERY ) {

      // Return optimal workspace

          IWORK(1) = IMINWRK
          WORK(1) = OPTWRK
          WORK(2) = MINWRK
          RWORK(1) = RMINWRK
          RETURN
      }

      // Quick return if the matrix is void.

      if ( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) {
      // .. all output is void.
         RETURN
      }

      BIG = DLAMCH('O')
      ASCALED = false;
      IWOFF = 1
      if ( ROWPRM ) {
            IWOFF = M
            // .. reordering the rows in decreasing sequence in the
            // ell-infinity norm - this enhances numerical robustness in
            // the case of differently scaled rows.
            for (p = 1; p <= M; p++) { // 1904
                // RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
                // [[DLANGE will return NaN if an entry of the p-th row is Nan]]
                RWORK(p) = DLANGE( 'M', 1, N, A(p,1), LDA, RDUMMY )
                // .. check for NaN's and Inf's
                if ( ( RWORK(p) .NE. RWORK(p) ) .OR. ( (RWORK(p)*ZERO) .NE. ZERO ) ) {
                    INFO = -8
                    xerbla('DGESVDQ', -INFO );
                    RETURN
                }
            } // 1904
            for (p = 1; p <= M - 1; p++) { // 1952
            q = IDAMAX( M-p+1, RWORK(p), 1 ) + p - 1
            IWORK(N+p) = q
            if ( p .NE. q ) {
               RTMP     = RWORK(p)
               RWORK(p) = RWORK(q)
               RWORK(q) = RTMP
            }
            } // 1952

            if ( RWORK(1) .EQ. ZERO ) {
               // Quick return: A is the M x N zero matrix.
               NUMRANK = 0
               dlaset('G', N, 1, ZERO, ZERO, S, N );
               if (WNTUS) CALL DLASET('G', M, N, ZERO, ONE, U, LDU);
               if (WNTUA) CALL DLASET('G', M, M, ZERO, ONE, U, LDU);
               if (WNTVA) CALL DLASET('G', N, N, ZERO, ONE, V, LDV);
               if ( WNTUF ) {
                   dlaset('G', N, 1, ZERO, ZERO, WORK, N );
                   dlaset('G', M, N, ZERO,  ONE, U, LDU );
               }
               for (p = 1; p <= N; p++) { // 5001
                   IWORK(p) = p
               } // 5001
               if ( ROWPRM ) {
                   for (p = N + 1; p <= N + M - 1; p++) { // 5002
                       IWORK(p) = p - N
                   } // 5002
               }
               if (CONDA) RWORK(1) = -1;
               RWORK(2) = -1
               RETURN
            }

            if ( RWORK(1) .GT. BIG / SQRT(DBLE(M)) ) {
                // .. to prevent overflow in the QR factorization, scale the
                // matrix by 1/sqrt(M) if too large entry detected
                dlascl('G',0,0,SQRT(DBLE(M)),ONE, M,N, A,LDA, IERR);
                ASCALED = true;
            }
            dlaswp(N, A, LDA, 1, M-1, IWORK(N+1), 1 );
      }

*    .. At this stage, preemptive scaling is done only to avoid column
*    norms overflows during the QR factorization. The SVD procedure should
*    have its own scaling to save the singular values from overflows and
*    underflows. That depends on the SVD procedure.

      if ( .NOT.ROWPRM ) {
          RTMP = DLANGE( 'M', M, N, A, LDA, RDUMMY )
          if ( ( RTMP .NE. RTMP ) .OR. ( (RTMP*ZERO) .NE. ZERO ) ) {
               INFO = -8
               xerbla('DGESVDQ', -INFO );
               RETURN
          }
          if ( RTMP .GT. BIG / SQRT(DBLE(M)) ) {
              // .. to prevent overflow in the QR factorization, scale the
              // matrix by 1/sqrt(M) if too large entry detected
              dlascl('G',0,0, SQRT(DBLE(M)),ONE, M,N, A,LDA, IERR);
              ASCALED = true;
          }
      }

      // .. QR factorization with column pivoting

      // A * P = Q * [ R ]
                  // [ 0 ]

      for (p = 1; p <= N; p++) { // 1963
         // .. all columns are free columns
         IWORK(p) = 0
      } // 1963
      dgeqp3(M, N, A, LDA, IWORK, WORK, WORK(N+1), LWORK-N, IERR );

*    If the user requested accuracy level allows truncation in the
*    computed upper triangular factor, the matrix R is examined and,
*    if possible, replaced with its leading upper trapezoidal part.

      EPSLN = DLAMCH('E')
      SFMIN = DLAMCH('S')
      // SMALL = SFMIN / EPSLN
      NR = N

      if ( ACCLA ) {

         // Standard absolute error bound suffices. All sigma_i with
         // sigma_i < N*EPS*||A||_F are flushed to zero. This is an
         // aggressive enforcement of lower numerical rank by introducing a
         // backward error of the order of N*EPS*||A||_F.
         NR = 1
         RTMP = SQRT(DBLE(N))*EPSLN
         for (p = 2; p <= N; p++) { // 3001
            IF ( ABS(A(p,p)) .LT. (RTMP*ABS(A(1,1))) ) GO TO 3002
               NR = NR + 1
         } // 3001
         } // 3002

      } else if ( ACCLM ) {
         // .. similarly as above, only slightly more gentle (less aggressive).
         // Sudden drop on the diagonal of R is used as the criterion for being
         // close-to-rank-deficient. The threshold is set to EPSLN=DLAMCH('E').
         // [[This can be made more flexible by replacing this hard-coded value
         // with a user specified threshold.]] Also, the values that underflow
         // will be truncated.
         NR = 1
         for (p = 2; p <= N; p++) { // 3401
            IF ( ( ABS(A(p,p)) .LT. (EPSLN*ABS(A(p-1,p-1))) ) .OR. ( ABS(A(p,p)) .LT. SFMIN ) ) GO TO 3402
            NR = NR + 1
         } // 3401
         } // 3402

      } else {
         // .. RRQR not authorized to determine numerical rank except in the
         // obvious case of zero pivots.
         // .. inspect R for exact zeros on the diagonal;
         // R(i,i)=0 => R(i:N,i:N)=0.
         NR = 1
         for (p = 2; p <= N; p++) { // 3501
            IF ( ABS(A(p,p)) .EQ. ZERO ) GO TO 3502
            NR = NR + 1
         } // 3501
         } // 3502

         if ( CONDA ) {
            // Estimate the scaled condition number of A. Use the fact that it is
            // the same as the scaled condition number of R.
               // .. V is used as workspace
               dlacpy('U', N, N, A, LDA, V, LDV );
               // Only the leading NR x NR submatrix of the triangular factor
               // is considered. Only if NR=N will this give a reliable error
               // bound. However, even for NR < N, this can be used on an
               // expert level and obtain useful information in the sense of
               // perturbation theory.
               for (p = 1; p <= NR; p++) { // 3053
                  RTMP = DNRM2( p, V(1,p), 1 )
                  dscal(p, ONE/RTMP, V(1,p), 1 );
               } // 3053
               if ( .NOT. ( LSVEC .OR. RSVEC ) ) {
                   dpocon('U', NR, V, LDV, ONE, RTMP, WORK, IWORK(N+IWOFF), IERR );
               } else {
                   dpocon('U', NR, V, LDV, ONE, RTMP, WORK(N+1), IWORK(N+IWOFF), IERR );
               }
               SCONDA = ONE / SQRT(RTMP)
            // For NR=N, SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1),
            // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
            // See the reference [1] for more details.
         }

      }

      if ( WNTUR ) {
          N1 = NR
      } else if ( WNTUS .OR. WNTUF) {
          N1 = N
      } else if ( WNTUA ) {
          N1 = M
      }

      if ( .NOT. ( RSVEC .OR. LSVEC ) ) {
*.......................................................................
         // .. only the singular values are requested
*.......................................................................
         if ( RTRANS ) {

          // .. compute the singular values of R**T = [A](1:NR,1:N)**T
            // .. set the lower triangle of [A] to [A](1:NR,1:N)**T and
            // the upper triangle of [A] to zero.
            DO 1146 p = 1, MIN( N, NR )
               for (q = p + 1; q <= N; q++) { // 1147
                  A(q,p) = A(p,q)
                  if (q .LE. NR) A(p,q) = ZERO;
               } // 1147
            } // 1146

            dgesvd('N', 'N', N, NR, A, LDA, S, U, LDU, V, LDV, WORK, LWORK, INFO );

         } else {

            // .. compute the singular values of R = [A](1:NR,1:N)

            if (NR .GT. 1) CALL DLASET( 'L', NR-1,NR-1, ZERO,ZERO, A(2,1), LDA );
            dgesvd('N', 'N', NR, N, A, LDA, S, U, LDU, V, LDV, WORK, LWORK, INFO );

         }

      } else if ( LSVEC .AND. ( .NOT. RSVEC) ) {
*.......................................................................
        // .. the singular values and the left singular vectors requested
*.......................................................................""""""""
         if ( RTRANS ) {
             // .. apply DGESVD to R**T
             // .. copy R**T into [U] and overwrite [U] with the right singular
             // vectors of R
            for (p = 1; p <= NR; p++) { // 1192
               for (q = p; q <= N; q++) { // 1193
                  U(q,p) = A(p,q)
               } // 1193
            } // 1192
            if (NR .GT. 1) CALL DLASET( 'U', NR-1,NR-1, ZERO,ZERO, U(1,2), LDU );
            // .. the left singular vectors not computed, the NR right singular
            // vectors overwrite [U](1:NR,1:NR) as transposed. These
            // will be pre-multiplied by Q to build the left singular vectors of A.
               dgesvd('N', 'O', N, NR, U, LDU, S, U, LDU, U, LDU, WORK(N+1), LWORK-N, INFO );

               for (p = 1; p <= NR; p++) { // 1119
                   for (q = p + 1; q <= NR; q++) { // 1120
                      RTMP   = U(q,p)
                      U(q,p) = U(p,q)
                      U(p,q) = RTMP
                   } // 1120
               } // 1119

         } else {
             // .. apply DGESVD to R
             // .. copy R into [U] and overwrite [U] with the left singular vectors
             dlacpy('U', NR, N, A, LDA, U, LDU );
             if (NR .GT. 1) CALL DLASET( 'L', NR-1, NR-1, ZERO, ZERO, U(2,1), LDU );
             // .. the right singular vectors not computed, the NR left singular
             // vectors overwrite [U](1:NR,1:NR)
                dgesvd('O', 'N', NR, N, U, LDU, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                // .. now [U](1:NR,1:NR) contains the NR left singular vectors of
                // R. These will be pre-multiplied by Q to build the left singular
                // vectors of A.
         }

            // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
         if ( ( NR .LT. M ) .AND. ( .NOT.WNTUF ) ) {
             dlaset('A', M-NR, NR, ZERO, ZERO, U(NR+1,1), LDU);
             if ( NR .LT. N1 ) {
                dlaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1), LDU );
                dlaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
             }
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         if (.NOT.WNTUF) CALL DORMQR( 'L', 'N', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );
         if (ROWPRM .AND. .NOT.WNTUF) CALL DLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 );

      } else if ( RSVEC .AND. ( .NOT. LSVEC ) ) {
*.......................................................................
        // .. the singular values and the right singular vectors requested
*.......................................................................
          if ( RTRANS ) {
             // .. apply DGESVD to R**T
             // .. copy R**T into V and overwrite V with the left singular vectors
            for (p = 1; p <= NR; p++) { // 1165
               for (q = p; q <= N; q++) { // 1166
                  V(q,p) = (A(p,q))
               } // 1166
            } // 1165
            if (NR .GT. 1) CALL DLASET( 'U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV );
            // .. the left singular vectors of R**T overwrite V, the right singular
            // vectors not computed
            if ( WNTVR .OR. ( NR .EQ. N ) ) {
               dgesvd('O', 'N', N, NR, V, LDV, S, U, LDU, U, LDU, WORK(N+1), LWORK-N, INFO );

               for (p = 1; p <= NR; p++) { // 1121
                   for (q = p + 1; q <= NR; q++) { // 1122
                      RTMP   = V(q,p)
                      V(q,p) = V(p,q)
                      V(p,q) = RTMP
                   } // 1122
               } // 1121

               if ( NR .LT. N ) {
                   for (p = 1; p <= NR; p++) { // 1103
                      for (q = NR + 1; q <= N; q++) { // 1104
                          V(p,q) = V(q,p)
                      } // 1104
                   } // 1103
               }
               dlapmt( false , NR, N, V, LDV, IWORK );
            } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:N,1:NR)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the QR factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                dlaset('G', N, N-NR, ZERO, ZERO, V(1,NR+1), LDV);
                dgesvd('O', 'N', N, N, V, LDV, S, U, LDU, U, LDU, WORK(N+1), LWORK-N, INFO );

                for (p = 1; p <= N; p++) { // 1123
                   for (q = p + 1; q <= N; q++) { // 1124
                      RTMP   = V(q,p)
                      V(q,p) = V(p,q)
                      V(p,q) = RTMP
                   } // 1124
                } // 1123
                dlapmt( false , N, N, V, LDV, IWORK );
            }

          } else {
             // .. aply DGESVD to R
             // .. copy R into V and overwrite V with the right singular vectors
             dlacpy('U', NR, N, A, LDA, V, LDV );
             if (NR .GT. 1) CALL DLASET( 'L', NR-1, NR-1, ZERO, ZERO, V(2,1), LDV );
             // .. the right singular vectors overwrite V, the NR left singular
             // vectors stored in U(1:NR,1:NR)
             if ( WNTVR .OR. ( NR .EQ. N ) ) {
                dgesvd('N', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                dlapmt( false , NR, N, V, LDV, IWORK );
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
             } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:NR,1:N)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the LQ factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                 dlaset('G', N-NR, N, ZERO,ZERO, V(NR+1,1), LDV);
                 dgesvd('N', 'O', N, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                 dlapmt( false , N, N, V, LDV, IWORK );
             }
             // .. now [V] contains the transposed matrix of the right singular
             // vectors of A.
          }

      } else {
*.......................................................................
        // .. FULL SVD requested
*.......................................................................
         if ( RTRANS ) {

             // .. apply DGESVD to R**T [[this option is left for R&D&T]]

            if ( WNTVR .OR. ( NR .EQ. N ) ) {
             // .. copy R**T into [V] and overwrite [V] with the left singular
             // vectors of R**T
            for (p = 1; p <= NR; p++) { // 1168
               for (q = p; q <= N; q++) { // 1169
                  V(q,p) = A(p,q)
               } // 1169
            } // 1168
            if (NR .GT. 1) CALL DLASET( 'U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV );

            // .. the left singular vectors of R**T overwrite [V], the NR right
            // singular vectors of R**T stored in [U](1:NR,1:NR) as transposed
               dgesvd('O', 'A', N, NR, V, LDV, S, V, LDV, U, LDU, WORK(N+1), LWORK-N, INFO );
               // .. assemble V
               for (p = 1; p <= NR; p++) { // 1115
                  for (q = p + 1; q <= NR; q++) { // 1116
                     RTMP   = V(q,p)
                     V(q,p) = V(p,q)
                     V(p,q) = RTMP
                  } // 1116
               } // 1115
               if ( NR .LT. N ) {
                   for (p = 1; p <= NR; p++) { // 1101
                      for (q = NR+1; q <= N; q++) { // 1102
                         V(p,q) = V(q,p)
                      } // 1102
                   } // 1101
               }
               dlapmt( false , NR, N, V, LDV, IWORK );

                for (p = 1; p <= NR; p++) { // 1117
                   for (q = p + 1; q <= NR; q++) { // 1118
                      RTMP   = U(q,p)
                      U(q,p) = U(p,q)
                      U(p,q) = RTMP
                   } // 1118
                } // 1117

                if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                  dlaset('A', M-NR,NR, ZERO,ZERO, U(NR+1,1), LDU);
                  if ( NR .LT. N1 ) {
                     dlaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     dlaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
                  }
               }

            } else {
                // .. need all N right singular vectors and NR < N
             // .. copy R**T into [V] and overwrite [V] with the left singular
             // vectors of R**T
                // [[The optimal ratio N/NR for using QRF instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'DGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = MAX( OPTRATIO, 2 )
                OPTRATIO = 2
                if ( OPTRATIO*NR .GT. N ) {
                   for (p = 1; p <= NR; p++) { // 1198
                      for (q = p; q <= N; q++) { // 1199
                         V(q,p) = A(p,q)
                      } // 1199
                   } // 1198
                   if (NR .GT. 1) CALL DLASET('U',NR-1,NR-1, ZERO,ZERO, V(1,2),LDV);

                   dlaset('A',N,N-NR,ZERO,ZERO,V(1,NR+1),LDV);
                   dgesvd('O', 'A', N, N, V, LDV, S, V, LDV, U, LDU, WORK(N+1), LWORK-N, INFO );

                   for (p = 1; p <= N; p++) { // 1113
                      for (q = p + 1; q <= N; q++) { // 1114
                         RTMP   = V(q,p)
                         V(q,p) = V(p,q)
                         V(p,q) = RTMP
                      } // 1114
                   } // 1113
                   dlapmt( false , N, N, V, LDV, IWORK );
               // .. assemble the left singular vector matrix U of dimensions
               // (M x N1), i.e. (M x N) or (M x M).

                   for (p = 1; p <= N; p++) { // 1111
                      for (q = p + 1; q <= N; q++) { // 1112
                         RTMP   = U(q,p)
                         U(q,p) = U(p,q)
                         U(p,q) = RTMP
                      } // 1112
                   } // 1111

                   if ( ( N .LT. M ) .AND. .NOT.(WNTUF)) {
                      dlaset('A',M-N,N,ZERO,ZERO,U(N+1,1),LDU);
                      if ( N .LT. N1 ) {
                        dlaset('A',N,N1-N,ZERO,ZERO,U(1,N+1),LDU);
                        dlaset('A',M-N,N1-N,ZERO,ONE, U(N+1,N+1), LDU );
                      }
                   }
                } else {
                   // .. copy R**T into [U] and overwrite [U] with the right
                   // singular vectors of R
                   for (p = 1; p <= NR; p++) { // 1196
                      for (q = p; q <= N; q++) { // 1197
                         U(q,NR+p) = A(p,q)
                      } // 1197
                   } // 1196
                   if (NR .GT. 1) CALL DLASET('U',NR-1,NR-1,ZERO,ZERO,U(1,NR+2),LDU);
                   dgeqrf(N, NR, U(1,NR+1), LDU, WORK(N+1), WORK(N+NR+1), LWORK-N-NR, IERR );
                   for (p = 1; p <= NR; p++) { // 1143
                       for (q = 1; q <= N; q++) { // 1144
                           V(q,p) = U(p,NR+q)
                       } // 1144
                   } // 1143
                  dlaset('U',NR-1,NR-1,ZERO,ZERO,V(1,2),LDV);
                  dgesvd('S', 'O', NR, NR, V, LDV, S, U, LDU, V,LDV, WORK(N+NR+1),LWORK-N-NR, INFO );
                  dlaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV);
                  dlaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV);
                  dlaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV);
                  dormqr('R','C', N, N, NR, U(1,NR+1), LDU, WORK(N+1),V,LDV,WORK(N+NR+1),LWORK-N-NR,IERR);
                  dlapmt( false , N, N, V, LDV, IWORK );
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x NR) or (M x N) or (M x M).
                  if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                     dlaset('A',M-NR,NR,ZERO,ZERO,U(NR+1,1),LDU);
                     if ( NR .LT. N1 ) {
                     dlaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     dlaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1),LDU);
                     }
                  }
                }
            }

         } else {

             // .. apply DGESVD to R [[this is the recommended option]]

             if ( WNTVR .OR. ( NR .EQ. N ) ) {
                 // .. copy R into [V] and overwrite V with the right singular vectors
                 dlacpy('U', NR, N, A, LDA, V, LDV );
                if (NR .GT. 1) CALL DLASET( 'L', NR-1,NR-1, ZERO,ZERO, V(2,1), LDV );
                // .. the right singular vectors of R overwrite [V], the NR left
                // singular vectors of R stored in [U](1:NR,1:NR)
                dgesvd('S', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                dlapmt( false , NR, N, V, LDV, IWORK );
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
               if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                  dlaset('A', M-NR,NR, ZERO,ZERO, U(NR+1,1), LDU);
                  if ( NR .LT. N1 ) {
                     dlaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     dlaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
                  }
               }

             } else {
               // .. need all N right singular vectors and NR < N
               // .. the requested number of the left singular vectors
                // is then N1 (N or M)
                // [[The optimal ratio N/NR for using LQ instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'DGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = MAX( OPTRATIO, 2 )
               OPTRATIO = 2
               if ( OPTRATIO * NR .GT. N ) {
                  dlacpy('U', NR, N, A, LDA, V, LDV );
                  if (NR .GT. 1) CALL DLASET('L', NR-1,NR-1, ZERO,ZERO, V(2,1),LDV);
               // .. the right singular vectors of R overwrite [V], the NR left
                  // singular vectors of R stored in [U](1:NR,1:NR)
                  dlaset('A', N-NR,N, ZERO,ZERO, V(NR+1,1),LDV);
                  dgesvd('S', 'O', N, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                  dlapmt( false , N, N, V, LDV, IWORK );
                  // .. now [V] contains the transposed matrix of the right
                  // singular vectors of A. The leading N left singular vectors
                  // are in [U](1:N,1:N)
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x N1), i.e. (M x N) or (M x M).
                  if ( ( N .LT. M ) .AND. .NOT.(WNTUF)) {
                      dlaset('A',M-N,N,ZERO,ZERO,U(N+1,1),LDU);
                      if ( N .LT. N1 ) {
                        dlaset('A',N,N1-N,ZERO,ZERO,U(1,N+1),LDU);
                        dlaset('A',M-N,N1-N,ZERO,ONE, U(N+1,N+1), LDU );
                      }
                  }
               } else {
                  dlacpy('U', NR, N, A, LDA, U(NR+1,1), LDU );
                  if (NR .GT. 1) CALL DLASET('L',NR-1,NR-1,ZERO,ZERO,U(NR+2,1),LDU);
                  dgelqf(NR, N, U(NR+1,1), LDU, WORK(N+1), WORK(N+NR+1), LWORK-N-NR, IERR );
                  dlacpy('L',NR,NR,U(NR+1,1),LDU,V,LDV);
                  if (NR .GT. 1) CALL DLASET('U',NR-1,NR-1,ZERO,ZERO,V(1,2),LDV);
                  dgesvd('S', 'O', NR, NR, V, LDV, S, U, LDU, V, LDV, WORK(N+NR+1), LWORK-N-NR, INFO );
                  dlaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV);
                  dlaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV);
                  dlaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV);
                  dormlq('R','N',N,N,NR,U(NR+1,1),LDU,WORK(N+1), V, LDV, WORK(N+NR+1),LWORK-N-NR,IERR);
                  dlapmt( false , N, N, V, LDV, IWORK );
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
                  if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                     dlaset('A',M-NR,NR,ZERO,ZERO,U(NR+1,1),LDU);
                     if ( NR .LT. N1 ) {
                     dlaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     dlaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
                     }
                  }
               }
             }
         // .. end of the "R**T or R" branch
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         if (.NOT. WNTUF) CALL DORMQR( 'L', 'N', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );
         if (ROWPRM .AND. .NOT.WNTUF) CALL DLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 );

      // ... end of the "full SVD" branch
      }

      // Check whether some singular values are returned as zeros, e.g.
      // due to underflow, and update the numerical rank.
      p = NR
      DO 4001 q = p, 1, -1
          IF ( S(q) .GT. ZERO ) GO TO 4002
          NR = NR - 1
      } // 4001
      } // 4002

      // .. if numerical rank deficiency is detected, the truncated
      // singular values are set to zero.
      if (NR .LT. N) CALL DLASET( 'G', N-NR,1, ZERO,ZERO, S(NR+1), N );
      // .. undo scaling; this may cause overflow in the largest singular
      // values.
      if (ASCALED) CALL DLASCL( 'G',0,0, ONE,SQRT(DBLE(M)), NR,1, S, N, IERR );
      if (CONDA) RWORK(1) = SCONDA;
      RWORK(2) = p - NR
      // .. p-NR is the number of singular values that are computed as
      // exact zeros in DGESVD() applied to the (possibly truncated)
      // full row rank triangular (trapezoidal) factor of A.
      NUMRANK = NR

      RETURN

      // End of DGESVDQ

      }
