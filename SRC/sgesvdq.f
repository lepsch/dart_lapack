      SUBROUTINE SGESVDQ( JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, INFO )
      // .. Scalar Arguments ..
      IMPLICIT    NONE
      String      JOBA, JOBP, JOBR, JOBU, JOBV;
      int         M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LWORK, LRWORK, INFO;
      // ..
      // .. Array Arguments ..
      REAL        A( LDA, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      REAL        S( * ), RWORK( * )
      int         IWORK( * );

*  =====================================================================

      // .. Parameters ..
      REAL        ZERO,         ONE
      const     ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      int         IERR, IWOFF, NR, N1, OPTRATIO, p, q;
      int         LWCON, LWQP3, LWRK_SGELQF, LWRK_SGESVD, LWRK_SGESVD2, LWRK_SGEQP3,  LWRK_SGEQRF, LWRK_SORMLQ, LWRK_SORMQR, LWRK_SORMQR2, LWLQF, LWQRF, LWSVD, LWSVD2, LWORQ, LWORQ2, LWUNLQ, MINWRK, MINWRK2, OPTWRK, OPTWRK2, IMINWRK, RMINWRK;
      bool        ACCLA,  ACCLM, ACCLH, ASCALED, CONDA, DNTWU,  DNTWV, LQUERY, LSVC0, LSVEC, ROWPRM,  RSVEC, RTRANS, WNTUA, WNTUF,  WNTUR, WNTUS, WNTVA,   WNTVR;
      REAL        BIG, EPSLN, RTMP, SCONDA, SFMIN
      // ..
      // .. Local Arrays
      REAL        RDUMMY(1)
      // ..
      // .. External Subroutines (BLAS, LAPACK)
      // EXTERNAL SGELQF, SGEQP3, SGEQRF, SGESVD, SLACPY, SLAPMT, SLASCL, SLASET, SLASWP, SSCAL,  SPOCON, SORMLQ, SORMQR, XERBLA
      // ..
      // .. External Functions (BLAS, LAPACK)
      bool       LSAME;
      int        ISAMAX;
      REAL        SLANGE, SNRM2, SLAMCH
      // EXTERNAL SLANGE, LSAME, ISAMAX, SNRM2, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      WNTUS  = LSAME( JOBU, 'S' ) || LSAME( JOBU, 'U' )
      WNTUR  = LSAME( JOBU, 'R' )
      WNTUA  = LSAME( JOBU, 'A' )
      WNTUF  = LSAME( JOBU, 'F' )
      LSVC0  = WNTUS || WNTUR || WNTUA
      LSVEC  = LSVC0 || WNTUF
      DNTWU  = LSAME( JOBU, 'N' )

      WNTVR  = LSAME( JOBV, 'R' )
      WNTVA  = LSAME( JOBV, 'A' ) || LSAME( JOBV, 'V' )
      RSVEC  = WNTVR || WNTVA
      DNTWV  = LSAME( JOBV, 'N' )

      ACCLA  = LSAME( JOBA, 'A' )
      ACCLM  = LSAME( JOBA, 'M' )
      CONDA  = LSAME( JOBA, 'E' )
      ACCLH  = LSAME( JOBA, 'H' ) || CONDA

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
      LQUERY = (LIWORK == -1 || LWORK == -1 || LRWORK == -1)
      INFO  = 0
      if ( .NOT. ( ACCLA || ACCLM || ACCLH ) ) {
         INFO = -1
      } else if ( .NOT.( ROWPRM || LSAME( JOBP, 'N' ) ) ) {
          INFO = -2
      } else if ( .NOT.( RTRANS || LSAME( JOBR, 'N' ) ) ) {
          INFO = -3
      } else if ( .NOT.( LSVEC || DNTWU ) ) {
         INFO = -4
      } else if ( WNTUR && WNTVA ) {
         INFO = -5
      } else if ( .NOT.( RSVEC || DNTWV )) {
         INFO = -5
      } else if ( M < 0 ) {
         INFO = -6
      } else if ( ( N < 0 ) || ( N.GT.M ) ) {
         INFO = -7
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -9
      } else if ( LDU < 1 || ( LSVC0 && LDU < M ) || ( WNTUF && LDU < N ) ) {
         INFO = -12
      } else if ( LDV < 1 || ( RSVEC && LDV < N ) || ( CONDA && LDV < N ) ) {
         INFO = -14
      } else if ( LIWORK < IMINWRK && .NOT. LQUERY ) {
         INFO = -17
      }


      if ( INFO == 0 ) {
         // .. compute the minimal and the optimal workspace lengths
         // [[The expressions for computing the minimal and the optimal
         // values of LWORK are written with a lot of redundancy and
         // can be simplified. However, this detailed form is easier for
         // maintenance and modifications of the code.]]

         // .. minimal workspace length for SGEQP3 of an M x N matrix
         LWQP3 = 3 * N + 1
         // .. minimal workspace length for SORMQR to build left singular vectors
         if ( WNTUS || WNTUR ) {
             LWORQ  = MAX( N  , 1 )
         } else if ( WNTUA ) {
             LWORQ = MAX( M , 1 )
         }
         // .. minimal workspace length for SPOCON of an N x N matrix
         LWCON = 3 * N
         // .. SGESVD of an N x N matrix
         LWSVD = MAX( 5 * N, 1 )
         if ( LQUERY ) {
             sgeqp3(M, N, A, LDA, IWORK, RDUMMY, RDUMMY, -1, IERR );
             LWRK_SGEQP3 = INT( RDUMMY(1) )
             if ( WNTUS || WNTUR ) {
                 sormqr('L', 'N', M, N, N, A, LDA, RDUMMY, U, LDU, RDUMMY, -1, IERR );
                 LWRK_SORMQR = INT( RDUMMY(1) )
             } else if ( WNTUA ) {
                 sormqr('L', 'N', M, M, N, A, LDA, RDUMMY, U, LDU, RDUMMY, -1, IERR );
                 LWRK_SORMQR = INT( RDUMMY(1) )
             } else {
                 LWRK_SORMQR = 0
             }
         }
         MINWRK = 2
         OPTWRK = 2
         if ( .NOT. (LSVEC || RSVEC )) {
             // .. minimal and optimal sizes of the workspace if
             // only the singular values are requested
             if ( CONDA ) {
                MINWRK = MAX( N+LWQP3, LWCON, LWSVD )
             } else {
                MINWRK = MAX( N+LWQP3, LWSVD )
             }
             if ( LQUERY ) {
                 sgesvd('N', 'N', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                 LWRK_SGESVD = INT( RDUMMY(1) )
                 if ( CONDA ) {
                    OPTWRK = MAX( N+LWRK_SGEQP3, N+LWCON, LWRK_SGESVD )
                 } else {
                    OPTWRK = MAX( N+LWRK_SGEQP3, LWRK_SGESVD )
                 }
             }
         } else if ( LSVEC && (.NOT.RSVEC) ) {
             // .. minimal and optimal sizes of the workspace if the
             // singular values and the left singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD, LWORQ )
             } else {
                 MINWRK = N + MAX( LWQP3, LWSVD, LWORQ )
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   sgesvd('N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                } else {
                   sgesvd('O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                }
                LWRK_SGESVD = INT( RDUMMY(1) )
                if ( CONDA ) {
                    OPTWRK = N + MAX( LWRK_SGEQP3, LWCON, LWRK_SGESVD, LWRK_SORMQR )
                } else {
                    OPTWRK = N + MAX( LWRK_SGEQP3, LWRK_SGESVD, LWRK_SORMQR )
                }
             }
         } else if ( RSVEC && (.NOT.LSVEC) ) {
             // .. minimal and optimal sizes of the workspace if the
             // singular values and the right singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD )
             } else {
                 MINWRK = N + MAX( LWQP3, LWSVD )
             }
             if ( LQUERY ) {
                 if ( RTRANS ) {
                     sgesvd('O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                 } else {
                     sgesvd('N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                 }
                 LWRK_SGESVD = INT( RDUMMY(1) )
                 if ( CONDA ) {
                     OPTWRK = N + MAX( LWRK_SGEQP3, LWCON, LWRK_SGESVD )
                 } else {
                     OPTWRK = N + MAX( LWRK_SGEQP3, LWRK_SGESVD )
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
                    // .. minimal workspace length for N x N/2 SGEQRF
                    LWQRF  = MAX( N/2, 1 )
                    // .. minimal workspace length for N/2 x N/2 SGESVD
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
                    // .. minimal workspace length for N/2 x N SGELQF
                    LWLQF  = MAX( N/2, 1 )
                    LWSVD2 = MAX( 5 * (N/2), 1 )
                    LWUNLQ = MAX( N , 1 )
                    MINWRK2 = MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWORQ )
                    if (CONDA) MINWRK2 = MAX( MINWRK2, LWCON );
                    MINWRK2 = N + MINWRK2
                    MINWRK = MAX( MINWRK, MINWRK2 )
                 }
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   sgesvd('O', 'A', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                   LWRK_SGESVD = INT( RDUMMY(1) )
                   OPTWRK = MAX(LWRK_SGEQP3,LWRK_SGESVD,LWRK_SORMQR)
                   if (CONDA) OPTWRK = MAX( OPTWRK, LWCON );
                   OPTWRK = N + OPTWRK
                   if ( WNTVA ) {
                       sgeqrf(N,N/2,U,LDU,RDUMMY,RDUMMY,-1,IERR);
                       LWRK_SGEQRF = INT( RDUMMY(1) )
                       sgesvd('S', 'O', N/2,N/2, V,LDV, S, U,LDU, V, LDV, RDUMMY, -1, IERR );
                       LWRK_SGESVD2 = INT( RDUMMY(1) )
                       sormqr('R', 'C', N, N, N/2, U, LDU, RDUMMY, V, LDV, RDUMMY, -1, IERR );
                       LWRK_SORMQR2 = INT( RDUMMY(1) )
                       OPTWRK2 = MAX( LWRK_SGEQP3, N/2+LWRK_SGEQRF, N/2+LWRK_SGESVD2, N/2+LWRK_SORMQR2 )
                       if (CONDA) OPTWRK2 = MAX( OPTWRK2, LWCON );
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   }
                } else {
                   sgesvd('S', 'O', N, N, A, LDA, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                   LWRK_SGESVD = INT( RDUMMY(1) )
                   OPTWRK = MAX(LWRK_SGEQP3,LWRK_SGESVD,LWRK_SORMQR)
                   if (CONDA) OPTWRK = MAX( OPTWRK, LWCON );
                   OPTWRK = N + OPTWRK
                   if ( WNTVA ) {
                      sgelqf(N/2,N,U,LDU,RDUMMY,RDUMMY,-1,IERR);
                      LWRK_SGELQF = INT( RDUMMY(1) )
                      sgesvd('S','O', N/2,N/2, V, LDV, S, U, LDU, V, LDV, RDUMMY, -1, IERR );
                      LWRK_SGESVD2 = INT( RDUMMY(1) )
                      sormlq('R', 'N', N, N, N/2, U, LDU, RDUMMY, V, LDV, RDUMMY,-1,IERR );
                      LWRK_SORMLQ = INT( RDUMMY(1) )
                      OPTWRK2 = MAX( LWRK_SGEQP3, N/2+LWRK_SGELQF, N/2+LWRK_SGESVD2, N/2+LWRK_SORMLQ )
                       if (CONDA) OPTWRK2 = MAX( OPTWRK2, LWCON );
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   }
                }
             }
         }

         MINWRK = MAX( 2, MINWRK )
         OPTWRK = MAX( 2, OPTWRK )
         IF ( LWORK < MINWRK && (.NOT.LQUERY) ) INFO = -19

      }

      if (INFO == 0 && LRWORK < RMINWRK && .NOT. LQUERY) {
         INFO = -21
      }
      if ( INFO != 0 ) {
         xerbla('SGESVDQ', -INFO );
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

      if ( ( M == 0 ) || ( N == 0 ) ) {
      // .. all output is void.
         RETURN
      }

      BIG = SLAMCH('O')
      ASCALED = false;
      IWOFF = 1
      if ( ROWPRM ) {
            IWOFF = M
            // .. reordering the rows in decreasing sequence in the
            // ell-infinity norm - this enhances numerical robustness in
            // the case of differently scaled rows.
            for (p = 1; p <= M; p++) { // 1904
                // RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
                // [[SLANGE will return NaN if an entry of the p-th row is Nan]]
                RWORK(p) = SLANGE( 'M', 1, N, A(p,1), LDA, RDUMMY )
                // .. check for NaN's and Inf's
                if ( ( RWORK(p) != RWORK(p) ) || ( (RWORK(p)*ZERO) != ZERO ) ) {
                    INFO = -8
                    xerbla('SGESVDQ', -INFO );
                    RETURN
                }
            } // 1904
            for (p = 1; p <= M - 1; p++) { // 1952
            q = ISAMAX( M-p+1, RWORK(p), 1 ) + p - 1
            IWORK(N+p) = q
            if ( p != q ) {
               RTMP     = RWORK(p)
               RWORK(p) = RWORK(q)
               RWORK(q) = RTMP
            }
            } // 1952

            if ( RWORK(1) == ZERO ) {
               // Quick return: A is the M x N zero matrix.
               NUMRANK = 0
               slaset('G', N, 1, ZERO, ZERO, S, N );
               if (WNTUS) CALL SLASET('G', M, N, ZERO, ONE, U, LDU);
               if (WNTUA) CALL SLASET('G', M, M, ZERO, ONE, U, LDU);
               if (WNTVA) CALL SLASET('G', N, N, ZERO, ONE, V, LDV);
               if ( WNTUF ) {
                   slaset('G', N, 1, ZERO, ZERO, WORK, N );
                   slaset('G', M, N, ZERO,  ONE, U, LDU );
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

            if ( RWORK(1) .GT. BIG / SQRT(REAL(M)) ) {
                // .. to prevent overflow in the QR factorization, scale the
                // matrix by 1/sqrt(M) if too large entry detected
                slascl('G',0,0,SQRT(REAL(M)),ONE, M,N, A,LDA, IERR);
                ASCALED = true;
            }
            slaswp(N, A, LDA, 1, M-1, IWORK(N+1), 1 );
      }

*    .. At this stage, preemptive scaling is done only to avoid column
*    norms overflows during the QR factorization. The SVD procedure should
*    have its own scaling to save the singular values from overflows and
*    underflows. That depends on the SVD procedure.

      if ( .NOT.ROWPRM ) {
          RTMP = SLANGE( 'M', M, N, A, LDA, RDUMMY )
          if ( ( RTMP != RTMP ) || ( (RTMP*ZERO) != ZERO ) ) {
               INFO = -8
               xerbla('SGESVDQ', -INFO );
               RETURN
          }
          if ( RTMP .GT. BIG / SQRT(REAL(M)) ) {
              // .. to prevent overflow in the QR factorization, scale the
              // matrix by 1/sqrt(M) if too large entry detected
              slascl('G',0,0, SQRT(REAL(M)),ONE, M,N, A,LDA, IERR);
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
      sgeqp3(M, N, A, LDA, IWORK, WORK, WORK(N+1), LWORK-N, IERR );

*    If the user requested accuracy level allows truncation in the
*    computed upper triangular factor, the matrix R is examined and,
*    if possible, replaced with its leading upper trapezoidal part.

      EPSLN = SLAMCH('E')
      SFMIN = SLAMCH('S')
      // SMALL = SFMIN / EPSLN
      NR = N

      if ( ACCLA ) {

         // Standard absolute error bound suffices. All sigma_i with
         // sigma_i < N*EPS*||A||_F are flushed to zero. This is an
         // aggressive enforcement of lower numerical rank by introducing a
         // backward error of the order of N*EPS*||A||_F.
         NR = 1
         RTMP = SQRT(REAL(N))*EPSLN
         for (p = 2; p <= N; p++) { // 3001
            IF ( ABS(A(p,p)) < (RTMP*ABS(A(1,1))) ) GO TO 3002
               NR = NR + 1
         } // 3001
         } // 3002

      } else if ( ACCLM ) {
         // .. similarly as above, only slightly more gentle (less aggressive).
         // Sudden drop on the diagonal of R is used as the criterion for being
         // close-to-rank-deficient. The threshold is set to EPSLN=SLAMCH('E').
         // [[This can be made more flexible by replacing this hard-coded value
         // with a user specified threshold.]] Also, the values that underflow
         // will be truncated.
         NR = 1
         for (p = 2; p <= N; p++) { // 3401
            IF ( ( ABS(A(p,p)) < (EPSLN*ABS(A(p-1,p-1))) ) || ( ABS(A(p,p)) < SFMIN ) ) GO TO 3402
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
            IF ( ABS(A(p,p)) == ZERO ) GO TO 3502
            NR = NR + 1
         } // 3501
         } // 3502

         if ( CONDA ) {
            // Estimate the scaled condition number of A. Use the fact that it is
            // the same as the scaled condition number of R.
               // .. V is used as workspace
               slacpy('U', N, N, A, LDA, V, LDV );
               // Only the leading NR x NR submatrix of the triangular factor
               // is considered. Only if NR=N will this give a reliable error
               // bound. However, even for NR < N, this can be used on an
               // expert level and obtain useful information in the sense of
               // perturbation theory.
               for (p = 1; p <= NR; p++) { // 3053
                  RTMP = SNRM2( p, V(1,p), 1 )
                  sscal(p, ONE/RTMP, V(1,p), 1 );
               } // 3053
               if ( .NOT. ( LSVEC || RSVEC ) ) {
                   spocon('U', NR, V, LDV, ONE, RTMP, WORK, IWORK(N+IWOFF), IERR );
               } else {
                   spocon('U', NR, V, LDV, ONE, RTMP, WORK(N+1), IWORK(N+IWOFF), IERR );
               }
               SCONDA = ONE / SQRT(RTMP)
            // For NR=N, SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1),
            // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
            // See the reference [1] for more details.
         }

      }

      if ( WNTUR ) {
          N1 = NR
      } else if ( WNTUS || WNTUF) {
          N1 = N
      } else if ( WNTUA ) {
          N1 = M
      }

      if ( .NOT. ( RSVEC || LSVEC ) ) {
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

            sgesvd('N', 'N', N, NR, A, LDA, S, U, LDU, V, LDV, WORK, LWORK, INFO );

         } else {

            // .. compute the singular values of R = [A](1:NR,1:N)

            if (NR .GT. 1) CALL SLASET( 'L', NR-1,NR-1, ZERO,ZERO, A(2,1), LDA );
            sgesvd('N', 'N', NR, N, A, LDA, S, U, LDU, V, LDV, WORK, LWORK, INFO );

         }

      } else if ( LSVEC && ( .NOT. RSVEC) ) {
*.......................................................................
        // .. the singular values and the left singular vectors requested
*.......................................................................""""""""
         if ( RTRANS ) {
             // .. apply SGESVD to R**T
             // .. copy R**T into [U] and overwrite [U] with the right singular
             // vectors of R
            for (p = 1; p <= NR; p++) { // 1192
               for (q = p; q <= N; q++) { // 1193
                  U(q,p) = A(p,q)
               } // 1193
            } // 1192
            if (NR .GT. 1) CALL SLASET( 'U', NR-1,NR-1, ZERO,ZERO, U(1,2), LDU );
            // .. the left singular vectors not computed, the NR right singular
            // vectors overwrite [U](1:NR,1:NR) as transposed. These
            // will be pre-multiplied by Q to build the left singular vectors of A.
               sgesvd('N', 'O', N, NR, U, LDU, S, U, LDU, U, LDU, WORK(N+1), LWORK-N, INFO );

               for (p = 1; p <= NR; p++) { // 1119
                   for (q = p + 1; q <= NR; q++) { // 1120
                      RTMP   = U(q,p)
                      U(q,p) = U(p,q)
                      U(p,q) = RTMP
                   } // 1120
               } // 1119

         } else {
             // .. apply SGESVD to R
             // .. copy R into [U] and overwrite [U] with the left singular vectors
             slacpy('U', NR, N, A, LDA, U, LDU );
             if (NR .GT. 1) CALL SLASET( 'L', NR-1, NR-1, ZERO, ZERO, U(2,1), LDU );
             // .. the right singular vectors not computed, the NR left singular
             // vectors overwrite [U](1:NR,1:NR)
                sgesvd('O', 'N', NR, N, U, LDU, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                // .. now [U](1:NR,1:NR) contains the NR left singular vectors of
                // R. These will be pre-multiplied by Q to build the left singular
                // vectors of A.
         }

            // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
         if ( ( NR < M ) && ( .NOT.WNTUF ) ) {
             slaset('A', M-NR, NR, ZERO, ZERO, U(NR+1,1), LDU);
             if ( NR < N1 ) {
                slaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1), LDU );
                slaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
             }
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         if (.NOT.WNTUF) CALL SORMQR( 'L', 'N', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );
         if (ROWPRM && .NOT.WNTUF) CALL SLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 );

      } else if ( RSVEC && ( .NOT. LSVEC ) ) {
*.......................................................................
        // .. the singular values and the right singular vectors requested
*.......................................................................
          if ( RTRANS ) {
             // .. apply SGESVD to R**T
             // .. copy R**T into V and overwrite V with the left singular vectors
            for (p = 1; p <= NR; p++) { // 1165
               for (q = p; q <= N; q++) { // 1166
                  V(q,p) = (A(p,q))
               } // 1166
            } // 1165
            if (NR .GT. 1) CALL SLASET( 'U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV );
            // .. the left singular vectors of R**T overwrite V, the right singular
            // vectors not computed
            if ( WNTVR || ( NR == N ) ) {
               sgesvd('O', 'N', N, NR, V, LDV, S, U, LDU, U, LDU, WORK(N+1), LWORK-N, INFO );

               for (p = 1; p <= NR; p++) { // 1121
                   for (q = p + 1; q <= NR; q++) { // 1122
                      RTMP   = V(q,p)
                      V(q,p) = V(p,q)
                      V(p,q) = RTMP
                   } // 1122
               } // 1121

               if ( NR < N ) {
                   for (p = 1; p <= NR; p++) { // 1103
                      for (q = NR + 1; q <= N; q++) { // 1104
                          V(p,q) = V(q,p)
                      } // 1104
                   } // 1103
               }
               slapmt( false , NR, N, V, LDV, IWORK );
            } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:N,1:NR)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the QR factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                slaset('G', N, N-NR, ZERO, ZERO, V(1,NR+1), LDV);
                sgesvd('O', 'N', N, N, V, LDV, S, U, LDU, U, LDU, WORK(N+1), LWORK-N, INFO );

                for (p = 1; p <= N; p++) { // 1123
                   for (q = p + 1; q <= N; q++) { // 1124
                      RTMP   = V(q,p)
                      V(q,p) = V(p,q)
                      V(p,q) = RTMP
                   } // 1124
                } // 1123
                slapmt( false , N, N, V, LDV, IWORK );
            }

          } else {
             // .. aply SGESVD to R
             // .. copy R into V and overwrite V with the right singular vectors
             slacpy('U', NR, N, A, LDA, V, LDV );
             if (NR .GT. 1) CALL SLASET( 'L', NR-1, NR-1, ZERO, ZERO, V(2,1), LDV );
             // .. the right singular vectors overwrite V, the NR left singular
             // vectors stored in U(1:NR,1:NR)
             if ( WNTVR || ( NR == N ) ) {
                sgesvd('N', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                slapmt( false , NR, N, V, LDV, IWORK );
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
             } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:NR,1:N)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the LQ factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                 slaset('G', N-NR, N, ZERO,ZERO, V(NR+1,1), LDV);
                 sgesvd('N', 'O', N, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                 slapmt( false , N, N, V, LDV, IWORK );
             }
             // .. now [V] contains the transposed matrix of the right singular
             // vectors of A.
          }

      } else {
*.......................................................................
        // .. FULL SVD requested
*.......................................................................
         if ( RTRANS ) {

             // .. apply SGESVD to R**T [[this option is left for R&D&T]]

            if ( WNTVR || ( NR == N ) ) {
             // .. copy R**T into [V] and overwrite [V] with the left singular
             // vectors of R**T
            for (p = 1; p <= NR; p++) { // 1168
               for (q = p; q <= N; q++) { // 1169
                  V(q,p) = A(p,q)
               } // 1169
            } // 1168
            if (NR .GT. 1) CALL SLASET( 'U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV );

            // .. the left singular vectors of R**T overwrite [V], the NR right
            // singular vectors of R**T stored in [U](1:NR,1:NR) as transposed
               sgesvd('O', 'A', N, NR, V, LDV, S, V, LDV, U, LDU, WORK(N+1), LWORK-N, INFO );
               // .. assemble V
               for (p = 1; p <= NR; p++) { // 1115
                  for (q = p + 1; q <= NR; q++) { // 1116
                     RTMP   = V(q,p)
                     V(q,p) = V(p,q)
                     V(p,q) = RTMP
                  } // 1116
               } // 1115
               if ( NR < N ) {
                   for (p = 1; p <= NR; p++) { // 1101
                      for (q = NR+1; q <= N; q++) { // 1102
                         V(p,q) = V(q,p)
                      } // 1102
                   } // 1101
               }
               slapmt( false , NR, N, V, LDV, IWORK );

                for (p = 1; p <= NR; p++) { // 1117
                   for (q = p + 1; q <= NR; q++) { // 1118
                      RTMP   = U(q,p)
                      U(q,p) = U(p,q)
                      U(p,q) = RTMP
                   } // 1118
                } // 1117

                if ( ( NR < M ) && .NOT.(WNTUF)) {
                  slaset('A', M-NR,NR, ZERO,ZERO, U(NR+1,1), LDU);
                  if ( NR < N1 ) {
                     slaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     slaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
                  }
               }

            } else {
                // .. need all N right singular vectors and NR < N
             // .. copy R**T into [V] and overwrite [V] with the left singular
             // vectors of R**T
                // [[The optimal ratio N/NR for using QRF instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'SGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = MAX( OPTRATIO, 2 )
                OPTRATIO = 2
                if ( OPTRATIO*NR .GT. N ) {
                   for (p = 1; p <= NR; p++) { // 1198
                      for (q = p; q <= N; q++) { // 1199
                         V(q,p) = A(p,q)
                      } // 1199
                   } // 1198
                   if (NR .GT. 1) CALL SLASET('U',NR-1,NR-1, ZERO,ZERO, V(1,2),LDV);

                   slaset('A',N,N-NR,ZERO,ZERO,V(1,NR+1),LDV);
                   sgesvd('O', 'A', N, N, V, LDV, S, V, LDV, U, LDU, WORK(N+1), LWORK-N, INFO );

                   for (p = 1; p <= N; p++) { // 1113
                      for (q = p + 1; q <= N; q++) { // 1114
                         RTMP   = V(q,p)
                         V(q,p) = V(p,q)
                         V(p,q) = RTMP
                      } // 1114
                   } // 1113
                   slapmt( false , N, N, V, LDV, IWORK );
               // .. assemble the left singular vector matrix U of dimensions
               // (M x N1), i.e. (M x N) or (M x M).

                   for (p = 1; p <= N; p++) { // 1111
                      for (q = p + 1; q <= N; q++) { // 1112
                         RTMP   = U(q,p)
                         U(q,p) = U(p,q)
                         U(p,q) = RTMP
                      } // 1112
                   } // 1111

                   if ( ( N < M ) && .NOT.(WNTUF)) {
                      slaset('A',M-N,N,ZERO,ZERO,U(N+1,1),LDU);
                      if ( N < N1 ) {
                        slaset('A',N,N1-N,ZERO,ZERO,U(1,N+1),LDU);
                        slaset('A',M-N,N1-N,ZERO,ONE, U(N+1,N+1), LDU );
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
                   if (NR .GT. 1) CALL SLASET('U',NR-1,NR-1,ZERO,ZERO,U(1,NR+2),LDU);
                   sgeqrf(N, NR, U(1,NR+1), LDU, WORK(N+1), WORK(N+NR+1), LWORK-N-NR, IERR );
                   for (p = 1; p <= NR; p++) { // 1143
                       for (q = 1; q <= N; q++) { // 1144
                           V(q,p) = U(p,NR+q)
                       } // 1144
                   } // 1143
                  slaset('U',NR-1,NR-1,ZERO,ZERO,V(1,2),LDV);
                  sgesvd('S', 'O', NR, NR, V, LDV, S, U, LDU, V,LDV, WORK(N+NR+1),LWORK-N-NR, INFO );
                  slaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV);
                  slaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV);
                  slaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV);
                  sormqr('R','C', N, N, NR, U(1,NR+1), LDU, WORK(N+1),V,LDV,WORK(N+NR+1),LWORK-N-NR,IERR);
                  slapmt( false , N, N, V, LDV, IWORK );
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x NR) or (M x N) or (M x M).
                  if ( ( NR < M ) && .NOT.(WNTUF)) {
                     slaset('A',M-NR,NR,ZERO,ZERO,U(NR+1,1),LDU);
                     if ( NR < N1 ) {
                     slaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     slaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1),LDU);
                     }
                  }
                }
            }

         } else {

             // .. apply SGESVD to R [[this is the recommended option]]

             if ( WNTVR || ( NR == N ) ) {
                 // .. copy R into [V] and overwrite V with the right singular vectors
                 slacpy('U', NR, N, A, LDA, V, LDV );
                if (NR .GT. 1) CALL SLASET( 'L', NR-1,NR-1, ZERO,ZERO, V(2,1), LDV );
                // .. the right singular vectors of R overwrite [V], the NR left
                // singular vectors of R stored in [U](1:NR,1:NR)
                sgesvd('S', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                slapmt( false , NR, N, V, LDV, IWORK );
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
               if ( ( NR < M ) && .NOT.(WNTUF)) {
                  slaset('A', M-NR,NR, ZERO,ZERO, U(NR+1,1), LDU);
                  if ( NR < N1 ) {
                     slaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     slaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
                  }
               }

             } else {
               // .. need all N right singular vectors and NR < N
               // .. the requested number of the left singular vectors
                // is then N1 (N or M)
                // [[The optimal ratio N/NR for using LQ instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'SGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = MAX( OPTRATIO, 2 )
               OPTRATIO = 2
               if ( OPTRATIO * NR .GT. N ) {
                  slacpy('U', NR, N, A, LDA, V, LDV );
                  if (NR .GT. 1) CALL SLASET('L', NR-1,NR-1, ZERO,ZERO, V(2,1),LDV);
               // .. the right singular vectors of R overwrite [V], the NR left
                  // singular vectors of R stored in [U](1:NR,1:NR)
                  slaset('A', N-NR,N, ZERO,ZERO, V(NR+1,1),LDV);
                  sgesvd('S', 'O', N, N, V, LDV, S, U, LDU, V, LDV, WORK(N+1), LWORK-N, INFO );
                  slapmt( false , N, N, V, LDV, IWORK );
                  // .. now [V] contains the transposed matrix of the right
                  // singular vectors of A. The leading N left singular vectors
                  // are in [U](1:N,1:N)
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x N1), i.e. (M x N) or (M x M).
                  if ( ( N < M ) && .NOT.(WNTUF)) {
                      slaset('A',M-N,N,ZERO,ZERO,U(N+1,1),LDU);
                      if ( N < N1 ) {
                        slaset('A',N,N1-N,ZERO,ZERO,U(1,N+1),LDU);
                        slaset('A',M-N,N1-N,ZERO,ONE, U(N+1,N+1), LDU );
                      }
                  }
               } else {
                  slacpy('U', NR, N, A, LDA, U(NR+1,1), LDU );
                  if (NR .GT. 1) CALL SLASET('L',NR-1,NR-1,ZERO,ZERO,U(NR+2,1),LDU);
                  sgelqf(NR, N, U(NR+1,1), LDU, WORK(N+1), WORK(N+NR+1), LWORK-N-NR, IERR );
                  slacpy('L',NR,NR,U(NR+1,1),LDU,V,LDV);
                  if (NR .GT. 1) CALL SLASET('U',NR-1,NR-1,ZERO,ZERO,V(1,2),LDV);
                  sgesvd('S', 'O', NR, NR, V, LDV, S, U, LDU, V, LDV, WORK(N+NR+1), LWORK-N-NR, INFO );
                  slaset('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV);
                  slaset('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV);
                  slaset('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV);
                  sormlq('R','N',N,N,NR,U(NR+1,1),LDU,WORK(N+1), V, LDV, WORK(N+NR+1),LWORK-N-NR,IERR);
                  slapmt( false , N, N, V, LDV, IWORK );
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
                  if ( ( NR < M ) && .NOT.(WNTUF)) {
                     slaset('A',M-NR,NR,ZERO,ZERO,U(NR+1,1),LDU);
                     if ( NR < N1 ) {
                     slaset('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU);
                     slaset('A',M-NR,N1-NR,ZERO,ONE, U(NR+1,NR+1), LDU );
                     }
                  }
               }
             }
         // .. end of the "R**T or R" branch
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         if (.NOT. WNTUF) CALL SORMQR( 'L', 'N', M, N1, N, A, LDA, WORK, U, LDU, WORK(N+1), LWORK-N, IERR );
         if (ROWPRM && .NOT.WNTUF) CALL SLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 );

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
      if (NR < N) CALL SLASET( 'G', N-NR,1, ZERO,ZERO, S(NR+1), N );
      // .. undo scaling; this may cause overflow in the largest singular
      // values.
      if (ASCALED) CALL SLASCL( 'G',0,0, ONE,SQRT(REAL(M)), NR,1, S, N, IERR );
      if (CONDA) RWORK(1) = SCONDA;
      RWORK(2) = p - NR
      // .. p-NR is the number of singular values that are computed as
      // exact zeros in SGESVD() applied to the (possibly truncated)
      // full row rank triangular (trapezoidal) factor of A.
      NUMRANK = NR

      RETURN

      // End of SGESVDQ

      }
