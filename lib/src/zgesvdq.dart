      void zgesvdq(JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK, CWORK, LCWORK, RWORK, LRWORK, INFO ) {
      // .. Scalar Arguments ..
      IMPLICIT    NONE;
      String      JOBA, JOBP, JOBR, JOBU, JOBV;
      int         M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LCWORK, LRWORK, INFO;
      // ..
      // .. Array Arguments ..
      Complex       A( LDA, * ), U( LDU, * ), V( LDV, * ), CWORK( * );
      double           S( * ), RWORK( * );
      int              IWORK( * );

// =====================================================================

      // .. Parameters ..
      double           ZERO,         ONE;
      const          ZERO = 0.0, ONE = 1.0 ;
      Complex       CZERO,                 CONE;
      const          CZERO = (0.0,0.0), CONE = (1.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int         IERR, NR, N1, OPTRATIO, p, q;
      int         LWCON, LWQP3, LWRK_ZGELQF, LWRK_ZGESVD, LWRK_ZGESVD2, LWRK_ZGEQP3, LWRK_ZGEQRF, LWRK_ZUNMLQ, LWRK_ZUNMQR, LWRK_ZUNMQR2, LWLQF, LWQRF, LWSVD, LWSVD2, LWUNQ, LWUNQ2, LWUNLQ, MINWRK, MINWRK2, OPTWRK, OPTWRK2, IMINWRK, RMINWRK;
      bool        ACCLA,  ACCLM, ACCLH, ASCALED, CONDA, DNTWU,  DNTWV, LQUERY, LSVC0, LSVEC, ROWPRM,  RSVEC, RTRANS, WNTUA, WNTUF,  WNTUR, WNTUS, WNTVA,   WNTVR;
      double           BIG, EPSLN, RTMP, SCONDA, SFMIN;
      Complex       CTMP;
      // ..
      // .. Local Arrays
      Complex         CDUMMY(1);
      double             RDUMMY(1);
      // ..
      // .. External Subroutines (BLAS, LAPACK)
      // EXTERNAL ZGELQF, ZGEQP3, ZGEQRF, ZGESVD, ZLACPY, ZLAPMT, ZLASCL, ZLASET, ZLASWP, ZDSCAL, DLASET, DLASCL, ZPOCON, ZUNMLQ, ZUNMQR, XERBLA
      // ..
      // .. External Functions (BLAS, LAPACK)
      //- bool        lsame;
      //- int                         idamax;
      //- double             ZLANGE,          DZNRM2, DLAMCH;
      // EXTERNAL lsame, ZLANGE,  idamax, DZNRM2, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, DBLE, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      WNTUS  = lsame( JOBU, 'S' ) || lsame( JOBU, 'U' );
      WNTUR  = lsame( JOBU, 'R' );
      WNTUA  = lsame( JOBU, 'A' );
      WNTUF  = lsame( JOBU, 'F' );
      LSVC0  = WNTUS || WNTUR || WNTUA;
      LSVEC  = LSVC0 || WNTUF;
      DNTWU  = lsame( JOBU, 'N' );

      WNTVR  = lsame( JOBV, 'R' );
      WNTVA  = lsame( JOBV, 'A' ) || lsame( JOBV, 'V' );
      RSVEC  = WNTVR || WNTVA;
      DNTWV  = lsame( JOBV, 'N' );

      ACCLA  = lsame( JOBA, 'A' );
      ACCLM  = lsame( JOBA, 'M' );
      CONDA  = lsame( JOBA, 'E' );
      ACCLH  = lsame( JOBA, 'H' ) || CONDA;

      ROWPRM = lsame( JOBP, 'P' );
      RTRANS = lsame( JOBR, 'T' );

      if ( ROWPRM ) {
         IMINWRK = max( 1, N + M - 1 );
         RMINWRK = max( 2, M, 5*N );
      } else {
         IMINWRK = max( 1, N );
         RMINWRK = max( 2, 5*N );
      }
      LQUERY = (LIWORK == -1 || LCWORK == -1 || LRWORK == -1);
      INFO  = 0;
      if ( !( ACCLA || ACCLM || ACCLH ) ) {
         INFO = -1;
      } else if ( !( ROWPRM || lsame( JOBP, 'N' ) ) ) {
          INFO = -2;
      } else if ( !( RTRANS || lsame( JOBR, 'N' ) ) ) {
          INFO = -3;
      } else if ( !( LSVEC || DNTWU ) ) {
         INFO = -4;
      } else if ( WNTUR && WNTVA ) {
         INFO = -5;
      } else if ( !( RSVEC || DNTWV )) {
         INFO = -5;
      } else if ( M < 0 ) {
         INFO = -6;
      } else if ( ( N < 0 ) || ( N > M ) ) {
         INFO = -7;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -9;
      } else if ( LDU < 1 || ( LSVC0 && LDU < M ) || ( WNTUF && LDU < N ) ) {
         INFO = -12;
      } else if ( LDV < 1 || ( RSVEC && LDV < N ) || ( CONDA && LDV < N ) ) {
         INFO = -14;
      } else if ( LIWORK < IMINWRK && !LQUERY ) {
         INFO = -17;
      }


      if ( INFO == 0 ) {
         // .. compute the minimal and the optimal workspace lengths
         // [[The expressions for computing the minimal and the optimal
         // values of LCWORK are written with a lot of redundancy and
         // can be simplified. However, this detailed form is easier for
         // maintenance and modifications of the code.]]

         // .. minimal workspace length for ZGEQP3 of an M x N matrix
         LWQP3 = N+1;
         // .. minimal workspace length for ZUNMQR to build left singular vectors
         if ( WNTUS || WNTUR ) {
             LWUNQ  = max( N  , 1 );
         } else if ( WNTUA ) {
             LWUNQ = max( M , 1 );
         }
         // .. minimal workspace length for ZPOCON of an N x N matrix
         LWCON = 2 * N;
         // .. ZGESVD of an N x N matrix
         LWSVD = max( 3 * N, 1 );
         if ( LQUERY ) {
             zgeqp3(M, N, A, LDA, IWORK, CDUMMY, CDUMMY, -1, RDUMMY, IERR );
             LWRK_ZGEQP3 = INT( CDUMMY(1) );
             if ( WNTUS || WNTUR ) {
                 zunmqr('L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR );
                 LWRK_ZUNMQR = INT( CDUMMY(1) );
             } else if ( WNTUA ) {
                 zunmqr('L', 'N', M, M, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR );
                 LWRK_ZUNMQR = INT( CDUMMY(1) );
             } else {
                 LWRK_ZUNMQR = 0;
             }
         }
         MINWRK = 2;
         OPTWRK = 2;
         if ( !(LSVEC || RSVEC ) ) {
             // .. minimal and optimal sizes of the complex workspace if
             // only the singular values are requested
             if ( CONDA ) {
                MINWRK = max( N+LWQP3, LWCON, LWSVD );
             } else {
                MINWRK = max( N+LWQP3, LWSVD );
             }
             if ( LQUERY ) {
                 zgesvd('N', 'N', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                 LWRK_ZGESVD = INT( CDUMMY(1) );
                 if ( CONDA ) {
                    OPTWRK = max( N+LWRK_ZGEQP3, N+LWCON, LWRK_ZGESVD );
                 } else {
                    OPTWRK = max( N+LWRK_ZGEQP3, LWRK_ZGESVD );
                 }
             }
         } else if ( LSVEC && ( !RSVEC) ) {
             // .. minimal and optimal sizes of the complex workspace if the
             // singular values and the left singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + max( LWQP3, LWCON, LWSVD, LWUNQ );
             } else {
                 MINWRK = N + max( LWQP3, LWSVD, LWUNQ );
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   zgesvd('N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                } else {
                   zgesvd('O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                }
                LWRK_ZGESVD = INT( CDUMMY(1) );
                if ( CONDA ) {
                    OPTWRK = N + max( LWRK_ZGEQP3, LWCON, LWRK_ZGESVD, LWRK_ZUNMQR );
                } else {
                    OPTWRK = N + max( LWRK_ZGEQP3, LWRK_ZGESVD, LWRK_ZUNMQR );
                }
             }
         } else if ( RSVEC && ( !LSVEC) ) {
             // .. minimal and optimal sizes of the complex workspace if the
             // singular values and the right singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + max( LWQP3, LWCON, LWSVD );
             } else {
                 MINWRK = N + max( LWQP3, LWSVD );
             }
             if ( LQUERY ) {
                 if ( RTRANS ) {
                     zgesvd('O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                 } else {
                     zgesvd('N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                 }
                 LWRK_ZGESVD = INT( CDUMMY(1) );
                 if ( CONDA ) {
                     OPTWRK = N + max( LWRK_ZGEQP3, LWCON, LWRK_ZGESVD );
                 } else {
                     OPTWRK = N + max( LWRK_ZGEQP3, LWRK_ZGESVD );
                 }
             }
         } else {
             // .. minimal and optimal sizes of the complex workspace if the
             // full SVD is requested
             if ( RTRANS ) {
                 MINWRK = max( LWQP3, LWSVD, LWUNQ );
                 if (CONDA) MINWRK = max( MINWRK, LWCON );
                 MINWRK = MINWRK + N;
                 if ( WNTVA ) {
                    // .. minimal workspace length for N x N/2 ZGEQRF
                    LWQRF  = max( N/2, 1 );
                    // .. minimal workspace length for N/2 x N/2 ZGESVD
                    LWSVD2 = max( 3 * (N/2), 1 );
                    LWUNQ2 = max( N, 1 );
                    MINWRK2 = max( LWQP3, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ );
                    if (CONDA) MINWRK2 = max( MINWRK2, LWCON );
                    MINWRK2 = N + MINWRK2;
                    MINWRK = max( MINWRK, MINWRK2 );
                 }
             } else {
                 MINWRK = max( LWQP3, LWSVD, LWUNQ );
                 if (CONDA) MINWRK = max( MINWRK, LWCON );
                 MINWRK = MINWRK + N;
                 if ( WNTVA ) {
                    // .. minimal workspace length for N/2 x N ZGELQF
                    LWLQF  = max( N/2, 1 );
                    LWSVD2 = max( 3 * (N/2), 1 );
                    LWUNLQ = max( N , 1 );
                    MINWRK2 = max( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ );
                    if (CONDA) MINWRK2 = max( MINWRK2, LWCON );
                    MINWRK2 = N + MINWRK2;
                    MINWRK = max( MINWRK, MINWRK2 );
                 }
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   zgesvd('O', 'A', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                   LWRK_ZGESVD = INT( CDUMMY(1) );
                   OPTWRK = max(LWRK_ZGEQP3,LWRK_ZGESVD,LWRK_ZUNMQR);
                   if (CONDA) OPTWRK = max( OPTWRK, LWCON );
                   OPTWRK = N + OPTWRK;
                   if ( WNTVA ) {
                       zgeqrf(N,N/2,U,LDU,CDUMMY,CDUMMY,-1,IERR);
                       LWRK_ZGEQRF = INT( CDUMMY(1) );
                       zgesvd('S', 'O', N/2,N/2, V,LDV, S, U,LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                       LWRK_ZGESVD2 = INT( CDUMMY(1) );
                       zunmqr('R', 'C', N, N, N/2, U, LDU, CDUMMY, V, LDV, CDUMMY, -1, IERR );
                       LWRK_ZUNMQR2 = INT( CDUMMY(1) );
                       OPTWRK2 = max( LWRK_ZGEQP3, N/2+LWRK_ZGEQRF, N/2+LWRK_ZGESVD2, N/2+LWRK_ZUNMQR2 );
                       if (CONDA) OPTWRK2 = max( OPTWRK2, LWCON );
                       OPTWRK2 = N + OPTWRK2;
                       OPTWRK = max( OPTWRK, OPTWRK2 );
                   }
                } else {
                   zgesvd('S', 'O', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                   LWRK_ZGESVD = INT( CDUMMY(1) );
                   OPTWRK = max(LWRK_ZGEQP3,LWRK_ZGESVD,LWRK_ZUNMQR);
                   if (CONDA) OPTWRK = max( OPTWRK, LWCON );
                   OPTWRK = N + OPTWRK;
                   if ( WNTVA ) {
                      zgelqf(N/2,N,U,LDU,CDUMMY,CDUMMY,-1,IERR);
                      LWRK_ZGELQF = INT( CDUMMY(1) );
                      zgesvd('S','O', N/2,N/2, V, LDV, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR );
                      LWRK_ZGESVD2 = INT( CDUMMY(1) );
                      zunmlq('R', 'N', N, N, N/2, U, LDU, CDUMMY, V, LDV, CDUMMY,-1,IERR );
                      LWRK_ZUNMLQ = INT( CDUMMY(1) );
                      OPTWRK2 = max( LWRK_ZGEQP3, N/2+LWRK_ZGELQF, N/2+LWRK_ZGESVD2, N/2+LWRK_ZUNMLQ );
                       if (CONDA) OPTWRK2 = max( OPTWRK2, LWCON );
                       OPTWRK2 = N + OPTWRK2;
                       OPTWRK = max( OPTWRK, OPTWRK2 );
                   }
                }
             }
         }

         MINWRK = max( 2, MINWRK );
         OPTWRK = max( 2, OPTWRK );
         if ( LCWORK < MINWRK && ( !LQUERY) ) INFO = -19;

      }

      if (INFO == 0 && LRWORK < RMINWRK && !LQUERY) {
         INFO = -21;
      }
      if ( INFO != 0 ) {
         xerbla('ZGESVDQ', -INFO );
         return;
      } else if ( LQUERY ) {

      // Return optimal workspace

          IWORK[1] = IMINWRK;
          CWORK[1] = OPTWRK;
          CWORK[2] = MINWRK;
          RWORK[1] = RMINWRK;
          return;
      }

      // Quick return if the matrix is void.

      if ( ( M == 0 ) || ( N == 0 ) ) {
      // .. all output is void.
         return;
      }

      BIG = DLAMCH('O');
      ASCALED = false;
      if ( ROWPRM ) {
            // .. reordering the rows in decreasing sequence in the
            // ell-infinity norm - this enhances numerical robustness in
            // the case of differently scaled rows.
            for (p = 1; p <= M; p++) { // 1904
                // RWORK(p) = ABS( A(p,IZAMAX(N,A(p,1),LDA)) )
                // [[ZLANGE will return NaN if an entry of the p-th row is Nan]]
                RWORK[p] = ZLANGE( 'M', 1, N, A(p,1), LDA, RDUMMY );
                // .. check for NaN's and Inf's
                if ( ( RWORK(p) != RWORK(p) ) || ( (RWORK(p)*ZERO) != ZERO ) ) {
                    INFO = -8;
                    xerbla('ZGESVDQ', -INFO );
                    return;
                }
            } // 1904
            for (p = 1; p <= M - 1; p++) { // 1952
            q = idamax( M-p+1, RWORK(p), 1 ) + p - 1;
            IWORK[N+p] = q;
            if ( p != q ) {
               RTMP     = RWORK(p);
               RWORK[p] = RWORK(q);
               RWORK[q] = RTMP;
            }
            } // 1952

            if ( RWORK(1) == ZERO ) {
               // Quick return: A is the M x N zero matrix.
               NUMRANK = 0;
               dlaset('G', N, 1, ZERO, ZERO, S, N );
               if (WNTUS) zlaset('G', M, N, CZERO, CONE, U, LDU);
               if (WNTUA) zlaset('G', M, M, CZERO, CONE, U, LDU);
               if (WNTVA) zlaset('G', N, N, CZERO, CONE, V, LDV);
               if ( WNTUF ) {
                   zlaset('G', N, 1, CZERO, CZERO, CWORK, N );
                   zlaset('G', M, N, CZERO, CONE, U, LDU );
               }
               for (p = 1; p <= N; p++) { // 5001
                   IWORK[p] = p;
               } // 5001
               if ( ROWPRM ) {
                   for (p = N + 1; p <= N + M - 1; p++) { // 5002
                       IWORK[p] = p - N;
                   } // 5002
               }
               if (CONDA) RWORK(1) = -1;
               RWORK[2] = -1;
               return;
            }

            if ( RWORK(1) > BIG / sqrt(M.toDouble()) ) {
                // .. to prevent overflow in the QR factorization, scale the
                // matrix by 1/sqrt(M) if too large entry detected
                zlascl('G',0,0,sqrt(M.toDouble()),ONE, M,N, A,LDA, IERR);
                ASCALED = true;
            }
            zlaswp(N, A, LDA, 1, M-1, IWORK(N+1), 1 );
      }

// .. At this stage, preemptive scaling is done only to avoid column
// norms overflows during the QR factorization. The SVD procedure should
// have its own scaling to save the singular values from overflows and
// underflows. That depends on the SVD procedure.

      if ( !ROWPRM ) {
          RTMP = ZLANGE( 'M', M, N, A, LDA, RWORK );
          if ( ( RTMP != RTMP ) || ( (RTMP*ZERO) != ZERO ) ) {
               INFO = -8;
               xerbla('ZGESVDQ', -INFO );
               return;
          }
          if ( RTMP > BIG / sqrt(M.toDouble()) ) {
              // .. to prevent overflow in the QR factorization, scale the
              // matrix by 1/sqrt(M) if too large entry detected
              zlascl('G',0,0, sqrt(M.toDouble()),ONE, M,N, A,LDA, IERR);
              ASCALED = true;
          }
      }

      // .. QR factorization with column pivoting

      // A * P = Q * [ R ]
                  // [ 0 ]

      for (p = 1; p <= N; p++) { // 1963
         // .. all columns are free columns
         IWORK[p] = 0;
      } // 1963
      zgeqp3(M, N, A, LDA, IWORK, CWORK, CWORK(N+1), LCWORK-N, RWORK, IERR );

// If the user requested accuracy level allows truncation in the
// computed upper triangular factor, the matrix R is examined and,
// if possible, replaced with its leading upper trapezoidal part.

      EPSLN = DLAMCH('E');
      SFMIN = DLAMCH('S');
      // SMALL = SFMIN / EPSLN
      NR = N;

      if ( ACCLA ) {

         // Standard absolute error bound suffices. All sigma_i with
         // sigma_i < N*EPS*||A||_F are flushed to zero. This is an
         // aggressive enforcement of lower numerical rank by introducing a
         // backward error of the order of N*EPS*||A||_F.
         NR = 1;
         RTMP = sqrt(N.toDouble())*EPSLN;
         for (p = 2; p <= N; p++) { // 3001
            if ( (A(p,p)).abs() < (RTMP*(A(1,1))) ).abs() GO TO 3002;
               NR = NR + 1;
         } // 3001
         } // 3002

      } else if ( ACCLM ) {
         // .. similarly as above, only slightly more gentle (less aggressive).
         // Sudden drop on the diagonal of R is used as the criterion for being
         // close-to-rank-deficient. The threshold is set to EPSLN=DLAMCH('E').
         // [[This can be made more flexible by replacing this hard-coded value
         // with a user specified threshold.]] Also, the values that underflow
         // will be truncated.
         NR = 1;
         for (p = 2; p <= N; p++) { // 3401
            if ( ( (A(p,p)).abs() < (EPSLN*(A(p-1,p-1))) ).abs() || ( (A(p,p)) < SFMIN ) ).abs() GO TO 3402;
            NR = NR + 1;
         } // 3401
         } // 3402

      } else {
         // .. RRQR not authorized to determine numerical rank except in the
         // obvious case of zero pivots.
         // .. inspect R for exact zeros on the diagonal;
         // R(i,i)=0 => R(i:N,i:N)=0.
         NR = 1;
         for (p = 2; p <= N; p++) { // 3501
            if ( (A(p,p)).abs() == ZERO ) GO TO 3502;
            NR = NR + 1;
         } // 3501
         } // 3502

         if ( CONDA ) {
            // Estimate the scaled condition number of A. Use the fact that it is
            // the same as the scaled condition number of R.
               // .. V is used as workspace
               zlacpy('U', N, N, A, LDA, V, LDV );
               // Only the leading NR x NR submatrix of the triangular factor
               // is considered. Only if NR=N will this give a reliable error
               // bound. However, even for NR < N, this can be used on an
               // expert level and obtain useful information in the sense of
               // perturbation theory.
               for (p = 1; p <= NR; p++) { // 3053
                  RTMP = DZNRM2( p, V(1,p), 1 );
                  zdscal(p, ONE/RTMP, V(1,p), 1 );
               } // 3053
               if ( !( LSVEC || RSVEC ) ) {
                   zpocon('U', NR, V, LDV, ONE, RTMP, CWORK, RWORK, IERR );
               } else {
                   zpocon('U', NR, V, LDV, ONE, RTMP, CWORK(N+1), RWORK, IERR );
               }
               SCONDA = ONE / sqrt(RTMP);
            // For NR=N, SCONDA is an estimate of sqrt(||(R^* * R)^(-1)||_1),
            // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
            // See the reference [1] for more details.
         }

      }

      if ( WNTUR ) {
          N1 = NR;
      } else if ( WNTUS || WNTUF) {
          N1 = N;
      } else if ( WNTUA ) {
          N1 = M;
      }

      if ( !( RSVEC || LSVEC ) ) {
// .......................................................................
         // .. only the singular values are requested
// .......................................................................
         if ( RTRANS ) {

          // .. compute the singular values of R**H = [A](1:NR,1:N)**H
            // .. set the lower triangle of [A] to [A](1:NR,1:N)**H and
            // the upper triangle of [A] to zero.
            for (p = 1; p <= min( N, NR ); p++) { // 1146
               A[p,p] = CONJG(A(p,p));
               for (q = p + 1; q <= N; q++) { // 1147
                  A[q,p] = CONJG(A(p,q));
                  if (q <= NR) A(p,q) = CZERO;
               } // 1147
            } // 1146

            zgesvd('N', 'N', N, NR, A, LDA, S, U, LDU, V, LDV, CWORK, LCWORK, RWORK, INFO );

         } else {

            // .. compute the singular values of R = [A](1:NR,1:N)

            if (NR > 1) zlaset( 'L', NR-1,NR-1, CZERO,CZERO, A(2,1), LDA );
            zgesvd('N', 'N', NR, N, A, LDA, S, U, LDU, V, LDV, CWORK, LCWORK, RWORK, INFO );

         }

      } else if ( LSVEC && ( !RSVEC) ) {
// .......................................................................
        // .. the singular values and the left singular vectors requested
// .......................................................................""""""""
         if ( RTRANS ) {
             // .. apply ZGESVD to R**H
             // .. copy R**H into [U] and overwrite [U] with the right singular
             // vectors of R
            for (p = 1; p <= NR; p++) { // 1192
               for (q = p; q <= N; q++) { // 1193
                  U[q,p] = CONJG(A(p,q));
               } // 1193
            } // 1192
            if (NR > 1) zlaset( 'U', NR-1,NR-1, CZERO,CZERO, U(1,2), LDU );
            // .. the left singular vectors not computed, the NR right singular
            // vectors overwrite [U](1:NR,1:NR) as conjugate transposed. These
            // will be pre-multiplied by Q to build the left singular vectors of A.
               zgesvd('N', 'O', N, NR, U, LDU, S, U, LDU, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO );

               for (p = 1; p <= NR; p++) { // 1119
                   U[p,p] = CONJG(U(p,p));
                   for (q = p + 1; q <= NR; q++) { // 1120
                      CTMP   = CONJG(U(q,p));
                      U[q,p] = CONJG(U(p,q));
                      U[p,q] = CTMP;
                   } // 1120
               } // 1119

         } else {
             // .. apply ZGESVD to R
             // .. copy R into [U] and overwrite [U] with the left singular vectors
             zlacpy('U', NR, N, A, LDA, U, LDU );
             if (NR > 1) zlaset( 'L', NR-1, NR-1, CZERO, CZERO, U(2,1), LDU );
             // .. the right singular vectors not computed, the NR left singular
             // vectors overwrite [U](1:NR,1:NR)
                zgesvd('O', 'N', NR, N, U, LDU, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO );
                // .. now [U](1:NR,1:NR) contains the NR left singular vectors of
                // R. These will be pre-multiplied by Q to build the left singular
                // vectors of A.
         }

            // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
         if ( ( NR < M ) && ( !WNTUF ) ) {
             zlaset('A', M-NR, NR, CZERO, CZERO, U(NR+1,1), LDU);
             if ( NR < N1 ) {
                zlaset('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1), LDU );
                zlaset('A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU );
             }
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         if ( !WNTUF) zunmqr( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LCWORK-N, IERR );
         if (ROWPRM && !WNTUF) zlaswp( N1, U, LDU, 1, M-1, IWORK(N+1), -1 );

      } else if ( RSVEC && ( !LSVEC ) ) {
// .......................................................................
        // .. the singular values and the right singular vectors requested
// .......................................................................
          if ( RTRANS ) {
             // .. apply ZGESVD to R**H
             // .. copy R**H into V and overwrite V with the left singular vectors
            for (p = 1; p <= NR; p++) { // 1165
               for (q = p; q <= N; q++) { // 1166
                  V[q,p] = CONJG(A(p,q));
               } // 1166
            } // 1165
            if (NR > 1) zlaset( 'U', NR-1,NR-1, CZERO,CZERO, V(1,2), LDV );
            // .. the left singular vectors of R**H overwrite V, the right singular
            // vectors not computed
            if ( WNTVR || ( NR == N ) ) {
               zgesvd('O', 'N', N, NR, V, LDV, S, U, LDU, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO );

               for (p = 1; p <= NR; p++) { // 1121
                   V[p,p] = CONJG(V(p,p));
                   for (q = p + 1; q <= NR; q++) { // 1122
                      CTMP   = CONJG(V(q,p));
                      V[q,p] = CONJG(V(p,q));
                      V[p,q] = CTMP;
                   } // 1122
               } // 1121

               if ( NR < N ) {
                   for (p = 1; p <= NR; p++) { // 1103
                      for (q = NR + 1; q <= N; q++) { // 1104
                          V[p,q] = CONJG(V(q,p));
                      } // 1104
                   } // 1103
               }
               zlapmt( false , NR, N, V, LDV, IWORK );
            } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:N,1:NR)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the QR factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                zlaset('G', N, N-NR, CZERO, CZERO, V(1,NR+1), LDV);
                zgesvd('O', 'N', N, N, V, LDV, S, U, LDU, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO );

                for (p = 1; p <= N; p++) { // 1123
                   V[p,p] = CONJG(V(p,p));
                   for (q = p + 1; q <= N; q++) { // 1124
                      CTMP   = CONJG(V(q,p));
                      V[q,p] = CONJG(V(p,q));
                      V[p,q] = CTMP;
                   } // 1124
                } // 1123
                zlapmt( false , N, N, V, LDV, IWORK );
            }

          } else {
             // .. aply ZGESVD to R
             // .. copy R into V and overwrite V with the right singular vectors
             zlacpy('U', NR, N, A, LDA, V, LDV );
             if (NR > 1) zlaset( 'L', NR-1, NR-1, CZERO, CZERO, V(2,1), LDV );
             // .. the right singular vectors overwrite V, the NR left singular
             // vectors stored in U(1:NR,1:NR)
             if ( WNTVR || ( NR == N ) ) {
                zgesvd('N', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO );
                zlapmt( false , NR, N, V, LDV, IWORK );
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
             } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:NR,1:N)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the LQ factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                 zlaset('G', N-NR, N, CZERO,CZERO, V(NR+1,1), LDV);
                 zgesvd('N', 'O', N, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO );
                 zlapmt( false , N, N, V, LDV, IWORK );
             }
             // .. now [V] contains the adjoint of the matrix of the right singular
             // vectors of A.
          }

      } else {
// .......................................................................
        // .. FULL SVD requested
// .......................................................................
         if ( RTRANS ) {

             // .. apply ZGESVD to R**H [[this option is left for R&D&T]]

            if ( WNTVR || ( NR == N ) ) {
             // .. copy R**H into [V] and overwrite [V] with the left singular
             // vectors of R**H
            for (p = 1; p <= NR; p++) { // 1168
               for (q = p; q <= N; q++) { // 1169
                  V[q,p] = CONJG(A(p,q));
               } // 1169
            } // 1168
            if (NR > 1) zlaset( 'U', NR-1,NR-1, CZERO,CZERO, V(1,2), LDV );

            // .. the left singular vectors of R**H overwrite [V], the NR right
            // singular vectors of R**H stored in [U](1:NR,1:NR) as conjugate
            // transposed
               zgesvd('O', 'A', N, NR, V, LDV, S, V, LDV, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO );
               // .. assemble V
               for (p = 1; p <= NR; p++) { // 1115
                  V[p,p] = CONJG(V(p,p));
                  for (q = p + 1; q <= NR; q++) { // 1116
                     CTMP   = CONJG(V(q,p));
                     V[q,p] = CONJG(V(p,q));
                     V[p,q] = CTMP;
                  } // 1116
               } // 1115
               if ( NR < N ) {
                   for (p = 1; p <= NR; p++) { // 1101
                      for (q = NR+1; q <= N; q++) { // 1102
                         V[p,q] = CONJG(V(q,p));
                      } // 1102
                   } // 1101
               }
               zlapmt( false , NR, N, V, LDV, IWORK );

                for (p = 1; p <= NR; p++) { // 1117
                   U[p,p] = CONJG(U(p,p));
                   for (q = p + 1; q <= NR; q++) { // 1118
                      CTMP   = CONJG(U(q,p));
                      U[q,p] = CONJG(U(p,q));
                      U[p,q] = CTMP;
                   } // 1118
                } // 1117

                if ( ( NR < M ) && !(WNTUF)) {
                  zlaset('A', M-NR,NR, CZERO,CZERO, U(NR+1,1), LDU);
                  if ( NR < N1 ) {
                     zlaset('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU);
                     zlaset('A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU );
                  }
               }

            } else {
                // .. need all N right singular vectors and NR < N
             // .. copy R**H into [V] and overwrite [V] with the left singular
             // vectors of R**H
                // [[The optimal ratio N/NR for using QRF instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'ZGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = max( OPTRATIO, 2 )
                OPTRATIO = 2;
                if ( OPTRATIO*NR > N ) {
                   for (p = 1; p <= NR; p++) { // 1198
                      for (q = p; q <= N; q++) { // 1199
                         V[q,p] = CONJG(A(p,q));
                      } // 1199
                   } // 1198
                   if (NR > 1) zlaset('U',NR-1,NR-1, CZERO,CZERO, V(1,2),LDV);

                   zlaset('A',N,N-NR,CZERO,CZERO,V(1,NR+1),LDV);
                   zgesvd('O', 'A', N, N, V, LDV, S, V, LDV, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO );

                   for (p = 1; p <= N; p++) { // 1113
                      V[p,p] = CONJG(V(p,p));
                      for (q = p + 1; q <= N; q++) { // 1114
                         CTMP   = CONJG(V(q,p));
                         V[q,p] = CONJG(V(p,q));
                         V[p,q] = CTMP;
                      } // 1114
                   } // 1113
                   zlapmt( false , N, N, V, LDV, IWORK );
               // .. assemble the left singular vector matrix U of dimensions
               // (M x N1), i.e. (M x N) or (M x M).

                   for (p = 1; p <= N; p++) { // 1111
                      U[p,p] = CONJG(U(p,p));
                      for (q = p + 1; q <= N; q++) { // 1112
                         CTMP   = CONJG(U(q,p));
                         U[q,p] = CONJG(U(p,q));
                         U[p,q] = CTMP;
                      } // 1112
                   } // 1111

                   if ( ( N < M ) && !(WNTUF)) {
                      zlaset('A',M-N,N,CZERO,CZERO,U(N+1,1),LDU);
                      if ( N < N1 ) {
                        zlaset('A',N,N1-N,CZERO,CZERO,U(1,N+1),LDU);
                        zlaset('A',M-N,N1-N,CZERO,CONE, U(N+1,N+1), LDU );
                      }
                   }
                } else {
                   // .. copy R**H into [U] and overwrite [U] with the right
                   // singular vectors of R
                   for (p = 1; p <= NR; p++) { // 1196
                      for (q = p; q <= N; q++) { // 1197
                         U[q,NR+p] = CONJG(A(p,q));
                      } // 1197
                   } // 1196
                   if (NR > 1) zlaset('U',NR-1,NR-1,CZERO,CZERO,U(1,NR+2),LDU);
                   zgeqrf(N, NR, U(1,NR+1), LDU, CWORK(N+1), CWORK(N+NR+1), LCWORK-N-NR, IERR );
                   for (p = 1; p <= NR; p++) { // 1143
                       for (q = 1; q <= N; q++) { // 1144
                           V[q,p] = CONJG(U(p,NR+q));
                       } // 1144
                   } // 1143
                  zlaset('U',NR-1,NR-1,CZERO,CZERO,V(1,2),LDV);
                  zgesvd('S', 'O', NR, NR, V, LDV, S, U, LDU, V,LDV, CWORK(N+NR+1),LCWORK-N-NR,RWORK, INFO );
                  zlaset('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV);
                  zlaset('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV);
                  zlaset('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV);
                  zunmqr('R','C', N, N, NR, U(1,NR+1), LDU, CWORK(N+1),V,LDV,CWORK(N+NR+1),LCWORK-N-NR,IERR);
                  zlapmt( false , N, N, V, LDV, IWORK );
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x NR) or (M x N) or (M x M).
                  if ( ( NR < M ) && !(WNTUF)) {
                     zlaset('A',M-NR,NR,CZERO,CZERO,U(NR+1,1),LDU);
                     if ( NR < N1 ) {
                     zlaset('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU);
                     zlaset('A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1),LDU);
                     }
                  }
                }
            }

         } else {

             // .. apply ZGESVD to R [[this is the recommended option]]

             if ( WNTVR || ( NR == N ) ) {
                 // .. copy R into [V] and overwrite V with the right singular vectors
                 zlacpy('U', NR, N, A, LDA, V, LDV );
                if (NR > 1) zlaset( 'L', NR-1,NR-1, CZERO,CZERO, V(2,1), LDV );
                // .. the right singular vectors of R overwrite [V], the NR left
                // singular vectors of R stored in [U](1:NR,1:NR)
                zgesvd('S', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO );
                zlapmt( false , NR, N, V, LDV, IWORK );
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
               if ( ( NR < M ) && !(WNTUF)) {
                  zlaset('A', M-NR,NR, CZERO,CZERO, U(NR+1,1), LDU);
                  if ( NR < N1 ) {
                     zlaset('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU);
                     zlaset('A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU );
                  }
               }

             } else {
               // .. need all N right singular vectors and NR < N
               // .. the requested number of the left singular vectors
                // is then N1 (N or M)
                // [[The optimal ratio N/NR for using LQ instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'ZGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = max( OPTRATIO, 2 )
               OPTRATIO = 2;
               if ( OPTRATIO * NR > N ) {
                  zlacpy('U', NR, N, A, LDA, V, LDV );
                  if (NR > 1) zlaset('L', NR-1,NR-1, CZERO,CZERO, V(2,1),LDV);
               // .. the right singular vectors of R overwrite [V], the NR left
                  // singular vectors of R stored in [U](1:NR,1:NR)
                  zlaset('A', N-NR,N, CZERO,CZERO, V(NR+1,1),LDV);
                  zgesvd('S', 'O', N, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO );
                  zlapmt( false , N, N, V, LDV, IWORK );
                  // .. now [V] contains the adjoint of the matrix of the right
                  // singular vectors of A. The leading N left singular vectors
                  // are in [U](1:N,1:N)
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x N1), i.e. (M x N) or (M x M).
                  if ( ( N < M ) && !(WNTUF)) {
                      zlaset('A',M-N,N,CZERO,CZERO,U(N+1,1),LDU);
                      if ( N < N1 ) {
                        zlaset('A',N,N1-N,CZERO,CZERO,U(1,N+1),LDU);
                        zlaset('A',M-N,N1-N,CZERO,CONE, U(N+1,N+1), LDU );
                      }
                  }
               } else {
                  zlacpy('U', NR, N, A, LDA, U(NR+1,1), LDU );
                  if (NR > 1) zlaset('L',NR-1,NR-1,CZERO,CZERO,U(NR+2,1),LDU);
                  zgelqf(NR, N, U(NR+1,1), LDU, CWORK(N+1), CWORK(N+NR+1), LCWORK-N-NR, IERR );
                  zlacpy('L',NR,NR,U(NR+1,1),LDU,V,LDV);
                  if (NR > 1) zlaset('U',NR-1,NR-1,CZERO,CZERO,V(1,2),LDV);
                  zgesvd('S', 'O', NR, NR, V, LDV, S, U, LDU, V, LDV, CWORK(N+NR+1), LCWORK-N-NR, RWORK, INFO );
                  zlaset('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV);
                  zlaset('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV);
                  zlaset('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV);
                  zunmlq('R','N',N,N,NR,U(NR+1,1),LDU,CWORK(N+1), V, LDV, CWORK(N+NR+1),LCWORK-N-NR,IERR);
                  zlapmt( false , N, N, V, LDV, IWORK );
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
                  if ( ( NR < M ) && !(WNTUF)) {
                     zlaset('A',M-NR,NR,CZERO,CZERO,U(NR+1,1),LDU);
                     if ( NR < N1 ) {
                     zlaset('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU);
                     zlaset('A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU );
                     }
                  }
               }
             }
         // .. end of the "R**H or R" branch
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         if ( !WNTUF) zunmqr( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LCWORK-N, IERR );
         if (ROWPRM && !WNTUF) zlaswp( N1, U, LDU, 1, M-1, IWORK(N+1), -1 );

      // ... end of the "full SVD" branch
      }

      // Check whether some singular values are returned as zeros, e.g.
      // due to underflow, and update the numerical rank.
      p = NR;
      for (q = p; q >= 1; q--) { // 4001
          if ( S(q) > ZERO ) GO TO 4002;
          NR = NR - 1;
      } // 4001
      } // 4002

      // .. if numerical rank deficiency is detected, the truncated
      // singular values are set to zero.
      if (NR < N) dlaset( 'G', N-NR,1, ZERO,ZERO, S(NR+1), N );
      // .. undo scaling; this may cause overflow in the largest singular
      // values.
      if (ASCALED) dlascl( 'G',0,0, ONE,sqrt(M.toDouble()), NR,1, S, N, IERR );
      if (CONDA) RWORK(1) = SCONDA;
      RWORK[2] = p - NR;
      // .. p-NR is the number of singular values that are computed as
      // exact zeros in ZGESVD() applied to the (possibly truncated)
      // full row rank triangular (trapezoidal) factor of A.
      NUMRANK = NR;

      return;
      }