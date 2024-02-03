      SUBROUTINE CGESVDQ( JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA, S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK, CWORK, LCWORK, RWORK, LRWORK, INFO )
      // .. Scalar Arguments ..
      IMPLICIT    NONE
      String      JOBA, JOBP, JOBR, JOBU, JOBV;
      int         M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LCWORK, LRWORK, INFO;
      // ..
      // .. Array Arguments ..
      COMPLEX     A( LDA, * ), U( LDU, * ), V( LDV, * ), CWORK( * )
      REAL        S( * ), RWORK( * )
      int         IWORK( * );

*  =====================================================================

      // .. Parameters ..
      REAL        ZERO,         ONE
      const     ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX     CZERO,                    CONE
      const     CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      int         IERR, NR, N1, OPTRATIO, p, q;
      int         LWCON, LWQP3, LWRK_CGELQF, LWRK_CGESVD, LWRK_CGESVD2, LWRK_CGEQP3, LWRK_CGEQRF, LWRK_CUNMLQ, LWRK_CUNMQR, LWRK_CUNMQR2, LWLQF, LWQRF, LWSVD, LWSVD2, LWUNQ, LWUNQ2, LWUNLQ, MINWRK, MINWRK2, OPTWRK, OPTWRK2, IMINWRK, RMINWRK;
      bool        ACCLA,  ACCLM, ACCLH, ASCALED, CONDA, DNTWU,  DNTWV, LQUERY, LSVC0, LSVEC, ROWPRM,  RSVEC, RTRANS, WNTUA, WNTUF,  WNTUR, WNTUS, WNTVA,   WNTVR;
      REAL        BIG, EPSLN, RTMP, SCONDA, SFMIN
      COMPLEX     CTMP
      // ..
      // .. Local Arrays
      COMPLEX     CDUMMY(1)
      REAL        RDUMMY(1)
      // ..
      // .. External Subroutines (BLAS, LAPACK)
      // EXTERNAL CGELQF, CGEQP3, CGEQRF, CGESVD, CLACPY, CLAPMT, CLASCL, CLASET, CLASWP, CSSCAL, SLASET, SLASCL, CPOCON, CUNMLQ, CUNMQR, XERBLA
      // ..
      // .. External Functions (BLAS, LAPACK)
      bool       LSAME;
      int        ISAMAX;
      REAL       CLANGE, SCNRM2, SLAMCH
      // EXTERNAL CLANGE, LSAME, ISAMAX, SCNRM2, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

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
         IMINWRK = MAX( 1, N + M - 1 )
         RMINWRK = MAX( 2, M, 5*N )
      } else {
         IMINWRK = MAX( 1, N )
         RMINWRK = MAX( 2, 5*N )
      }
      LQUERY = (LIWORK .EQ. -1 .OR. LCWORK .EQ. -1 .OR. LRWORK .EQ. -1)
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

      // Compute workspace
         // .. compute the minimal and the optimal workspace lengths
         // [[The expressions for computing the minimal and the optimal
         // values of LCWORK are written with a lot of redundancy and
         // can be simplified. However, this detailed form is easier for
         // maintenance and modifications of the code.]]

         // .. minimal workspace length for CGEQP3 of an M x N matrix
         LWQP3 = N+1
         // .. minimal workspace length for CUNMQR to build left singular vectors
         if ( WNTUS .OR. WNTUR ) {
             LWUNQ  = MAX( N  , 1 )
         } else if ( WNTUA ) {
             LWUNQ = MAX( M , 1 )
         }
         // .. minimal workspace length for CPOCON of an N x N matrix
         LWCON = 2 * N
         // .. CGESVD of an N x N matrix
         LWSVD = MAX( 3 * N, 1 )
         if ( LQUERY ) {
             CALL CGEQP3( M, N, A, LDA, IWORK, CDUMMY, CDUMMY, -1, RDUMMY, IERR )
             LWRK_CGEQP3 = INT( CDUMMY(1) )
             if ( WNTUS .OR. WNTUR ) {
                 CALL CUNMQR( 'L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR )
                 LWRK_CUNMQR = INT( CDUMMY(1) )
             } else if ( WNTUA ) {
                 CALL CUNMQR( 'L', 'N', M, M, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR )
                 LWRK_CUNMQR = INT( CDUMMY(1) )
             } else {
                 LWRK_CUNMQR = 0
             }
         }
         MINWRK = 2
         OPTWRK = 2
         if ( .NOT. (LSVEC .OR. RSVEC )) {
             // .. minimal and optimal sizes of the complex workspace if
             // only the singular values are requested
             if ( CONDA ) {
                MINWRK = MAX( N+LWQP3, LWCON, LWSVD )
             } else {
                MINWRK = MAX( N+LWQP3, LWSVD )
             }
             if ( LQUERY ) {
                 CALL CGESVD( 'N', 'N', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                 LWRK_CGESVD = INT( CDUMMY(1) )
                 if ( CONDA ) {
                    OPTWRK = MAX( N+LWRK_CGEQP3, N+LWCON, LWRK_CGESVD )
                 } else {
                    OPTWRK = MAX( N+LWRK_CGEQP3, LWRK_CGESVD )
                 }
             }
         } else if ( LSVEC .AND. (.NOT.RSVEC) ) {
             // .. minimal and optimal sizes of the complex workspace if the
             // singular values and the left singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD, LWUNQ )
             } else {
                 MINWRK = N + MAX( LWQP3, LWSVD, LWUNQ )
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   CALL CGESVD( 'N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                } else {
                   CALL CGESVD( 'O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                }
                LWRK_CGESVD = INT( CDUMMY(1) )
                if ( CONDA ) {
                    OPTWRK = N + MAX( LWRK_CGEQP3, LWCON, LWRK_CGESVD, LWRK_CUNMQR )
                } else {
                    OPTWRK = N + MAX( LWRK_CGEQP3, LWRK_CGESVD, LWRK_CUNMQR )
                }
             }
         } else if ( RSVEC .AND. (.NOT.LSVEC) ) {
             // .. minimal and optimal sizes of the complex workspace if the
             // singular values and the right singular vectors are requested
             if ( CONDA ) {
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD )
             } else {
                 MINWRK = N + MAX( LWQP3, LWSVD )
             }
             if ( LQUERY ) {
                 if ( RTRANS ) {
                     CALL CGESVD( 'O', 'N', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                 } else {
                     CALL CGESVD( 'N', 'O', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                 }
                 LWRK_CGESVD = INT( CDUMMY(1) )
                 if ( CONDA ) {
                     OPTWRK = N + MAX( LWRK_CGEQP3, LWCON, LWRK_CGESVD )
                 } else {
                     OPTWRK = N + MAX( LWRK_CGEQP3, LWRK_CGESVD )
                 }
             }
         } else {
             // .. minimal and optimal sizes of the complex workspace if the
             // full SVD is requested
             if ( RTRANS ) {
                 MINWRK = MAX( LWQP3, LWSVD, LWUNQ )
                 IF ( CONDA ) MINWRK = MAX( MINWRK, LWCON )
                 MINWRK = MINWRK + N
                 if ( WNTVA ) {
                    // .. minimal workspace length for N x N/2 CGEQRF
                    LWQRF  = MAX( N/2, 1 )
                    // .. minimal workspace length for N/2 x N/2 CGESVD
                    LWSVD2 = MAX( 3 * (N/2), 1 )
                    LWUNQ2 = MAX( N, 1 )
                    MINWRK2 = MAX( LWQP3, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ )
                    IF ( CONDA ) MINWRK2 = MAX( MINWRK2, LWCON )
                    MINWRK2 = N + MINWRK2
                    MINWRK = MAX( MINWRK, MINWRK2 )
                 }
             } else {
                 MINWRK = MAX( LWQP3, LWSVD, LWUNQ )
                 IF ( CONDA ) MINWRK = MAX( MINWRK, LWCON )
                 MINWRK = MINWRK + N
                 if ( WNTVA ) {
                    // .. minimal workspace length for N/2 x N CGELQF
                    LWLQF  = MAX( N/2, 1 )
                    LWSVD2 = MAX( 3 * (N/2), 1 )
                    LWUNLQ = MAX( N , 1 )
                    MINWRK2 = MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ )
                    IF ( CONDA ) MINWRK2 = MAX( MINWRK2, LWCON )
                    MINWRK2 = N + MINWRK2
                    MINWRK = MAX( MINWRK, MINWRK2 )
                 }
             }
             if ( LQUERY ) {
                if ( RTRANS ) {
                   CALL CGESVD( 'O', 'A', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                   LWRK_CGESVD = INT( CDUMMY(1) )
                   OPTWRK = MAX(LWRK_CGEQP3,LWRK_CGESVD,LWRK_CUNMQR)
                   IF ( CONDA ) OPTWRK = MAX( OPTWRK, LWCON )
                   OPTWRK = N + OPTWRK
                   if ( WNTVA ) {
                       CALL CGEQRF(N,N/2,U,LDU,CDUMMY,CDUMMY,-1,IERR)
                       LWRK_CGEQRF = INT( CDUMMY(1) )
                       CALL CGESVD( 'S', 'O', N/2,N/2, V,LDV, S, U,LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                       LWRK_CGESVD2 = INT( CDUMMY(1) )
                       CALL CUNMQR( 'R', 'C', N, N, N/2, U, LDU, CDUMMY, V, LDV, CDUMMY, -1, IERR )
                       LWRK_CUNMQR2 = INT( CDUMMY(1) )
                       OPTWRK2 = MAX( LWRK_CGEQP3, N/2+LWRK_CGEQRF, N/2+LWRK_CGESVD2, N/2+LWRK_CUNMQR2 )
                       IF ( CONDA ) OPTWRK2 = MAX( OPTWRK2, LWCON )
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   }
                } else {
                   CALL CGESVD( 'S', 'O', N, N, A, LDA, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                   LWRK_CGESVD = INT( CDUMMY(1) )
                   OPTWRK = MAX(LWRK_CGEQP3,LWRK_CGESVD,LWRK_CUNMQR)
                   IF ( CONDA ) OPTWRK = MAX( OPTWRK, LWCON )
                   OPTWRK = N + OPTWRK
                   if ( WNTVA ) {
                      CALL CGELQF(N/2,N,U,LDU,CDUMMY,CDUMMY,-1,IERR)
                      LWRK_CGELQF = INT( CDUMMY(1) )
                      CALL CGESVD( 'S','O', N/2,N/2, V, LDV, S, U, LDU, V, LDV, CDUMMY, -1, RDUMMY, IERR )
                      LWRK_CGESVD2 = INT( CDUMMY(1) )
                      CALL CUNMLQ( 'R', 'N', N, N, N/2, U, LDU, CDUMMY, V, LDV, CDUMMY,-1,IERR )
                      LWRK_CUNMLQ = INT( CDUMMY(1) )
                      OPTWRK2 = MAX( LWRK_CGEQP3, N/2+LWRK_CGELQF, N/2+LWRK_CGESVD2, N/2+LWRK_CUNMLQ )
                       IF ( CONDA ) OPTWRK2 = MAX( OPTWRK2, LWCON )
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   }
                }
             }
         }

         MINWRK = MAX( 2, MINWRK )
         OPTWRK = MAX( 2, OPTWRK )
         IF ( LCWORK .LT. MINWRK .AND. (.NOT.LQUERY) ) INFO = -19

      }

      if (INFO .EQ. 0 .AND. LRWORK .LT. RMINWRK .AND. .NOT. LQUERY) {
         INFO = -21
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CGESVDQ', -INFO )
         RETURN
      } else if ( LQUERY ) {

      // Return optimal workspace

          IWORK(1) = IMINWRK
          CWORK(1) = OPTWRK
          CWORK(2) = MINWRK
          RWORK(1) = RMINWRK
          RETURN
      }

      // Quick return if the matrix is void.

      if ( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) {
      // .. all output is void.
         RETURN
      }

      BIG = SLAMCH('O')
      ASCALED = .FALSE.
      if ( ROWPRM ) {
            // .. reordering the rows in decreasing sequence in the
            // ell-infinity norm - this enhances numerical robustness in
            // the case of differently scaled rows.
            DO 1904 p = 1, M
                // RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
                // [[CLANGE will return NaN if an entry of the p-th row is Nan]]
                RWORK(p) = CLANGE( 'M', 1, N, A(p,1), LDA, RDUMMY )
                // .. check for NaN's and Inf's
                if ( ( RWORK(p) .NE. RWORK(p) ) .OR. ( (RWORK(p)*ZERO) .NE. ZERO ) ) {
                    INFO = - 8
                    CALL XERBLA( 'CGESVDQ', -INFO )
                    RETURN
                }
 1904       CONTINUE
            DO 1952 p = 1, M - 1
            q = ISAMAX( M-p+1, RWORK(p), 1 ) + p - 1
            IWORK(N+p) = q
            if ( p .NE. q ) {
               RTMP     = RWORK(p)
               RWORK(p) = RWORK(q)
               RWORK(q) = RTMP
            }
 1952       CONTINUE

            if ( RWORK(1) .EQ. ZERO ) {
               // Quick return: A is the M x N zero matrix.
               NUMRANK = 0
               CALL SLASET( 'G', N, 1, ZERO, ZERO, S, N )
               IF ( WNTUS ) CALL CLASET('G', M, N, CZERO, CONE, U, LDU)
               IF ( WNTUA ) CALL CLASET('G', M, M, CZERO, CONE, U, LDU)
               IF ( WNTVA ) CALL CLASET('G', N, N, CZERO, CONE, V, LDV)
               if ( WNTUF ) {
                   CALL CLASET( 'G', N, 1, CZERO, CZERO, CWORK, N )
                   CALL CLASET( 'G', M, N, CZERO, CONE, U, LDU )
               }
               DO 5001 p = 1, N
                   IWORK(p) = p
 5001          CONTINUE
               if ( ROWPRM ) {
                   DO 5002 p = N + 1, N + M - 1
                       IWORK(p) = p - N
 5002              CONTINUE
               }
               IF ( CONDA ) RWORK(1) = -1
               RWORK(2) = -1
               RETURN
            }

            if ( RWORK(1) .GT. BIG / SQRT(REAL(M)) ) {
                // .. to prevent overflow in the QR factorization, scale the
                // matrix by 1/sqrt(M) if too large entry detected
                CALL CLASCL('G',0,0,SQRT(REAL(M)),ONE, M,N, A,LDA, IERR)
                ASCALED = .TRUE.
            }
            CALL CLASWP( N, A, LDA, 1, M-1, IWORK(N+1), 1 )
      }

*    .. At this stage, preemptive scaling is done only to avoid column
*    norms overflows during the QR factorization. The SVD procedure should
*    have its own scaling to save the singular values from overflows and
*    underflows. That depends on the SVD procedure.

      if ( .NOT.ROWPRM ) {
          RTMP = CLANGE( 'M', M, N, A, LDA, RWORK )
          if ( ( RTMP .NE. RTMP ) .OR. ( (RTMP*ZERO) .NE. ZERO ) ) {
               INFO = - 8
               CALL XERBLA( 'CGESVDQ', -INFO )
               RETURN
          }
          if ( RTMP .GT. BIG / SQRT(REAL(M)) ) {
              // .. to prevent overflow in the QR factorization, scale the
              // matrix by 1/sqrt(M) if too large entry detected
              CALL CLASCL('G',0,0, SQRT(REAL(M)),ONE, M,N, A,LDA, IERR)
              ASCALED = .TRUE.
          }
      }

      // .. QR factorization with column pivoting

      // A * P = Q * [ R ]
                  // [ 0 ]

      DO 1963 p = 1, N
         // .. all columns are free columns
         IWORK(p) = 0
 1963 CONTINUE
      CALL CGEQP3( M, N, A, LDA, IWORK, CWORK, CWORK(N+1), LCWORK-N, RWORK, IERR )

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
         DO 3001 p = 2, N
            IF ( ABS(A(p,p)) .LT. (RTMP*ABS(A(1,1))) ) GO TO 3002
               NR = NR + 1
 3001    CONTINUE
 3002    CONTINUE

      ELSEIF ( ACCLM ) THEN
         // .. similarly as above, only slightly more gentle (less aggressive).
         // Sudden drop on the diagonal of R is used as the criterion for being
         // close-to-rank-deficient. The threshold is set to EPSLN=SLAMCH('E').
         // [[This can be made more flexible by replacing this hard-coded value
         // with a user specified threshold.]] Also, the values that underflow
         // will be truncated.
         NR = 1
         DO 3401 p = 2, N
            IF ( ( ABS(A(p,p)) .LT. (EPSLN*ABS(A(p-1,p-1))) ) .OR. ( ABS(A(p,p)) .LT. SFMIN ) ) GO TO 3402
            NR = NR + 1
 3401    CONTINUE
 3402    CONTINUE

      } else {
         // .. RRQR not authorized to determine numerical rank except in the
         // obvious case of zero pivots.
         // .. inspect R for exact zeros on the diagonal;
         // R(i,i)=0 => R(i:N,i:N)=0.
         NR = 1
         DO 3501 p = 2, N
            IF ( ABS(A(p,p)) .EQ. ZERO ) GO TO 3502
            NR = NR + 1
 3501    CONTINUE
 3502    CONTINUE

         if ( CONDA ) {
            // Estimate the scaled condition number of A. Use the fact that it is
            // the same as the scaled condition number of R.
               // .. V is used as workspace
               CALL CLACPY( 'U', N, N, A, LDA, V, LDV )
               // Only the leading NR x NR submatrix of the triangular factor
               // is considered. Only if NR=N will this give a reliable error
               // bound. However, even for NR < N, this can be used on an
               // expert level and obtain useful information in the sense of
               // perturbation theory.
               DO 3053 p = 1, NR
                  RTMP = SCNRM2( p, V(1,p), 1 )
                  CALL CSSCAL( p, ONE/RTMP, V(1,p), 1 )
 3053          CONTINUE
               if ( .NOT. ( LSVEC .OR. RSVEC ) ) {
                   CALL CPOCON( 'U', NR, V, LDV, ONE, RTMP, CWORK, RWORK, IERR )
               } else {
                   CALL CPOCON( 'U', NR, V, LDV, ONE, RTMP, CWORK(N+1), RWORK, IERR )
               }
               SCONDA = ONE / SQRT(RTMP)
            // For NR=N, SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1),
            // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
            // See the reference [1] for more details.
         }

      ENDIF

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

          // .. compute the singular values of R**H = [A](1:NR,1:N)**H
            // .. set the lower triangle of [A] to [A](1:NR,1:N)**H and
            // the upper triangle of [A] to zero.
            DO 1146 p = 1, MIN( N, NR )
               A(p,p) = CONJG(A(p,p))
               DO 1147 q = p + 1, N
                  A(q,p) = CONJG(A(p,q))
                  IF ( q .LE. NR ) A(p,q) = CZERO
 1147          CONTINUE
 1146       CONTINUE

            CALL CGESVD( 'N', 'N', N, NR, A, LDA, S, U, LDU, V, LDV, CWORK, LCWORK, RWORK, INFO )

         } else {

            // .. compute the singular values of R = [A](1:NR,1:N)

            IF ( NR .GT. 1 ) CALL CLASET( 'L', NR-1,NR-1, CZERO,CZERO, A(2,1), LDA )             CALL CGESVD( 'N', 'N', NR, N, A, LDA, S, U, LDU, V, LDV, CWORK, LCWORK, RWORK, INFO )

         }

      } else if ( LSVEC .AND. ( .NOT. RSVEC) ) {
*.......................................................................
        // .. the singular values and the left singular vectors requested
*.......................................................................""""""""
         if ( RTRANS ) {
             // .. apply CGESVD to R**H
             // .. copy R**H into [U] and overwrite [U] with the right singular
             // vectors of R
            DO 1192 p = 1, NR
               DO 1193 q = p, N
                  U(q,p) = CONJG(A(p,q))
 1193          CONTINUE
 1192       CONTINUE
            IF ( NR .GT. 1 ) CALL CLASET( 'U', NR-1,NR-1, CZERO,CZERO, U(1,2), LDU )
            // .. the left singular vectors not computed, the NR right singular
            // vectors overwrite [U](1:NR,1:NR) as conjugate transposed. These
            // will be pre-multiplied by Q to build the left singular vectors of A.
               CALL CGESVD( 'N', 'O', N, NR, U, LDU, S, U, LDU, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO )

               DO 1119 p = 1, NR
                   U(p,p) = CONJG(U(p,p))
                   DO 1120 q = p + 1, NR
                      CTMP   = CONJG(U(q,p))
                      U(q,p) = CONJG(U(p,q))
                      U(p,q) = CTMP
 1120              CONTINUE
 1119          CONTINUE

         } else {
             // .. apply CGESVD to R
             // .. copy R into [U] and overwrite [U] with the left singular vectors
             CALL CLACPY( 'U', NR, N, A, LDA, U, LDU )
             IF ( NR .GT. 1 ) CALL CLASET( 'L', NR-1, NR-1, CZERO, CZERO, U(2,1), LDU )
             // .. the right singular vectors not computed, the NR left singular
             // vectors overwrite [U](1:NR,1:NR)
                CALL CGESVD( 'O', 'N', NR, N, U, LDU, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO )
                // .. now [U](1:NR,1:NR) contains the NR left singular vectors of
                // R. These will be pre-multiplied by Q to build the left singular
                // vectors of A.
         }

            // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
         if ( ( NR .LT. M ) .AND. ( .NOT.WNTUF ) ) {
             CALL CLASET('A', M-NR, NR, CZERO, CZERO, U(NR+1,1), LDU)
             if ( NR .LT. N1 ) {
                CALL CLASET( 'A',NR,N1-NR,CZERO,CZERO,U(1,NR+1), LDU )
                CALL CLASET( 'A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU )
             }
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         IF ( .NOT.WNTUF ) CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LCWORK-N, IERR )
         IF ( ROWPRM .AND. .NOT.WNTUF ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 )

      } else if ( RSVEC .AND. ( .NOT. LSVEC ) ) {
*.......................................................................
        // .. the singular values and the right singular vectors requested
*.......................................................................
          if ( RTRANS ) {
             // .. apply CGESVD to R**H
             // .. copy R**H into V and overwrite V with the left singular vectors
            DO 1165 p = 1, NR
               DO 1166 q = p, N
                  V(q,p) = CONJG(A(p,q))
 1166          CONTINUE
 1165       CONTINUE
            IF ( NR .GT. 1 ) CALL CLASET( 'U', NR-1,NR-1, CZERO,CZERO, V(1,2), LDV )
            // .. the left singular vectors of R**H overwrite V, the right singular
            // vectors not computed
            if ( WNTVR .OR. ( NR .EQ. N ) ) {
               CALL CGESVD( 'O', 'N', N, NR, V, LDV, S, U, LDU, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO )

               DO 1121 p = 1, NR
                   V(p,p) = CONJG(V(p,p))
                   DO 1122 q = p + 1, NR
                      CTMP   = CONJG(V(q,p))
                      V(q,p) = CONJG(V(p,q))
                      V(p,q) = CTMP
 1122              CONTINUE
 1121          CONTINUE

               if ( NR .LT. N ) {
                   DO 1103 p = 1, NR
                      DO 1104 q = NR + 1, N
                          V(p,q) = CONJG(V(q,p))
 1104                 CONTINUE
 1103              CONTINUE
               }
               CALL CLAPMT( .FALSE., NR, N, V, LDV, IWORK )
            } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:N,1:NR)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the QR factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                CALL CLASET('G', N, N-NR, CZERO, CZERO, V(1,NR+1), LDV)
                CALL CGESVD( 'O', 'N', N, N, V, LDV, S, U, LDU, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO )

                DO 1123 p = 1, N
                   V(p,p) = CONJG(V(p,p))
                   DO 1124 q = p + 1, N
                      CTMP   = CONJG(V(q,p))
                      V(q,p) = CONJG(V(p,q))
                      V(p,q) = CTMP
 1124              CONTINUE
 1123           CONTINUE
                CALL CLAPMT( .FALSE., N, N, V, LDV, IWORK )
            }

          } else {
             // .. aply CGESVD to R
             // .. copy R into V and overwrite V with the right singular vectors
             CALL CLACPY( 'U', NR, N, A, LDA, V, LDV )
             IF ( NR .GT. 1 ) CALL CLASET( 'L', NR-1, NR-1, CZERO, CZERO, V(2,1), LDV )
             // .. the right singular vectors overwrite V, the NR left singular
             // vectors stored in U(1:NR,1:NR)
             if ( WNTVR .OR. ( NR .EQ. N ) ) {
                CALL CGESVD( 'N', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO )
                CALL CLAPMT( .FALSE., NR, N, V, LDV, IWORK )
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
             } else {
                // .. need all N right singular vectors and NR < N
                // [!] This is simple implementation that augments [V](1:NR,1:N)
                // by padding a zero block. In the case NR << N, a more efficient
                // way is to first use the LQ factorization. For more details
                // how to implement this, see the " FULL SVD " branch.
                 CALL CLASET('G', N-NR, N, CZERO,CZERO, V(NR+1,1), LDV)
                 CALL CGESVD( 'N', 'O', N, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO )
                 CALL CLAPMT( .FALSE., N, N, V, LDV, IWORK )
             }
             // .. now [V] contains the adjoint of the matrix of the right singular
             // vectors of A.
          }

      } else {
*.......................................................................
        // .. FULL SVD requested
*.......................................................................
         if ( RTRANS ) {

             // .. apply CGESVD to R**H [[this option is left for R&D&T]]

            if ( WNTVR .OR. ( NR .EQ. N ) ) {
             // .. copy R**H into [V] and overwrite [V] with the left singular
             // vectors of R**H
            DO 1168 p = 1, NR
               DO 1169 q = p, N
                  V(q,p) = CONJG(A(p,q))
 1169          CONTINUE
 1168       CONTINUE
            IF ( NR .GT. 1 ) CALL CLASET( 'U', NR-1,NR-1, CZERO,CZERO, V(1,2), LDV )

            // .. the left singular vectors of R**H overwrite [V], the NR right
            // singular vectors of R**H stored in [U](1:NR,1:NR) as conjugate
            // transposed
               CALL CGESVD( 'O', 'A', N, NR, V, LDV, S, V, LDV, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO )
               // .. assemble V
               DO 1115 p = 1, NR
                  V(p,p) = CONJG(V(p,p))
                  DO 1116 q = p + 1, NR
                     CTMP   = CONJG(V(q,p))
                     V(q,p) = CONJG(V(p,q))
                     V(p,q) = CTMP
 1116             CONTINUE
 1115          CONTINUE
               if ( NR .LT. N ) {
                   DO 1101 p = 1, NR
                      DO 1102 q = NR+1, N
                         V(p,q) = CONJG(V(q,p))
 1102                 CONTINUE
 1101              CONTINUE
               }
               CALL CLAPMT( .FALSE., NR, N, V, LDV, IWORK )

                DO 1117 p = 1, NR
                   U(p,p) = CONJG(U(p,p))
                   DO 1118 q = p + 1, NR
                      CTMP   = CONJG(U(q,p))
                      U(q,p) = CONJG(U(p,q))
                      U(p,q) = CTMP
 1118              CONTINUE
 1117           CONTINUE

                if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                  CALL CLASET('A', M-NR,NR, CZERO,CZERO, U(NR+1,1), LDU)
                  if ( NR .LT. N1 ) {
                     CALL CLASET('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU)
                     CALL CLASET( 'A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU )
                  }
               }

            } else {
                // .. need all N right singular vectors and NR < N
             // .. copy R**H into [V] and overwrite [V] with the left singular
             // vectors of R**H
                // [[The optimal ratio N/NR for using QRF instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'CGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = MAX( OPTRATIO, 2 )
                OPTRATIO = 2
                if ( OPTRATIO*NR .GT. N ) {
                   DO 1198 p = 1, NR
                      DO 1199 q = p, N
                         V(q,p) = CONJG(A(p,q))
 1199                 CONTINUE
 1198              CONTINUE
                   IF ( NR .GT. 1 ) CALL CLASET('U',NR-1,NR-1, CZERO,CZERO, V(1,2),LDV)

                   CALL CLASET('A',N,N-NR,CZERO,CZERO,V(1,NR+1),LDV)
                   CALL CGESVD( 'O', 'A', N, N, V, LDV, S, V, LDV, U, LDU, CWORK(N+1), LCWORK-N, RWORK, INFO )

                   DO 1113 p = 1, N
                      V(p,p) = CONJG(V(p,p))
                      DO 1114 q = p + 1, N
                         CTMP   = CONJG(V(q,p))
                         V(q,p) = CONJG(V(p,q))
                         V(p,q) = CTMP
 1114                 CONTINUE
 1113              CONTINUE
                   CALL CLAPMT( .FALSE., N, N, V, LDV, IWORK )
               // .. assemble the left singular vector matrix U of dimensions
               // (M x N1), i.e. (M x N) or (M x M).

                   DO 1111 p = 1, N
                      U(p,p) = CONJG(U(p,p))
                      DO 1112 q = p + 1, N
                         CTMP   = CONJG(U(q,p))
                         U(q,p) = CONJG(U(p,q))
                         U(p,q) = CTMP
 1112                 CONTINUE
 1111              CONTINUE

                   if ( ( N .LT. M ) .AND. .NOT.(WNTUF)) {
                      CALL CLASET('A',M-N,N,CZERO,CZERO,U(N+1,1),LDU)
                      if ( N .LT. N1 ) {
                        CALL CLASET('A',N,N1-N,CZERO,CZERO,U(1,N+1),LDU)
                        CALL CLASET('A',M-N,N1-N,CZERO,CONE, U(N+1,N+1), LDU )
                      }
                   }
                } else {
                   // .. copy R**H into [U] and overwrite [U] with the right
                   // singular vectors of R
                   DO 1196 p = 1, NR
                      DO 1197 q = p, N
                         U(q,NR+p) = CONJG(A(p,q))
 1197                 CONTINUE
 1196              CONTINUE
                   IF ( NR .GT. 1 ) CALL CLASET('U',NR-1,NR-1,CZERO,CZERO,U(1,NR+2),LDU)                    CALL CGEQRF( N, NR, U(1,NR+1), LDU, CWORK(N+1), CWORK(N+NR+1), LCWORK-N-NR, IERR )
                   DO 1143 p = 1, NR
                       DO 1144 q = 1, N
                           V(q,p) = CONJG(U(p,NR+q))
 1144                  CONTINUE
 1143              CONTINUE
                  CALL CLASET('U',NR-1,NR-1,CZERO,CZERO,V(1,2),LDV)
                  CALL CGESVD( 'S', 'O', NR, NR, V, LDV, S, U, LDU, V,LDV, CWORK(N+NR+1),LCWORK-N-NR,RWORK, INFO )
                  CALL CLASET('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV)
                  CALL CLASET('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV)
                  CALL CLASET('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV)
                  CALL CUNMQR('R','C', N, N, NR, U(1,NR+1), LDU, CWORK(N+1),V,LDV,CWORK(N+NR+1),LCWORK-N-NR,IERR)
                  CALL CLAPMT( .FALSE., N, N, V, LDV, IWORK )
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x NR) or (M x N) or (M x M).
                  if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                     CALL CLASET('A',M-NR,NR,CZERO,CZERO,U(NR+1,1),LDU)
                     if ( NR .LT. N1 ) {
                     CALL CLASET('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU)
                     CALL CLASET( 'A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1),LDU)
                     }
                  }
                }
            }

         } else {

             // .. apply CGESVD to R [[this is the recommended option]]

             if ( WNTVR .OR. ( NR .EQ. N ) ) {
                 // .. copy R into [V] and overwrite V with the right singular vectors
                 CALL CLACPY( 'U', NR, N, A, LDA, V, LDV )
                IF ( NR .GT. 1 ) CALL CLASET( 'L', NR-1,NR-1, CZERO,CZERO, V(2,1), LDV )
                // .. the right singular vectors of R overwrite [V], the NR left
                // singular vectors of R stored in [U](1:NR,1:NR)
                CALL CGESVD( 'S', 'O', NR, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO )
                CALL CLAPMT( .FALSE., NR, N, V, LDV, IWORK )
                // .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
               if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                  CALL CLASET('A', M-NR,NR, CZERO,CZERO, U(NR+1,1), LDU)
                  if ( NR .LT. N1 ) {
                     CALL CLASET('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU)
                     CALL CLASET( 'A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU )
                  }
               }

             } else {
               // .. need all N right singular vectors and NR < N
               // .. the requested number of the left singular vectors
                // is then N1 (N or M)
                // [[The optimal ratio N/NR for using LQ instead of padding
                  // with zeros. Here hard coded to 2; it must be at least
                  // two due to work space constraints.]]
                // OPTRATIO = ILAENV(6, 'CGESVD', 'S' // 'O', NR,N,0,0)
                // OPTRATIO = MAX( OPTRATIO, 2 )
               OPTRATIO = 2
               if ( OPTRATIO * NR .GT. N ) {
                  CALL CLACPY( 'U', NR, N, A, LDA, V, LDV )
                  IF ( NR .GT. 1 ) CALL CLASET('L', NR-1,NR-1, CZERO,CZERO, V(2,1),LDV)
               // .. the right singular vectors of R overwrite [V], the NR left
                  // singular vectors of R stored in [U](1:NR,1:NR)
                  CALL CLASET('A', N-NR,N, CZERO,CZERO, V(NR+1,1),LDV)
                  CALL CGESVD( 'S', 'O', N, N, V, LDV, S, U, LDU, V, LDV, CWORK(N+1), LCWORK-N, RWORK, INFO )
                  CALL CLAPMT( .FALSE., N, N, V, LDV, IWORK )
                  // .. now [V] contains the adjoint of the matrix of the right
                  // singular vectors of A. The leading N left singular vectors
                  // are in [U](1:N,1:N)
                  // .. assemble the left singular vector matrix U of dimensions
                  // (M x N1), i.e. (M x N) or (M x M).
                  if ( ( N .LT. M ) .AND. .NOT.(WNTUF)) {
                      CALL CLASET('A',M-N,N,CZERO,CZERO,U(N+1,1),LDU)
                      if ( N .LT. N1 ) {
                        CALL CLASET('A',N,N1-N,CZERO,CZERO,U(1,N+1),LDU)
                        CALL CLASET( 'A',M-N,N1-N,CZERO,CONE, U(N+1,N+1), LDU )
                      }
                  }
               } else {
                  CALL CLACPY( 'U', NR, N, A, LDA, U(NR+1,1), LDU )
                  IF ( NR .GT. 1 ) CALL CLASET('L',NR-1,NR-1,CZERO,CZERO,U(NR+2,1),LDU)                   CALL CGELQF( NR, N, U(NR+1,1), LDU, CWORK(N+1), CWORK(N+NR+1), LCWORK-N-NR, IERR )
                  CALL CLACPY('L',NR,NR,U(NR+1,1),LDU,V,LDV)
                  IF ( NR .GT. 1 ) CALL CLASET('U',NR-1,NR-1,CZERO,CZERO,V(1,2),LDV)                   CALL CGESVD( 'S', 'O', NR, NR, V, LDV, S, U, LDU, V, LDV, CWORK(N+NR+1), LCWORK-N-NR, RWORK, INFO )
                  CALL CLASET('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV)
                  CALL CLASET('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV)
                  CALL CLASET('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV)
                  CALL CUNMLQ('R','N',N,N,NR,U(NR+1,1),LDU,CWORK(N+1), V, LDV, CWORK(N+NR+1),LCWORK-N-NR,IERR)
                  CALL CLAPMT( .FALSE., N, N, V, LDV, IWORK )
                // .. assemble the left singular vector matrix U of dimensions
               // (M x NR) or (M x N) or (M x M).
                  if ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) {
                     CALL CLASET('A',M-NR,NR,CZERO,CZERO,U(NR+1,1),LDU)
                     if ( NR .LT. N1 ) {
                     CALL CLASET('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU)
                     CALL CLASET( 'A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1), LDU )
                     }
                  }
               }
             }
         // .. end of the "R**H or R" branch
         }

            // The Q matrix from the first QRF is built into the left singular
            // vectors matrix U.

         IF ( .NOT. WNTUF ) CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LCWORK-N, IERR )
         IF ( ROWPRM .AND. .NOT.WNTUF ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 )

      // ... end of the "full SVD" branch
      }

      // Check whether some singular values are returned as zeros, e.g.
      // due to underflow, and update the numerical rank.
      p = NR
      DO 4001 q = p, 1, -1
          IF ( S(q) .GT. ZERO ) GO TO 4002
          NR = NR - 1
 4001 CONTINUE
 4002 CONTINUE

      // .. if numerical rank deficiency is detected, the truncated
      // singular values are set to zero.
      IF ( NR .LT. N ) CALL SLASET( 'G', N-NR,1, ZERO,ZERO, S(NR+1), N )
      // .. undo scaling; this may cause overflow in the largest singular
      // values.
      IF ( ASCALED ) CALL SLASCL( 'G',0,0, ONE,SQRT(REAL(M)), NR,1, S, N, IERR )
      IF ( CONDA ) RWORK(1) = SCONDA
      RWORK(2) = p - NR
      // .. p-NR is the number of singular values that are computed as
      // exact zeros in CGESVD() applied to the (possibly truncated)
      // full row rank triangular (trapezoidal) factor of A.
      NUMRANK = NR

      RETURN

      // End of CGESVDQ

      }
