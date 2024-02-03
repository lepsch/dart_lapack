      SUBROUTINE CGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, M, N, A, LDA, SVA, U, LDU, V, LDV, CWORK, LWORK, RWORK, LRWORK, IWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      IMPLICIT    NONE;
      int         INFO, LDA, LDU, LDV, LWORK, LRWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX     A( LDA, * ), U( LDU, * ), V( LDV, * ), CWORK( LWORK );
      REAL        SVA( N ), RWORK( LRWORK );
      int         IWORK( * );
      String      JOBA, JOBP, JOBR, JOBT, JOBU, JOBV;
      // ..

*  ===========================================================================

      // .. Local Parameters ..
      REAL        ZERO,         ONE;
      const     ZERO = 0.0, ONE = 1.0 ;
      COMPLEX     CZERO,                    CONE;
      const     CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX CTEMP;
      REAL    AAPP,   AAQQ,   AATMAX, AATMIN, BIG,    BIG1,   COND_OK, CONDR1, CONDR2, ENTRA,  ENTRAT, EPSLN,  MAXPRJ, SCALEM, SCONDA, SFMIN,  SMALL,  TEMP1,  USCAL1, USCAL2, XSC;
      int     IERR,   N1,     NR,     NUMRANK,        p, q,   WARNING;
      bool    ALMORT, DEFR,   ERREST, GOSCAL,  JRACC,  KILL,   LQUERY, LSVEC,  L2ABER, L2KILL, L2PERT,  L2RANK, L2TRAN, NOSCAL, ROWPIV, RSVEC,  TRANSP;

      int     OPTWRK, MINWRK, MINRWRK, MINIWRK;
      int     LWCON,  LWLQF, LWQP3, LWQRF, LWUNMLQ, LWUNMQR, LWUNMQRM, LWSVDJ, LWSVDJV, LRWQP3, LRWCON, LRWSVDJ, IWOFF;
      int     LWRK_CGELQF, LWRK_CGEQP3,  LWRK_CGEQP3N, LWRK_CGEQRF,   LWRK_CGESVJ, LWRK_CGESVJV, LWRK_CGESVJU, LWRK_CUNMLQ, LWRK_CUNMQR, LWRK_CUNMQRM;
      // ..
      // .. Local Arrays
      COMPLEX CDUMMY(1);
      REAL    RDUMMY(1);

      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG, ALOG, MAX, MIN, REAL, NINT, SQRT
      // ..
      // .. External Functions ..
      REAL      SLAMCH, SCNRM2;
      int       ISAMAX, ICAMAX;
      bool      LSAME;
      // EXTERNAL ISAMAX, ICAMAX, LSAME, SLAMCH, SCNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ, CCOPY,  CGELQF, CGEQP3, CGEQRF, CLACPY, CLAPMR, CLASCL, SLASCL, CLASET, CLASSQ, CLASWP, CUNGQR, CUNMLQ, CUNMQR, CPOCON, SSCAL,  CSSCAL, CSWAP,  CTRSM,  CLACGV, XERBLA

      // EXTERNAL CGESVJ
      // ..

      // Test the input arguments

      LSVEC  = LSAME( JOBU, 'U' ) || LSAME( JOBU, 'F' );
      JRACC  = LSAME( JOBV, 'J' );
      RSVEC  = LSAME( JOBV, 'V' ) || JRACC;
      ROWPIV = LSAME( JOBA, 'F' ) || LSAME( JOBA, 'G' );
      L2RANK = LSAME( JOBA, 'R' );
      L2ABER = LSAME( JOBA, 'A' );
      ERREST = LSAME( JOBA, 'E' ) || LSAME( JOBA, 'G' );
      L2TRAN = LSAME( JOBT, 'T' ) && ( M == N );
      L2KILL = LSAME( JOBR, 'R' );
      DEFR   = LSAME( JOBR, 'N' );
      L2PERT = LSAME( JOBP, 'P' );

      LQUERY = ( LWORK == -1 ) || ( LRWORK == -1 );

      if ( !(ROWPIV || L2RANK || L2ABER || ERREST || LSAME( JOBA, 'C' ) )) {
         INFO = - 1;
      } else if ( !( LSVEC || LSAME( JOBU, 'N' ) || ( LSAME( JOBU, 'W' ) && RSVEC && L2TRAN ) ) ) {
         INFO = - 2;
      } else if ( !( RSVEC || LSAME( JOBV, 'N' ) || ( LSAME( JOBV, 'W' ) && LSVEC && L2TRAN ) ) ) {
         INFO = - 3;
      } else if ( !( L2KILL || DEFR ) ) {
         INFO = - 4;
      } else if ( !( LSAME(JOBT,'T') || LSAME(JOBT,'N') ) ) {
         INFO = - 5;
      } else if ( !( L2PERT || LSAME( JOBP, 'N' ) ) ) {
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
      } else {
         // #:)
         INFO = 0;
      }

      if ( INFO == 0 ) {
          // .. compute the minimal and the optimal workspace lengths
          // [[The expressions for computing the minimal and the optimal
          // values of LCWORK, LRWORK are written with a lot of redundancy and
          // can be simplified. However, this verbose form is useful for
          // maintenance and modifications of the code.]]

         // .. minimal workspace length for CGEQP3 of an M x N matrix,
          // CGEQRF of an N x N matrix, CGELQF of an N x N matrix,
          // CUNMLQ for computing N x N matrix, CUNMQR for computing N x N
          // matrix, CUNMQR for computing M x N matrix, respectively.
          LWQP3 = N+1;
          LWQRF = MAX( 1, N );
          LWLQF = MAX( 1, N );
          LWUNMLQ  = MAX( 1, N );
          LWUNMQR  = MAX( 1, N );
          LWUNMQRM = MAX( 1, M );
         // .. minimal workspace length for CPOCON of an N x N matrix
          LWCON = 2 * N;
         // .. minimal workspace length for CGESVJ of an N x N matrix,
          // without and with explicit accumulation of Jacobi rotations
          LWSVDJ  = MAX( 2 * N, 1 );
          LWSVDJV = MAX( 2 * N, 1 );
          // .. minimal REAL workspace length for CGEQP3, CPOCON, CGESVJ
          LRWQP3  = 2 * N;
          LRWCON  = N;
          LRWSVDJ = N;
          if ( LQUERY ) {
              cgeqp3(M, N, A, LDA, IWORK, CDUMMY, CDUMMY, -1,  RDUMMY, IERR );
              LWRK_CGEQP3 = INT( CDUMMY(1) );
              cgeqrf(N, N, A, LDA, CDUMMY, CDUMMY,-1, IERR );
              LWRK_CGEQRF = INT( CDUMMY(1) );
              cgelqf(N, N, A, LDA, CDUMMY, CDUMMY,-1, IERR );
              LWRK_CGELQF = INT( CDUMMY(1) );
          }
          MINWRK  = 2;
          OPTWRK  = 2;
          MINIWRK = N;
          if ( !(LSVEC || RSVEC ) ) {
              // .. minimal and optimal sizes of the complex workspace if
              // only the singular values are requested
              if ( ERREST ) {
                  MINWRK = MAX( N+LWQP3, N**2+LWCON, N+LWQRF, LWSVDJ );
              } else {
                  MINWRK = MAX( N+LWQP3, N+LWQRF, LWSVDJ );
              }
              if ( LQUERY ) {
                  cgesvj('L', 'N', 'N', N, N, A, LDA, SVA, N, V,  LDV, CDUMMY, -1, RDUMMY, -1, IERR );
                  LWRK_CGESVJ = INT( CDUMMY(1) );
                  if ( ERREST ) {
                      OPTWRK = MAX( N+LWRK_CGEQP3, N**2+LWCON,  N+LWRK_CGEQRF, LWRK_CGESVJ );
                  } else {
                      OPTWRK = MAX( N+LWRK_CGEQP3, N+LWRK_CGEQRF,  LWRK_CGESVJ );
                  }
              }
              if ( L2TRAN || ROWPIV ) {
                  if ( ERREST ) {
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWCON, LRWSVDJ );
                  } else {
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ );
                  }
              } else {
                  if ( ERREST ) {
                     MINRWRK = MAX( 7, LRWQP3, LRWCON, LRWSVDJ );
                  } else {
                     MINRWRK = MAX( 7, LRWQP3, LRWSVDJ );
                  }
              }
              if (ROWPIV || L2TRAN) MINIWRK = MINIWRK + M;
          } else if ( RSVEC && ( !LSVEC) ) {
             // .. minimal and optimal sizes of the complex workspace if the
             // singular values and the right singular vectors are requested
             if ( ERREST ) {
                 MINWRK = MAX( N+LWQP3, LWCON, LWSVDJ, N+LWLQF,   2*N+LWQRF, N+LWSVDJ, N+LWUNMLQ );
             } else {
                 MINWRK = MAX( N+LWQP3, LWSVDJ, N+LWLQF, 2*N+LWQRF,  N+LWSVDJ, N+LWUNMLQ );
             }
             if ( LQUERY ) {
                 cgesvj('L', 'U', 'N', N,N, U, LDU, SVA, N, A, LDA, CDUMMY, -1, RDUMMY, -1, IERR );
                 LWRK_CGESVJ = INT( CDUMMY(1) );
                 cunmlq('L', 'C', N, N, N, A, LDA, CDUMMY, V, LDV, CDUMMY, -1, IERR );
                 LWRK_CUNMLQ = INT( CDUMMY(1) );
                 if ( ERREST ) {
                 OPTWRK = MAX( N+LWRK_CGEQP3, LWCON, LWRK_CGESVJ,  N+LWRK_CGELQF, 2*N+LWRK_CGEQRF, N+LWRK_CGESVJ,  N+LWRK_CUNMLQ );
                 } else {
                 OPTWRK = MAX( N+LWRK_CGEQP3, LWRK_CGESVJ,N+LWRK_CGELQF, 2*N+LWRK_CGEQRF, N+LWRK_CGESVJ, N+LWRK_CUNMLQ );
                 }
             }
             if ( L2TRAN || ROWPIV ) {
                  if ( ERREST ) {
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON );
                  } else {
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ );
                  }
             } else {
                  if ( ERREST ) {
                     MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON );
                  } else {
                     MINRWRK = MAX( 7, LRWQP3, LRWSVDJ );
                  }
             }
             if (ROWPIV || L2TRAN) MINIWRK = MINIWRK + M;
          } else if ( LSVEC && ( !RSVEC) ) {
             // .. minimal and optimal sizes of the complex workspace if the
             // singular values and the left singular vectors are requested
             if ( ERREST ) {
                 MINWRK = N + MAX( LWQP3,LWCON,N+LWQRF,LWSVDJ,LWUNMQRM );
             } else {
                 MINWRK = N + MAX( LWQP3, N+LWQRF, LWSVDJ, LWUNMQRM );
             }
             if ( LQUERY ) {
                 cgesvj('L', 'U', 'N', N,N, U, LDU, SVA, N, A, LDA, CDUMMY, -1, RDUMMY, -1, IERR );
                 LWRK_CGESVJ = INT( CDUMMY(1) );
                 cunmqr('L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR );
                 LWRK_CUNMQRM = INT( CDUMMY(1) );
                 if ( ERREST ) {
                 OPTWRK = N + MAX( LWRK_CGEQP3, LWCON, N+LWRK_CGEQRF, LWRK_CGESVJ, LWRK_CUNMQRM );
                 } else {
                 OPTWRK = N + MAX( LWRK_CGEQP3, N+LWRK_CGEQRF, LWRK_CGESVJ, LWRK_CUNMQRM );
                 }
             }
             if ( L2TRAN || ROWPIV ) {
                 if ( ERREST ) {
                    MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON );
                 } else {
                    MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ );
                 }
             } else {
                 if ( ERREST ) {
                    MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON );
                 } else {
                    MINRWRK = MAX( 7, LRWQP3, LRWSVDJ );
                 }
             }
             if (ROWPIV || L2TRAN) MINIWRK = MINIWRK + M;
          } else {
             // .. minimal and optimal sizes of the complex workspace if the
             // full SVD is requested
             if ( !JRACC ) {
                 if ( ERREST ) {
                    MINWRK = MAX( N+LWQP3, N+LWCON,  2*N+N**2+LWCON,  2*N+LWQRF,         2*N+LWQP3, 2*N+N**2+N+LWLQF,  2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWSVDJ, 2*N+N**2+N+LWSVDJV, 2*N+N**2+N+LWUNMQR,2*N+N**2+N+LWUNMLQ, N+N**2+LWSVDJ,   N+LWUNMQRM );
                 } else {
                    MINWRK = MAX( N+LWQP3,        2*N+N**2+LWCON,  2*N+LWQRF,         2*N+LWQP3, 2*N+N**2+N+LWLQF,  2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWSVDJ, 2*N+N**2+N+LWSVDJV, 2*N+N**2+N+LWUNMQR,2*N+N**2+N+LWUNMLQ, N+N**2+LWSVDJ,      N+LWUNMQRM );
                 }
                 MINIWRK = MINIWRK + N;
                 if (ROWPIV || L2TRAN) MINIWRK = MINIWRK + M;
             } else {
                 if ( ERREST ) {
                    MINWRK = MAX( N+LWQP3, N+LWCON, 2*N+LWQRF,  2*N+N**2+LWSVDJV, 2*N+N**2+N+LWUNMQR, N+LWUNMQRM );
                 } else {
                    MINWRK = MAX( N+LWQP3, 2*N+LWQRF,  2*N+N**2+LWSVDJV, 2*N+N**2+N+LWUNMQR, N+LWUNMQRM );
                 }
                 if (ROWPIV || L2TRAN) MINIWRK = MINIWRK + M;
             }
             if ( LQUERY ) {
                 cunmqr('L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR );
                 LWRK_CUNMQRM = INT( CDUMMY(1) );
                 cunmqr('L', 'N', N, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR );
                 LWRK_CUNMQR = INT( CDUMMY(1) );
                 if ( !JRACC ) {
                     cgeqp3(N,N, A, LDA, IWORK, CDUMMY,CDUMMY, -1, RDUMMY, IERR );
                     LWRK_CGEQP3N = INT( CDUMMY(1) );
                     cgesvj('L', 'U', 'N', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR );
                     LWRK_CGESVJ = INT( CDUMMY(1) );
                     cgesvj('U', 'U', 'N', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR );
                     LWRK_CGESVJU = INT( CDUMMY(1) );
                     cgesvj('L', 'U', 'V', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR );
                     LWRK_CGESVJV = INT( CDUMMY(1) );
                     cunmlq('L', 'C', N, N, N, A, LDA, CDUMMY, V, LDV, CDUMMY, -1, IERR );
                     LWRK_CUNMLQ = INT( CDUMMY(1) );
                     if ( ERREST ) {
                       OPTWRK = MAX( N+LWRK_CGEQP3, N+LWCON,  2*N+N**2+LWCON, 2*N+LWRK_CGEQRF, 2*N+LWRK_CGEQP3N, 2*N+N**2+N+LWRK_CGELQF, 2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWRK_CGESVJ, 2*N+N**2+N+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR, 2*N+N**2+N+LWRK_CUNMLQ, N+N**2+LWRK_CGESVJU, N+LWRK_CUNMQRM );
                     } else {
                       OPTWRK = MAX( N+LWRK_CGEQP3,   2*N+N**2+LWCON, 2*N+LWRK_CGEQRF, 2*N+LWRK_CGEQP3N, 2*N+N**2+N+LWRK_CGELQF, 2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWRK_CGESVJ, 2*N+N**2+N+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR, 2*N+N**2+N+LWRK_CUNMLQ, N+N**2+LWRK_CGESVJU, N+LWRK_CUNMQRM );
                     }
                 } else {
                     cgesvj('L', 'U', 'V', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR );
                     LWRK_CGESVJV = INT( CDUMMY(1) );
                     cunmqr('L', 'N', N, N, N, CDUMMY, N, CDUMMY, V, LDV, CDUMMY, -1, IERR );
                     LWRK_CUNMQR = INT( CDUMMY(1) );
                     cunmqr('L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR );
                     LWRK_CUNMQRM = INT( CDUMMY(1) );
                     if ( ERREST ) {
                        OPTWRK = MAX( N+LWRK_CGEQP3, N+LWCON,    2*N+LWRK_CGEQRF, 2*N+N**2, 2*N+N**2+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR,N+LWRK_CUNMQRM );
                     } else {
                        OPTWRK = MAX( N+LWRK_CGEQP3, 2*N+LWRK_CGEQRF,   2*N+N**2, 2*N+N**2+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR, N+LWRK_CUNMQRM );
                     }
                 }
             }
             if ( L2TRAN || ROWPIV ) {
                 MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON );
             } else {
                 MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON );
             }
          }
          MINWRK = MAX( 2, MINWRK );
          OPTWRK = MAX( OPTWRK, MINWRK );
          IF ( LWORK < MINWRK && ( !LQUERY) ) INFO = - 17;
          IF ( LRWORK < MINRWRK && ( !LQUERY) ) INFO = - 19;
      }

      if ( INFO != 0 ) {
        // #:(
         xerbla('CGEJSV', - INFO );
         RETURN;
      } else if ( LQUERY ) {
          CWORK(1) = OPTWRK;
          CWORK(2) = MINWRK;
          RWORK(1) = MINRWRK;
          IWORK(1) = MAX( 4, MINIWRK );
          RETURN;
      }

      // Quick return for void matrix (Y3K safe)
* #:)
      if ( ( M == 0 ) || ( N == 0 ) ) {
         IWORK(1:4) = 0;
         RWORK(1:7) = 0;
         RETURN;
      }

      // Determine whether the matrix U should be M x N or M x M

      if ( LSVEC ) {
         N1 = N;
         IF ( LSAME( JOBU, 'F' ) ) N1 = M;
      }

      // Set numerical parameters

*!    NOTE: Make sure SLAMCH() does not fail on the target architecture.

      EPSLN = SLAMCH('Epsilon');
      SFMIN = SLAMCH('SafeMinimum');
      SMALL = SFMIN / EPSLN;
      BIG   = SLAMCH('O');
      // BIG   = ONE / SFMIN

      // Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N

*(!)  If necessary, scale SVA() to protect the largest norm from
      // overflow. It is possible that this scaling pushes the smallest
      // column norm left from the underflow threshold (extreme case).

      SCALEM  = ONE / SQRT(REAL(M)*REAL(N));
      NOSCAL  = true;
      GOSCAL  = true;
      for (p = 1; p <= N; p++) { // 1874
         AAPP = ZERO;
         AAQQ = ONE;
         classq(M, A(1,p), 1, AAPP, AAQQ );
         if ( AAPP > BIG ) {
            INFO = - 9;
            xerbla('CGEJSV', -INFO );
            RETURN;
         }
         AAQQ = SQRT(AAQQ);
         if ( ( AAPP < (BIG / AAQQ) ) && NOSCAL  ) {
            SVA(p)  = AAPP * AAQQ;
         } else {
            NOSCAL  = false;
            SVA(p)  = AAPP * ( AAQQ * SCALEM );
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
         AAPP = MAX( AAPP, SVA(p) );
         IF ( SVA(p) != ZERO ) AAQQ = MIN( AAQQ, SVA(p) );
      } // 4781

      // Quick return for zero M x N matrix
* #:)
      if ( AAPP == ZERO ) {
         if (LSVEC) CALL CLASET( 'G', M, N1, CZERO, CONE, U, LDU );
         if (RSVEC) CALL CLASET( 'G', N, N,  CZERO, CONE, V, LDV );
         RWORK(1) = ONE;
         RWORK(2) = ONE;
         if (ERREST) RWORK(3) = ONE;
         if ( LSVEC && RSVEC ) {
            RWORK(4) = ONE;
            RWORK(5) = ONE;
         }
         if ( L2TRAN ) {
            RWORK(6) = ZERO;
            RWORK(7) = ZERO;
         }
         IWORK(1) = 0;
         IWORK(2) = 0;
         IWORK(3) = 0;
         IWORK(4) = -1;
         RETURN;
      }

      // Issue warning if denormalized column norms detected. Override the
      // high relative accuracy request. Issue licence to kill nonzero columns
      // (set them to zero) whose norm is less than sigma_max / BIG (roughly).
* #:(
      WARNING = 0;
      if ( AAQQ <= SFMIN ) {
         L2RANK = true;
         L2KILL = true;
         WARNING = 1;
      }

      // Quick return for one-column matrix
* #:)
      if ( N == 1 ) {

         if ( LSVEC ) {
            clascl('G',0,0,SVA(1),SCALEM, M,1,A(1,1),LDA,IERR );
            clacpy('A', M, 1, A, LDA, U, LDU );
            // computing all M left singular vectors of the M x 1 matrix
            if ( N1 != N  ) {
              cgeqrf(M, N, U,LDU, CWORK, CWORK(N+1),LWORK-N,IERR );
              cungqr(M,N1,1, U,LDU,CWORK,CWORK(N+1),LWORK-N,IERR );
              ccopy(M, A(1,1), 1, U(1,1), 1 );
            }
         }
         if ( RSVEC ) {
             V(1,1) = CONE;
         }
         if ( SVA(1) < (BIG*SCALEM) ) {
            SVA(1)  = SVA(1) / SCALEM;
            SCALEM  = ONE;
         }
         RWORK(1) = ONE / SCALEM;
         RWORK(2) = ONE;
         if ( SVA(1) != ZERO ) {
            IWORK(1) = 1;
            if ( ( SVA(1) / SCALEM) >= SFMIN ) {
               IWORK(2) = 1;
            } else {
               IWORK(2) = 0;
            }
         } else {
            IWORK(1) = 0;
            IWORK(2) = 0;
         }
         IWORK(3) = 0;
         IWORK(4) = -1;
         if (ERREST) RWORK(3) = ONE;
         if ( LSVEC && RSVEC ) {
            RWORK(4) = ONE;
            RWORK(5) = ONE;
         }
         if ( L2TRAN ) {
            RWORK(6) = ZERO;
            RWORK(7) = ZERO;
         }
         RETURN;

      }

      TRANSP = false;

      AATMAX = -ONE;
      AATMIN =  BIG;
      if ( ROWPIV || L2TRAN ) {

      // Compute the row norms, needed to determine row pivoting sequence
      // (in the case of heavily row weighted A, row pivoting is strongly
      // advised) and to collect information needed to compare the
      // structures of A * A^* and A^* * A (in the case L2TRAN == true ).

         if ( L2TRAN ) {
            for (p = 1; p <= M; p++) { // 1950
               XSC   = ZERO;
               TEMP1 = ONE;
               classq(N, A(p,1), LDA, XSC, TEMP1 );
               // CLASSQ gets both the ell_2 and the ell_infinity norm
               // in one pass through the vector
               RWORK(M+p)  = XSC * SCALEM;
               RWORK(p)    = XSC * (SCALEM*SQRT(TEMP1));
               AATMAX = MAX( AATMAX, RWORK(p) );
               IF (RWORK(p) != ZERO)  AATMIN = MIN(AATMIN,RWORK(p));
            } // 1950
         } else {
            for (p = 1; p <= M; p++) { // 1904
               RWORK(M+p) = SCALEM*ABS( A(p,ICAMAX(N,A(p,1),LDA)) );
               AATMAX = MAX( AATMAX, RWORK(M+p) );
               AATMIN = MIN( AATMIN, RWORK(M+p) );
            } // 1904
         }

      }

      // For square matrix A try to determine whether A^*  would be better
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
         ENTRA = - ENTRA / ALOG(REAL(N));

         // Now, SVA().^2/Trace(A^* * A) is a point in the probability simplex.
         // It is derived from the diagonal of  A^* * A.  Do the same with the
         // diagonal of A * A^*, compute the entropy of the corresponding
         // probability distribution. Note that A * A^* and A^* * A have the
         // same trace.

         ENTRAT = ZERO;
         for (p = 1; p <= M; p++) { // 1114
            BIG1 = ( ( RWORK(p) / XSC )**2 ) * TEMP1;
            if (BIG1 != ZERO) ENTRAT = ENTRAT + BIG1 * ALOG(BIG1);
         } // 1114
         ENTRAT = - ENTRAT / ALOG(REAL(M));

         // Analyze the entropies and decide A or A^*. Smaller entropy
         // usually means better input for the algorithm.

         TRANSP = ( ENTRAT < ENTRA );

         // If A^* is better than A, take the adjoint of A. This is allowed
         // only for square matrices, M=N.
         if ( TRANSP ) {
            // In an optimal implementation, this trivial transpose
            // should be replaced with faster transpose.
            for (p = 1; p <= N - 1; p++) { // 1115
               A(p,p) = CONJG(A(p,p));
               for (q = p + 1; q <= N; q++) { // 1116
                   CTEMP = CONJG(A(q,p));
                  A(q,p) = CONJG(A(p,q));
                  A(p,q) = CTEMP;
               } // 1116
            } // 1115
            A(N,N) = CONJG(A(N,N));
            for (p = 1; p <= N; p++) { // 1117
               RWORK(M+p) = SVA(p);
               SVA(p) = RWORK(p);
               // previously computed row 2-norms are now column 2-norms
               // of the transposed matrix
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
      // than SQRT(BIG) -- the matrix is scaled so that its maximal column
      // has Euclidean norm equal to SQRT(BIG/N). The only reason to keep
      // SQRT(BIG) instead of BIG is the fact that CGEJSV uses LAPACK and
      // BLAS routines that, in some implementations, are not capable of
      // working in the full interval [SFMIN,BIG] and that they may provoke
      // overflows in the intermediate results. If the singular values spread
      // from SFMIN to BIG, then CGESVJ will compute them. So, in that case,
      // one should use CGESVJ instead of CGEJSV.
      BIG1   = SQRT( BIG );
      TEMP1  = SQRT( BIG / REAL(N) );
      // >> for future updates: allow bigger range, i.e. the largest column
      // will be allowed up to BIG/N and CGESVJ will do the rest. However, for
      // this all other (LAPACK) components must allow such a range.
      // TEMP1  = BIG/REAL(N)
      // TEMP1  = BIG * EPSLN  this should 'almost' work with current LAPACK components
      slascl('G', 0, 0, AAPP, TEMP1, N, 1, SVA, N, IERR );
      if ( AAQQ > (AAPP * SFMIN) ) {
          AAQQ = ( AAQQ / AAPP ) * TEMP1;
      } else {
          AAQQ = ( AAQQ * TEMP1 ) / AAPP;
      }
      TEMP1 = TEMP1 * SCALEM;
      clascl('G', 0, 0, AAPP, TEMP1, M, N, A, LDA, IERR );

      // To undo scaling at the end of this procedure, multiply the
      // computed singular values with USCAL2 / USCAL1.

      USCAL1 = TEMP1;
      USCAL2 = AAPP;

      if ( L2KILL ) {
         // L2KILL enforces computation of nonzero singular values in
         // the restricted range of condition number of the initial A,
         // sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN).
         XSC = SQRT( SFMIN );
      } else {
         XSC = SMALL;

         // Now, if the condition number of A is too big,
         // sigma_max(A) / sigma_min(A) > SQRT(BIG/N) * EPSLN / SFMIN,
         // as a precaution measure, the full SVD is computed using CGESVJ
         // with accumulated Jacobi rotations. This provides numerically
         // more robust computation, at the cost of slightly increased run
         // time. Depending on the concrete implementation of BLAS and LAPACK
         // (i.e. how they behave in presence of extreme ill-conditioning) the
         // implementor may decide to remove this switch.
         if ( ( AAQQ < SQRT(SFMIN) ) && LSVEC && RSVEC ) {
            JRACC = true;
         }

      }
      if ( AAQQ < XSC ) {
         for (p = 1; p <= N; p++) { // 700
            if ( SVA(p) < XSC ) {
               claset('A', M, 1, CZERO, CZERO, A(1,p), LDA );
               SVA(p) = ZERO;
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
         if ( ( LSVEC && RSVEC ) && !( JRACC ) ) {
              IWOFF = 2*N;
         } else {
              IWOFF = N;
         }
         for (p = 1; p <= M - 1; p++) { // 1952
            q = ISAMAX( M-p+1, RWORK(M+p), 1 ) + p - 1;
            IWORK(IWOFF+p) = q;
            if ( p != q ) {
               TEMP1      = RWORK(M+p);
               RWORK(M+p) = RWORK(M+q);
               RWORK(M+q) = TEMP1;
            }
         } // 1952
         claswp(N, A, LDA, 1, M-1, IWORK(IWOFF+1), 1 );
      }

      // End of the preparation phase (scaling, optional sorting and
      // transposing, optional flushing of small columns).

      // Preconditioning

      // If the full SVD is needed, the right singular vectors are computed
      // from a matrix equation, and for that we need theoretical analysis
      // of the Businger-Golub pivoting. So we use CGEQP3 as the first RR QRF.
      // In all other cases the first RR QRF can be chosen by other criteria
      // (eg speed by replacing global with restricted window pivoting, such
      // as in xGEQPX from TOMS # 782). Good results will be obtained using
      // xGEQPX with properly (!) chosen numerical parameters.
      // Any improvement of CGEQP3 improves overall performance of CGEJSV.

      // A * P1 = Q1 * [ R1^* 0]^*:
      for (p = 1; p <= N; p++) { // 1963
         // .. all columns are free columns
         IWORK(p) = 0;
      } // 1963
      cgeqp3(M, N, A, LDA, IWORK, CWORK, CWORK(N+1), LWORK-N, RWORK, IERR );

      // The upper triangular matrix R1 from the first QRF is inspected for
      // rank deficiency and possibilities for deflation, or possible
      // ill-conditioning. Depending on the user specified flag L2RANK,
      // the procedure explores possibilities to reduce the numerical
      // rank by inspecting the computed upper triangular factor. If
      // L2RANK or L2ABER are up, then CGEJSV will compute the SVD of
      // A + dA, where ||dA|| <= f(M,N)*EPSLN.

      NR = 1;
      if ( L2ABER ) {
         // Standard absolute error bound suffices. All sigma_i with
         // sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
         // aggressive enforcement of lower numerical rank by introducing a
         // backward error of the order of N*EPSLN*||A||.
         TEMP1 = SQRT(REAL(N))*EPSLN;
         for (p = 2; p <= N; p++) { // 3001
            if ( ABS(A(p,p)) >= (TEMP1*ABS(A(1,1))) ) {
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
         TEMP1 = SQRT(SFMIN);
         for (p = 2; p <= N; p++) { // 3401
            IF ( ( ABS(A(p,p)) < (EPSLN*ABS(A(p-1,p-1))) ) || ( ABS(A(p,p)) < SMALL ) || ( L2KILL && (ABS(A(p,p)) < TEMP1) ) ) GO TO 3402;
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
         TEMP1  = SQRT(SFMIN);
         for (p = 2; p <= N; p++) { // 3301
            IF ( ( ABS(A(p,p)) < SMALL ) || ( L2KILL && (ABS(A(p,p)) < TEMP1) ) ) GO TO 3302;
            NR = NR + 1;
         } // 3301
         } // 3302

      }

      ALMORT = false;
      if ( NR == N ) {
         MAXPRJ = ONE;
         for (p = 2; p <= N; p++) { // 3051
            TEMP1  = ABS(A(p,p)) / SVA(IWORK(p));
            MAXPRJ = MIN( MAXPRJ, TEMP1 );
         } // 3051
         IF ( MAXPRJ**2 >= ONE - REAL(N)*EPSLN ) ALMORT = true;
      }


      SCONDA = - ONE;
      CONDR1 = - ONE;
      CONDR2 = - ONE;

      if ( ERREST ) {
         if ( N == NR ) {
            if ( RSVEC ) {
               // .. V is available as workspace
               clacpy('U', N, N, A, LDA, V, LDV );
               for (p = 1; p <= N; p++) { // 3053
                  TEMP1 = SVA(IWORK(p));
                  csscal(p, ONE/TEMP1, V(1,p), 1 );
               } // 3053
               if ( LSVEC ) {
                   cpocon('U', N, V, LDV, ONE, TEMP1, CWORK(N+1), RWORK, IERR );
               } else {
                   cpocon('U', N, V, LDV, ONE, TEMP1, CWORK, RWORK, IERR );
               }

            } else if ( LSVEC ) {
               // .. U is available as workspace
               clacpy('U', N, N, A, LDA, U, LDU );
               for (p = 1; p <= N; p++) { // 3054
                  TEMP1 = SVA(IWORK(p));
                  csscal(p, ONE/TEMP1, U(1,p), 1 );
               } // 3054
               cpocon('U', N, U, LDU, ONE, TEMP1, CWORK(N+1), RWORK, IERR );
            } else {
               clacpy('U', N, N, A, LDA, CWORK, N );
*[]            CALL CLACPY( 'U', N, N, A, LDA, CWORK(N+1), N )
               // Change: here index shifted by N to the left, CWORK(1:N)
               // not needed for SIGMA only computation
               for (p = 1; p <= N; p++) { // 3052
                  TEMP1 = SVA(IWORK(p));
*[]               CALL CSSCAL( p, ONE/TEMP1, CWORK(N+(p-1)*N+1), 1 )
                  csscal(p, ONE/TEMP1, CWORK((p-1)*N+1), 1 );
               } // 3052
            // .. the columns of R are scaled to have unit Euclidean lengths.
*[]               CALL CPOCON( 'U', N, CWORK(N+1), N, ONE, TEMP1,
*[]     $              CWORK(N+N*N+1), RWORK, IERR )
               cpocon('U', N, CWORK, N, ONE, TEMP1, CWORK(N*N+1), RWORK, IERR );

            }
            if ( TEMP1 != ZERO ) {
               SCONDA = ONE / SQRT(TEMP1);
            } else {
               SCONDA = - ONE;
            }
            // SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1).
            // N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
         } else {
            SCONDA = - ONE;
         }
      }

      L2PERT = L2PERT && ( ABS( A(1,1)/A(NR,NR) ) > SQRT(BIG1) );
      // If there is no violent scaling, artificial perturbation is not needed.

      // Phase 3:

      if ( !( RSVEC || LSVEC ) ) {

          // Singular Values only

          // .. transpose A(1:NR,1:N)
         DO 1946 p = 1, MIN( N-1, NR );
            ccopy(N-p, A(p,p+1), LDA, A(p+1,p), 1 );
            clacgv(N-p+1, A(p,p), 1 );
         } // 1946
         if (NR == N) A(N,N) = CONJG(A(N,N));

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
               // XSC = SQRT(SMALL)
               XSC = EPSLN / REAL(N);
               for (q = 1; q <= NR; q++) { // 4947
                  CTEMP = CMPLX(XSC*ABS(A(q,q)),ZERO);
                  for (p = 1; p <= N; p++) { // 4949
                     IF ( ( (p > q) && (ABS(A(p,q)) <= TEMP1) ) || ( p < q ) );
      // $                     A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) A(p,q) = CTEMP
                  } // 4949
               } // 4947
            } else {
               claset('U', NR-1,NR-1, CZERO,CZERO, A(1,2),LDA );
            }

             // .. second preconditioning using the QR factorization

            cgeqrf(N,NR, A,LDA, CWORK, CWORK(N+1),LWORK-N, IERR );

            // .. and transpose upper to lower triangular
            for (p = 1; p <= NR - 1; p++) { // 1948
               ccopy(NR-p, A(p,p+1), LDA, A(p+1,p), 1 );
               clacgv(NR-p+1, A(p,p), 1 );
            } // 1948

         }

            // Row-cyclic Jacobi SVD algorithm with column pivoting

            // .. again some perturbation (a "background noise") is added
            // to drown denormals
            if ( L2PERT ) {
               // XSC = SQRT(SMALL)
               XSC = EPSLN / REAL(N);
               for (q = 1; q <= NR; q++) { // 1947
                  CTEMP = CMPLX(XSC*ABS(A(q,q)),ZERO);
                  for (p = 1; p <= NR; p++) { // 1949
                     IF ( ( (p > q) && (ABS(A(p,q)) <= TEMP1) ) || ( p < q ) );
      // $                   A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) A(p,q) = CTEMP
                  } // 1949
               } // 1947
            } else {
               claset('U', NR-1, NR-1, CZERO, CZERO, A(1,2), LDA );
            }

            // .. and one-sided Jacobi rotations are started on a lower
            // triangular matrix (plus perturbation which is ignored in
            // the part which destroys triangular form (confusing?!))

            cgesvj('L', 'N', 'N', NR, NR, A, LDA, SVA, N, V, LDV, CWORK, LWORK, RWORK, LRWORK, INFO );

            SCALEM  = RWORK(1);
            NUMRANK = NINT(RWORK(2));


      } else if ( ( RSVEC && ( !LSVEC ) && ( !JRACC ) ) || ( JRACC && ( !LSVEC ) && ( NR != N ) ) ) {

         // -> Singular Values and Right Singular Vectors <-

         if ( ALMORT ) {

            // .. in this case NR equals N
            for (p = 1; p <= NR; p++) { // 1998
               ccopy(N-p+1, A(p,p), LDA, V(p,p), 1 );
               clacgv(N-p+1, V(p,p), 1 );
            } // 1998
            claset('U', NR-1,NR-1, CZERO, CZERO, V(1,2), LDV );

            cgesvj('L','U','N', N, NR, V, LDV, SVA, NR, A, LDA, CWORK, LWORK, RWORK, LRWORK, INFO );
            SCALEM  = RWORK(1);
            NUMRANK = NINT(RWORK(2));

         } else {

         // .. two more QR factorizations ( one QRF is not enough, two require
         // accumulated product of Jacobi rotations, three are perfect )

            claset('L', NR-1,NR-1, CZERO, CZERO, A(2,1), LDA );
            cgelqf(NR,N, A, LDA, CWORK, CWORK(N+1), LWORK-N, IERR);
            clacpy('L', NR, NR, A, LDA, V, LDV );
            claset('U', NR-1,NR-1, CZERO, CZERO, V(1,2), LDV );
            cgeqrf(NR, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR );
            for (p = 1; p <= NR; p++) { // 8998
               ccopy(NR-p+1, V(p,p), LDV, V(p,p), 1 );
               clacgv(NR-p+1, V(p,p), 1 );
            } // 8998
            claset('U', NR-1, NR-1, CZERO, CZERO, V(1,2), LDV);

            cgesvj('L', 'U','N', NR, NR, V,LDV, SVA, NR, U, LDU, CWORK(N+1), LWORK-N, RWORK, LRWORK, INFO );
            SCALEM  = RWORK(1);
            NUMRANK = NINT(RWORK(2));
            if ( NR < N ) {
               claset('A',N-NR, NR, CZERO,CZERO, V(NR+1,1),  LDV );
               claset('A',NR, N-NR, CZERO,CZERO, V(1,NR+1),  LDV );
               claset('A',N-NR,N-NR,CZERO,CONE, V(NR+1,NR+1),LDV );
            }

         cunmlq('L', 'C', N, N, NR, A, LDA, CWORK, V, LDV, CWORK(N+1), LWORK-N, IERR );

         }
          // .. permute the rows of V
          // DO 8991 p = 1, N
             // CALL CCOPY( N, V(p,1), LDV, A(IWORK(p),1), LDA )
* 8991    CONTINUE
          // CALL CLACPY( 'All', N, N, A, LDA, V, LDV )
         clapmr( false , N, N, V, LDV, IWORK );

          if ( TRANSP ) {
            clacpy('A', N, N, V, LDV, U, LDU );
          }

      } else if ( JRACC && ( !LSVEC) && ( NR == N ) ) {

         claset('L', N-1,N-1, CZERO, CZERO, A(2,1), LDA );

         cgesvj('U','N','V', N, N, A, LDA, SVA, N, V, LDV, CWORK, LWORK, RWORK, LRWORK, INFO );
          SCALEM  = RWORK(1);
          NUMRANK = NINT(RWORK(2));
          clapmr( false , N, N, V, LDV, IWORK );

      } else if ( LSVEC && ( !RSVEC ) ) {

         // .. Singular Values and Left Singular Vectors                 ..

         // .. second preconditioning step to avoid need to accumulate
         // Jacobi rotations in the Jacobi iterations.
         for (p = 1; p <= NR; p++) { // 1965
            ccopy(N-p+1, A(p,p), LDA, U(p,p), 1 );
            clacgv(N-p+1, U(p,p), 1 );
         } // 1965
         claset('U', NR-1, NR-1, CZERO, CZERO, U(1,2), LDU );

         cgeqrf(N, NR, U, LDU, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR );

         for (p = 1; p <= NR - 1; p++) { // 1967
            ccopy(NR-p, U(p,p+1), LDU, U(p+1,p), 1 );
            clacgv(N-p+1, U(p,p), 1 );
         } // 1967
         claset('U', NR-1, NR-1, CZERO, CZERO, U(1,2), LDU );

         cgesvj('L', 'U', 'N', NR,NR, U, LDU, SVA, NR, A, LDA, CWORK(N+1), LWORK-N, RWORK, LRWORK, INFO );
         SCALEM  = RWORK(1);
         NUMRANK = NINT(RWORK(2));

         if ( NR < M ) {
            claset('A',  M-NR, NR,CZERO, CZERO, U(NR+1,1), LDU );
            if ( NR < N1 ) {
               claset('A',NR, N1-NR, CZERO, CZERO, U(1,NR+1),LDU );
               claset('A',M-NR,N1-NR,CZERO,CONE,U(NR+1,NR+1),LDU );
            }
         }

         cunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR );

         if (ROWPIV) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 );

         for (p = 1; p <= N1; p++) { // 1974
            XSC = ONE / SCNRM2( M, U(1,p), 1 );
            csscal(M, XSC, U(1,p), 1 );
         } // 1974

         if ( TRANSP ) {
            clacpy('A', N, N, U, LDU, V, LDV );
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
            // optimized implementation of CGEJSV.

            for (p = 1; p <= NR; p++) { // 1968
               ccopy(N-p+1, A(p,p), LDA, V(p,p), 1 );
               clacgv(N-p+1, V(p,p), 1 );
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
               XSC = SQRT(SMALL);
               for (q = 1; q <= NR; q++) { // 2969
                  CTEMP = CMPLX(XSC*ABS( V(q,q) ),ZERO);
                  for (p = 1; p <= N; p++) { // 2968
                     IF ( ( p > q ) && ( ABS(V(p,q)) <= TEMP1 ) || ( p < q ) );
      // $                   V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) V(p,q) = CTEMP
                     if (p < q) V(p,q) = - V(p,q);
                  } // 2968
               } // 2969
            } else {
               claset('U', NR-1, NR-1, CZERO, CZERO, V(1,2), LDV );
            }

            // Estimate the row scaled condition number of R1
            // (If R1 is rectangular, N > NR, then the condition number
            // of the leading NR x NR submatrix is estimated.)

            clacpy('L', NR, NR, V, LDV, CWORK(2*N+1), NR );
            for (p = 1; p <= NR; p++) { // 3950
               TEMP1 = SCNRM2(NR-p+1,CWORK(2*N+(p-1)*NR+p),1);
               csscal(NR-p+1,ONE/TEMP1,CWORK(2*N+(p-1)*NR+p),1);
            } // 3950
            cpocon('L',NR,CWORK(2*N+1),NR,ONE,TEMP1, CWORK(2*N+NR*NR+1),RWORK,IERR);
            CONDR1 = ONE / SQRT(TEMP1);
            // .. here need a second opinion on the condition number
            // .. then assume worst case scenario
            // R1 is OK for inverse <=> CONDR1 < REAL(N)
            // more conservative    <=> CONDR1 < SQRT(REAL(N))

            COND_OK = SQRT(SQRT(REAL(NR)));
*[TP]       COND_OK is a tuning parameter.

            if ( CONDR1 < COND_OK ) {
               // .. the second QRF without pivoting. Note: in an optimized
               // implementation, this QRF should be implemented as the QRF
               // of a lower triangular matrix.
               // R1^* = Q2 * R2
               cgeqrf(N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR );

               if ( L2PERT ) {
                  XSC = SQRT(SMALL)/EPSLN;
                  for (p = 2; p <= NR; p++) { // 3959
                     for (q = 1; q <= p - 1; q++) { // 3958
                        CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), ZERO);
                        IF ( ABS(V(q,p)) <= TEMP1 );
      // $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) V(q,p) = CTEMP
                     } // 3958
                  } // 3959
               }

               if (NR != N) CALL CLACPY( 'A', N, NR, V, LDV, CWORK(2*N+1), N );
               // .. save ...

            // .. this transposed copy should be better than naive
               for (p = 1; p <= NR - 1; p++) { // 1969
                  ccopy(NR-p, V(p,p+1), LDV, V(p+1,p), 1 );
                  clacgv(NR-p+1, V(p,p), 1 );
               } // 1969
               V(NR,NR)=CONJG(V(NR,NR));

               CONDR2 = CONDR1;

            } else {

               // .. ill-conditioned case: second QRF with pivoting
               // Note that windowed pivoting would be equally good
               // numerically, and more run-time efficient. So, in
               // an optimal implementation, the next call to CGEQP3
               // should be replaced with eg. CALL CGEQPX (ACM TOMS #782)
               // with properly (carefully) chosen parameters.

               // R1^* * P2 = Q2 * R2
               for (p = 1; p <= NR; p++) { // 3003
                  IWORK(N+p) = 0;
               } // 3003
               cgeqp3(N, NR, V, LDV, IWORK(N+1), CWORK(N+1), CWORK(2*N+1), LWORK-2*N, RWORK, IERR );
**               CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1),
**     $              LWORK-2*N, IERR )
               if ( L2PERT ) {
                  XSC = SQRT(SMALL);
                  for (p = 2; p <= NR; p++) { // 3969
                     for (q = 1; q <= p - 1; q++) { // 3968
                        CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), ZERO);
                        IF ( ABS(V(q,p)) <= TEMP1 );
      // $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) V(q,p) = CTEMP
                     } // 3968
                  } // 3969
               }

               clacpy('A', N, NR, V, LDV, CWORK(2*N+1), N );

               if ( L2PERT ) {
                  XSC = SQRT(SMALL);
                  for (p = 2; p <= NR; p++) { // 8970
                     for (q = 1; q <= p - 1; q++) { // 8971
                        CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), ZERO);
                         // V(p,q) = - TEMP1*( V(q,p) / ABS(V(q,p)) )
                        V(p,q) = - CTEMP;
                     } // 8971
                  } // 8970
               } else {
                  claset('L',NR-1,NR-1,CZERO,CZERO,V(2,1),LDV );
               }
               // Now, compute R2 = L3 * Q3, the LQ factorization.
               cgelqf(NR, NR, V, LDV, CWORK(2*N+N*NR+1), CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR );
               // .. and estimate the condition number
               clacpy('L',NR,NR,V,LDV,CWORK(2*N+N*NR+NR+1),NR );
               for (p = 1; p <= NR; p++) { // 4950
                  TEMP1 = SCNRM2( p, CWORK(2*N+N*NR+NR+p), NR );
                  csscal(p, ONE/TEMP1, CWORK(2*N+N*NR+NR+p), NR );
               } // 4950
               cpocon('L',NR,CWORK(2*N+N*NR+NR+1),NR,ONE,TEMP1, CWORK(2*N+N*NR+NR+NR*NR+1),RWORK,IERR );
               CONDR2 = ONE / SQRT(TEMP1);


               if ( CONDR2 >= COND_OK ) {
                  // .. save the Householder vectors used for Q3
                  // (this overwrites the copy of R2, as it will not be
                  // needed in this branch, but it does not overwrite the
                  // Huseholder vectors of Q2.).
                  clacpy('U', NR, NR, V, LDV, CWORK(2*N+1), N );
                  // .. and the rest of the information on Q3 is in
                  // WORK(2*N+N*NR+1:2*N+N*NR+N)
               }

            }

            if ( L2PERT ) {
               XSC = SQRT(SMALL);
               for (q = 2; q <= NR; q++) { // 4968
                  CTEMP = XSC * V(q,q);
                  for (p = 1; p <= q - 1; p++) { // 4969
                      // V(p,q) = - TEMP1*( V(p,q) / ABS(V(p,q)) )
                     V(p,q) = - CTEMP;
                  } // 4969
               } // 4968
            } else {
               claset('U', NR-1,NR-1, CZERO,CZERO, V(1,2), LDV );
            }

         // Second preconditioning finished; continue with Jacobi SVD
         // The input matrix is lower triangular.

         // Recover the right singular vectors as solution of a well
         // conditioned triangular matrix equation.

            if ( CONDR1 < COND_OK ) {

               cgesvj('L','U','N',NR,NR,V,LDV,SVA,NR,U, LDU, CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,RWORK, LRWORK, INFO );
               SCALEM  = RWORK(1);
               NUMRANK = NINT(RWORK(2));
               for (p = 1; p <= NR; p++) { // 3970
                  ccopy(NR, V(1,p), 1, U(1,p), 1 );
                  csscal(NR, SVA(p),    V(1,p), 1 );
               } // 3970

         // .. pick the right matrix equation and solve it

               if ( NR == N ) {
* :))             .. best case, R1 is inverted. The solution of this matrix
                  // equation is Q2*V2 = the product of the Jacobi rotations
                  // used in CGESVJ, premultiplied with the orthogonal matrix
                  // from the second QR factorization.
                  ctrsm('L','U','N','N', NR,NR,CONE, A,LDA, V,LDV);
               } else {
                  // .. R1 is well conditioned, but non-square. Adjoint of R2
                  // is inverted to get the product of the Jacobi rotations
                  // used in CGESVJ. The Q-factor from the second QR
                  // factorization is then built in explicitly.
                  ctrsm('L','U','C','N',NR,NR,CONE,CWORK(2*N+1), N,V,LDV);
                  if ( NR < N ) {
                  claset('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV);
                  claset('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV);
                  claset('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV);
                  }
                  cunmqr('L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR);
               }

            } else if ( CONDR2 < COND_OK ) {

               // The matrix R2 is inverted. The solution of the matrix equation
               // is Q3^* * V3 = the product of the Jacobi rotations (applied to
               // the lower triangular L3 from the LQ factorization of
               // R2=L3*Q3), pre-multiplied with the transposed Q3.
               cgesvj('L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, RWORK, LRWORK, INFO );
               SCALEM  = RWORK(1);
               NUMRANK = NINT(RWORK(2));
               for (p = 1; p <= NR; p++) { // 3870
                  ccopy(NR, V(1,p), 1, U(1,p), 1 );
                  csscal(NR, SVA(p),    U(1,p), 1 );
               } // 3870
               ctrsm('L','U','N','N',NR,NR,CONE,CWORK(2*N+1),N, U,LDU);
               // .. apply the permutation from the second QR factorization
               for (q = 1; q <= NR; q++) { // 873
                  for (p = 1; p <= NR; p++) { // 872
                     CWORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q);
                  } // 872
                  for (p = 1; p <= NR; p++) { // 874
                     U(p,q) = CWORK(2*N+N*NR+NR+p);
                  } // 874
               } // 873
               if ( NR < N ) {
                  claset('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV );
                  claset('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV );
                  claset('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV);
               }
               cunmqr('L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR );
            } else {
               // Last line of defense.
* #:(          This is a rather pathological case: no scaled condition
               // improvement after two pivoted QR factorizations. Other
               // possibility is that the rank revealing QR factorization
               // or the condition estimator has failed, or the COND_OK
               // is set very close to ONE (which is unnecessary). Normally,
               // this branch should never be executed, but in rare cases of
               // failure of the RRQR or condition estimator, the last line of
               // defense ensures that CGEJSV completes the task.
               // Compute the full SVD of L3 using CGESVJ with explicit
               // accumulation of Jacobi rotations.
               cgesvj('L', 'U', 'V', NR, NR, V, LDV, SVA, NR, U, LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, RWORK, LRWORK, INFO );
               SCALEM  = RWORK(1);
               NUMRANK = NINT(RWORK(2));
               if ( NR < N ) {
                  claset('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV );
                  claset('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV );
                  claset('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV);
               }
               cunmqr('L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR );

               cunmlq('L', 'C', NR, NR, NR, CWORK(2*N+1), N, CWORK(2*N+N*NR+1), U, LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR );
               for (q = 1; q <= NR; q++) { // 773
                  for (p = 1; p <= NR; p++) { // 772
                     CWORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q);
                  } // 772
                  for (p = 1; p <= NR; p++) { // 774
                     U(p,q) = CWORK(2*N+N*NR+NR+p);
                  } // 774
               } // 773

            }

            // Permute the rows of V using the (column) permutation from the
            // first QRF. Also, scale the columns to make them unit in
            // Euclidean norm. This applies to all cases.

            TEMP1 = SQRT(REAL(N)) * EPSLN;
            for (q = 1; q <= N; q++) { // 1972
               for (p = 1; p <= N; p++) { // 972
                  CWORK(2*N+N*NR+NR+IWORK(p)) = V(p,q);
               } // 972
               for (p = 1; p <= N; p++) { // 973
                  V(p,q) = CWORK(2*N+N*NR+NR+p);
               } // 973
               XSC = ONE / SCNRM2( N, V(1,q), 1 );
               IF ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) CALL CSSCAL( N, XSC, V(1,q), 1 );
            } // 1972
            // At this moment, V contains the right singular vectors of A.
            // Next, assemble the left singular vector matrix U (M x N).
            if ( NR < M ) {
               claset('A', M-NR, NR, CZERO, CZERO, U(NR+1,1), LDU);
               if ( NR < N1 ) {
                  claset('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU);
                  claset('A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1),LDU);
               }
            }

            // The Q matrix from the first QRF is built into the left singular
            // matrix U. This applies to all cases.

            cunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR );

            // The columns of U are normalized. The cost is O(M*N) flops.
            TEMP1 = SQRT(REAL(M)) * EPSLN;
            for (p = 1; p <= NR; p++) { // 1973
               XSC = ONE / SCNRM2( M, U(1,p), 1 );
               IF ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) CALL CSSCAL( M, XSC, U(1,p), 1 );
            } // 1973

            // If the initial QRF is computed with row pivoting, the left
            // singular vectors must be adjusted.

            if (ROWPIV) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 );

         } else {

         // .. the initial matrix A has almost orthogonal columns and
         // the second QRF is not needed

            clacpy('U', N, N, A, LDA, CWORK(N+1), N );
            if ( L2PERT ) {
               XSC = SQRT(SMALL);
               for (p = 2; p <= N; p++) { // 5970
                  CTEMP = XSC * CWORK( N + (p-1)*N + p );
                  for (q = 1; q <= p - 1; q++) { // 5971
                      // CWORK(N+(q-1)*N+p)=-TEMP1 * ( CWORK(N+(p-1)*N+q) /
      // $                                        ABS(CWORK(N+(p-1)*N+q)) )
                     CWORK(N+(q-1)*N+p)=-CTEMP;
                  } // 5971
               } // 5970
            } else {
               claset('L',N-1,N-1,CZERO,CZERO,CWORK(N+2),N );
            }

            cgesvj('U', 'U', 'N', N, N, CWORK(N+1), N, SVA, N, U, LDU, CWORK(N+N*N+1), LWORK-N-N*N, RWORK, LRWORK, INFO );

            SCALEM  = RWORK(1);
            NUMRANK = NINT(RWORK(2));
            for (p = 1; p <= N; p++) { // 6970
               ccopy(N, CWORK(N+(p-1)*N+1), 1, U(1,p), 1 );
               csscal(N, SVA(p), CWORK(N+(p-1)*N+1), 1 );
            } // 6970

            ctrsm('L', 'U', 'N', 'N', N, N, CONE, A, LDA, CWORK(N+1), N );
            for (p = 1; p <= N; p++) { // 6972
               ccopy(N, CWORK(N+p), N, V(IWORK(p),1), LDV );
            } // 6972
            TEMP1 = SQRT(REAL(N))*EPSLN;
            for (p = 1; p <= N; p++) { // 6971
               XSC = ONE / SCNRM2( N, V(1,p), 1 );
               IF ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) CALL CSSCAL( N, XSC, V(1,p), 1 );
            } // 6971

            // Assemble the left singular vector matrix U (M x N).

            if ( N < M ) {
               claset('A',  M-N, N, CZERO, CZERO, U(N+1,1), LDU );
               if ( N < N1 ) {
                  claset('A',N,  N1-N, CZERO, CZERO,  U(1,N+1),LDU);
                  claset('A',M-N,N1-N, CZERO, CONE,U(N+1,N+1),LDU);
               }
            }
            cunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR );
            TEMP1 = SQRT(REAL(M))*EPSLN;
            for (p = 1; p <= N1; p++) { // 6973
               XSC = ONE / SCNRM2( M, U(1,p), 1 );
               IF ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) CALL CSSCAL( M, XSC, U(1,p), 1 );
            } // 6973

            if (ROWPIV) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 );

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
         // in presence of extreme values, e.g. when the singular values spread from
         // the underflow to the overflow threshold.

         for (p = 1; p <= NR; p++) { // 7968
            ccopy(N-p+1, A(p,p), LDA, V(p,p), 1 );
            clacgv(N-p+1, V(p,p), 1 );
         } // 7968

         if ( L2PERT ) {
            XSC = SQRT(SMALL/EPSLN);
            for (q = 1; q <= NR; q++) { // 5969
               CTEMP = CMPLX(XSC*ABS( V(q,q) ),ZERO);
               for (p = 1; p <= N; p++) { // 5968
                  IF ( ( p > q ) && ( ABS(V(p,q)) <= TEMP1 ) || ( p < q ) );
      // $                V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) V(p,q) = CTEMP
                  if (p < q) V(p,q) = - V(p,q);
               } // 5968
            } // 5969
         } else {
            claset('U', NR-1, NR-1, CZERO, CZERO, V(1,2), LDV );
         }
          cgeqrf(N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR );
         clacpy('L', N, NR, V, LDV, CWORK(2*N+1), N );

         for (p = 1; p <= NR; p++) { // 7969
            ccopy(NR-p+1, V(p,p), LDV, U(p,p), 1 );
            clacgv(NR-p+1, U(p,p), 1 );
         } // 7969

         if ( L2PERT ) {
            XSC = SQRT(SMALL/EPSLN);
            for (q = 2; q <= NR; q++) { // 9970
               for (p = 1; p <= q - 1; p++) { // 9971
                  CTEMP = CMPLX(XSC * MIN(ABS(U(p,p)),ABS(U(q,q))), ZERO);
                   // U(p,q) = - TEMP1 * ( U(q,p) / ABS(U(q,p)) )
                  U(p,q) = - CTEMP;
               } // 9971
            } // 9970
         } else {
            claset('U', NR-1, NR-1, CZERO, CZERO, U(1,2), LDU );
         }
          cgesvj('L', 'U', 'V', NR, NR, U, LDU, SVA, N, V, LDV, CWORK(2*N+N*NR+1), LWORK-2*N-N*NR, RWORK, LRWORK, INFO );
         SCALEM  = RWORK(1);
         NUMRANK = NINT(RWORK(2));

         if ( NR < N ) {
            claset('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV );
            claset('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV );
            claset('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV );
         }
          cunmqr('L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR );

            // Permute the rows of V using the (column) permutation from the
            // first QRF. Also, scale the columns to make them unit in
            // Euclidean norm. This applies to all cases.

            TEMP1 = SQRT(REAL(N)) * EPSLN;
            for (q = 1; q <= N; q++) { // 7972
               for (p = 1; p <= N; p++) { // 8972
                  CWORK(2*N+N*NR+NR+IWORK(p)) = V(p,q);
               } // 8972
               for (p = 1; p <= N; p++) { // 8973
                  V(p,q) = CWORK(2*N+N*NR+NR+p);
               } // 8973
               XSC = ONE / SCNRM2( N, V(1,q), 1 );
               IF ( (XSC < (ONE-TEMP1)) || (XSC > (ONE+TEMP1)) ) CALL CSSCAL( N, XSC, V(1,q), 1 );
            } // 7972

            // At this moment, V contains the right singular vectors of A.
            // Next, assemble the left singular vector matrix U (M x N).

         if ( NR < M ) {
            claset('A',  M-NR, NR, CZERO, CZERO, U(NR+1,1), LDU );
            if ( NR < N1 ) {
               claset('A',NR,  N1-NR, CZERO, CZERO,  U(1,NR+1),LDU);
               claset('A',M-NR,N1-NR, CZERO, CONE,U(NR+1,NR+1),LDU);
            }
         }

         cunmqr('L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR );

            if (ROWPIV) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 );


         }
         if ( TRANSP ) {
            // .. swap U and V because the procedure worked on A^*
            for (p = 1; p <= N; p++) { // 6974
               cswap(N, U(1,p), 1, V(1,p), 1 );
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
            SVA(p) = ZERO;
         } // 3004
      }

      RWORK(1) = USCAL2 * SCALEM;
      RWORK(2) = USCAL1;
      if (ERREST) RWORK(3) = SCONDA;
      if ( LSVEC && RSVEC ) {
         RWORK(4) = CONDR1;
         RWORK(5) = CONDR2;
      }
      if ( L2TRAN ) {
         RWORK(6) = ENTRA;
         RWORK(7) = ENTRAT;
      }

      IWORK(1) = NR;
      IWORK(2) = NUMRANK;
      IWORK(3) = WARNING;
      if ( TRANSP ) {
          IWORK(4) =  1;
      } else {
          IWORK(4) = -1;
      }


      RETURN;
      // ..
      // .. END OF CGEJSV
      // ..
      }
