      SUBROUTINE CGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, M, N, A, LDA, SVA, U, LDU, V, LDV, CWORK, LWORK, RWORK, LRWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      IMPLICIT    NONE
      int         INFO, LDA, LDU, LDV, LWORK, LRWORK, M, N;
*     ..
*     .. Array Arguments ..
      COMPLEX     A( LDA, * ), U( LDU, * ), V( LDV, * ), CWORK( LWORK )
      REAL        SVA( N ), RWORK( LRWORK )
      int         IWORK( * );
      String      JOBA, JOBP, JOBR, JOBT, JOBU, JOBV;
*     ..
*
*  ===========================================================================
*
*     .. Local Parameters ..
      REAL        ZERO,         ONE
      PARAMETER ( ZERO = 0.0E0, ONE = 1.0E0 )
      COMPLEX     CZERO,                    CONE
      PARAMETER ( CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX CTEMP
      REAL    AAPP,   AAQQ,   AATMAX, AATMIN, BIG,    BIG1,   COND_OK, CONDR1, CONDR2, ENTRA,  ENTRAT, EPSLN,  MAXPRJ, SCALEM, SCONDA, SFMIN,  SMALL,  TEMP1,  USCAL1, USCAL2, XSC
      int     IERR,   N1,     NR,     NUMRANK,        p, q,   WARNING;
      bool    ALMORT, DEFR,   ERREST, GOSCAL,  JRACC,  KILL,   LQUERY, LSVEC,  L2ABER, L2KILL, L2PERT,  L2RANK, L2TRAN, NOSCAL, ROWPIV, RSVEC,  TRANSP;
*
      int     OPTWRK, MINWRK, MINRWRK, MINIWRK;
      int     LWCON,  LWLQF, LWQP3, LWQRF, LWUNMLQ, LWUNMQR, LWUNMQRM, LWSVDJ, LWSVDJV, LRWQP3, LRWCON, LRWSVDJ, IWOFF       int     LWRK_CGELQF, LWRK_CGEQP3,  LWRK_CGEQP3N, LWRK_CGEQRF,   LWRK_CGESVJ, LWRK_CGESVJV, LWRK_CGESVJU, LWRK_CUNMLQ, LWRK_CUNMQR, LWRK_CUNMQRM;
*     ..
*     .. Local Arrays
      COMPLEX CDUMMY(1)
      REAL    RDUMMY(1)
*
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG, ALOG, MAX, MIN, REAL, NINT, SQRT
*     ..
*     .. External Functions ..
      REAL      SLAMCH, SCNRM2
      int       ISAMAX, ICAMAX;
      bool      LSAME;
      // EXTERNAL ISAMAX, ICAMAX, LSAME, SLAMCH, SCNRM2
*     ..
*     .. External Subroutines ..
      // EXTERNAL SLASSQ, CCOPY,  CGELQF, CGEQP3, CGEQRF, CLACPY, CLAPMR, CLASCL, SLASCL, CLASET, CLASSQ, CLASWP, CUNGQR, CUNMLQ, CUNMQR, CPOCON, SSCAL,  CSSCAL, CSWAP,  CTRSM,  CLACGV, XERBLA
*
      // EXTERNAL CGESVJ
*     ..
*
*     Test the input arguments
*
      LSVEC  = LSAME( JOBU, 'U' ) .OR. LSAME( JOBU, 'F' )
      JRACC  = LSAME( JOBV, 'J' )
      RSVEC  = LSAME( JOBV, 'V' ) .OR. JRACC
      ROWPIV = LSAME( JOBA, 'F' ) .OR. LSAME( JOBA, 'G' )
      L2RANK = LSAME( JOBA, 'R' )
      L2ABER = LSAME( JOBA, 'A' )
      ERREST = LSAME( JOBA, 'E' ) .OR. LSAME( JOBA, 'G' )
      L2TRAN = LSAME( JOBT, 'T' ) .AND. ( M .EQ. N )
      L2KILL = LSAME( JOBR, 'R' )
      DEFR   = LSAME( JOBR, 'N' )
      L2PERT = LSAME( JOBP, 'P' )
*
      LQUERY = ( LWORK .EQ. -1 ) .OR. ( LRWORK .EQ. -1 )
*
      IF ( .NOT.(ROWPIV .OR. L2RANK .OR. L2ABER .OR. ERREST .OR. LSAME( JOBA, 'C' ) )) THEN
         INFO = - 1
      ELSE IF ( .NOT.( LSVEC .OR. LSAME( JOBU, 'N' ) .OR. ( LSAME( JOBU, 'W' ) .AND. RSVEC .AND. L2TRAN ) ) ) THEN
         INFO = - 2
      ELSE IF ( .NOT.( RSVEC .OR. LSAME( JOBV, 'N' ) .OR. ( LSAME( JOBV, 'W' ) .AND. LSVEC .AND. L2TRAN ) ) ) THEN
         INFO = - 3
      ELSE IF ( .NOT. ( L2KILL .OR. DEFR ) )    THEN
         INFO = - 4
      ELSE IF ( .NOT. ( LSAME(JOBT,'T') .OR. LSAME(JOBT,'N') ) ) THEN
         INFO = - 5
      ELSE IF ( .NOT. ( L2PERT .OR. LSAME( JOBP, 'N' ) ) ) THEN
         INFO = - 6
      ELSE IF ( M .LT. 0 ) THEN
         INFO = - 7
      ELSE IF ( ( N .LT. 0 ) .OR. ( N .GT. M ) ) THEN
         INFO = - 8
      ELSE IF ( LDA .LT. M ) THEN
         INFO = - 10
      ELSE IF ( LSVEC .AND. ( LDU .LT. M ) ) THEN
         INFO = - 13
      ELSE IF ( RSVEC .AND. ( LDV .LT. N ) ) THEN
         INFO = - 15
      ELSE
*        #:)
         INFO = 0
      END IF
*
      IF ( INFO .EQ. 0 ) THEN
*         .. compute the minimal and the optimal workspace lengths
*         [[The expressions for computing the minimal and the optimal
*         values of LCWORK, LRWORK are written with a lot of redundancy and
*         can be simplified. However, this verbose form is useful for
*         maintenance and modifications of the code.]]
*
*        .. minimal workspace length for CGEQP3 of an M x N matrix,
*         CGEQRF of an N x N matrix, CGELQF of an N x N matrix,
*         CUNMLQ for computing N x N matrix, CUNMQR for computing N x N
*         matrix, CUNMQR for computing M x N matrix, respectively.
          LWQP3 = N+1
          LWQRF = MAX( 1, N )
          LWLQF = MAX( 1, N )
          LWUNMLQ  = MAX( 1, N )
          LWUNMQR  = MAX( 1, N )
          LWUNMQRM = MAX( 1, M )
*        .. minimal workspace length for CPOCON of an N x N matrix
          LWCON = 2 * N
*        .. minimal workspace length for CGESVJ of an N x N matrix,
*         without and with explicit accumulation of Jacobi rotations
          LWSVDJ  = MAX( 2 * N, 1 )
          LWSVDJV = MAX( 2 * N, 1 )
*         .. minimal REAL workspace length for CGEQP3, CPOCON, CGESVJ
          LRWQP3  = 2 * N
          LRWCON  = N
          LRWSVDJ = N
          IF ( LQUERY ) THEN
              CALL CGEQP3( M, N, A, LDA, IWORK, CDUMMY, CDUMMY, -1,  RDUMMY, IERR )
              LWRK_CGEQP3 = INT( CDUMMY(1) )
              CALL CGEQRF( N, N, A, LDA, CDUMMY, CDUMMY,-1, IERR )
              LWRK_CGEQRF = INT( CDUMMY(1) )
              CALL CGELQF( N, N, A, LDA, CDUMMY, CDUMMY,-1, IERR )
              LWRK_CGELQF = INT( CDUMMY(1) )
          END IF
          MINWRK  = 2
          OPTWRK  = 2
          MINIWRK = N
          IF ( .NOT. (LSVEC .OR. RSVEC ) ) THEN
*             .. minimal and optimal sizes of the complex workspace if
*             only the singular values are requested
              IF ( ERREST ) THEN
                  MINWRK = MAX( N+LWQP3, N**2+LWCON, N+LWQRF, LWSVDJ )
              ELSE
                  MINWRK = MAX( N+LWQP3, N+LWQRF, LWSVDJ )
              END IF
              IF ( LQUERY ) THEN
                  CALL CGESVJ( 'L', 'N', 'N', N, N, A, LDA, SVA, N, V,  LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                  LWRK_CGESVJ = INT( CDUMMY(1) )
                  IF ( ERREST ) THEN
                      OPTWRK = MAX( N+LWRK_CGEQP3, N**2+LWCON,  N+LWRK_CGEQRF, LWRK_CGESVJ )
                  ELSE
                      OPTWRK = MAX( N+LWRK_CGEQP3, N+LWRK_CGEQRF,  LWRK_CGESVJ )
                  END IF
              END IF
              IF ( L2TRAN .OR. ROWPIV ) THEN
                  IF ( ERREST ) THEN
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWCON, LRWSVDJ )
                  ELSE
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ )
                  END IF
              ELSE
                  IF ( ERREST ) THEN
                     MINRWRK = MAX( 7, LRWQP3, LRWCON, LRWSVDJ )
                  ELSE
                     MINRWRK = MAX( 7, LRWQP3, LRWSVDJ )
                  END IF
              END IF
              IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
          ELSE IF ( RSVEC .AND. (.NOT.LSVEC) ) THEN
*            .. minimal and optimal sizes of the complex workspace if the
*            singular values and the right singular vectors are requested
             IF ( ERREST ) THEN
                 MINWRK = MAX( N+LWQP3, LWCON, LWSVDJ, N+LWLQF,   2*N+LWQRF, N+LWSVDJ, N+LWUNMLQ )
             ELSE
                 MINWRK = MAX( N+LWQP3, LWSVDJ, N+LWLQF, 2*N+LWQRF,  N+LWSVDJ, N+LWUNMLQ )
             END IF
             IF ( LQUERY ) THEN
                 CALL CGESVJ( 'L', 'U', 'N', N,N, U, LDU, SVA, N, A, LDA, CDUMMY, -1, RDUMMY, -1, IERR )
                 LWRK_CGESVJ = INT( CDUMMY(1) )
                 CALL CUNMLQ( 'L', 'C', N, N, N, A, LDA, CDUMMY, V, LDV, CDUMMY, -1, IERR )
                 LWRK_CUNMLQ = INT( CDUMMY(1) )
                 IF ( ERREST ) THEN
                 OPTWRK = MAX( N+LWRK_CGEQP3, LWCON, LWRK_CGESVJ,  N+LWRK_CGELQF, 2*N+LWRK_CGEQRF, N+LWRK_CGESVJ,  N+LWRK_CUNMLQ )
                 ELSE
                 OPTWRK = MAX( N+LWRK_CGEQP3, LWRK_CGESVJ,N+LWRK_CGELQF, 2*N+LWRK_CGEQRF, N+LWRK_CGESVJ, N+LWRK_CUNMLQ )
                 END IF
             END IF
             IF ( L2TRAN .OR. ROWPIV ) THEN
                  IF ( ERREST ) THEN
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON )
                  ELSE
                     MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ )
                  END IF
             ELSE
                  IF ( ERREST ) THEN
                     MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON )
                  ELSE
                     MINRWRK = MAX( 7, LRWQP3, LRWSVDJ )
                  END IF
             END IF
             IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
          ELSE IF ( LSVEC .AND. (.NOT.RSVEC) ) THEN
*            .. minimal and optimal sizes of the complex workspace if the
*            singular values and the left singular vectors are requested
             IF ( ERREST ) THEN
                 MINWRK = N + MAX( LWQP3,LWCON,N+LWQRF,LWSVDJ,LWUNMQRM )
             ELSE
                 MINWRK = N + MAX( LWQP3, N+LWQRF, LWSVDJ, LWUNMQRM )
             END IF
             IF ( LQUERY ) THEN
                 CALL CGESVJ( 'L', 'U', 'N', N,N, U, LDU, SVA, N, A, LDA, CDUMMY, -1, RDUMMY, -1, IERR )
                 LWRK_CGESVJ = INT( CDUMMY(1) )
                 CALL CUNMQR( 'L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR )
                 LWRK_CUNMQRM = INT( CDUMMY(1) )
                 IF ( ERREST ) THEN
                 OPTWRK = N + MAX( LWRK_CGEQP3, LWCON, N+LWRK_CGEQRF, LWRK_CGESVJ, LWRK_CUNMQRM )
                 ELSE
                 OPTWRK = N + MAX( LWRK_CGEQP3, N+LWRK_CGEQRF, LWRK_CGESVJ, LWRK_CUNMQRM )
                 END IF
             END IF
             IF ( L2TRAN .OR. ROWPIV ) THEN
                 IF ( ERREST ) THEN
                    MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON )
                 ELSE
                    MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ )
                 END IF
             ELSE
                 IF ( ERREST ) THEN
                    MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON )
                 ELSE
                    MINRWRK = MAX( 7, LRWQP3, LRWSVDJ )
                 END IF
             END IF
             IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
          ELSE
*            .. minimal and optimal sizes of the complex workspace if the
*            full SVD is requested
             IF ( .NOT. JRACC ) THEN
                 IF ( ERREST ) THEN
                    MINWRK = MAX( N+LWQP3, N+LWCON,  2*N+N**2+LWCON,  2*N+LWQRF,         2*N+LWQP3, 2*N+N**2+N+LWLQF,  2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWSVDJ, 2*N+N**2+N+LWSVDJV, 2*N+N**2+N+LWUNMQR,2*N+N**2+N+LWUNMLQ, N+N**2+LWSVDJ,   N+LWUNMQRM )
                 ELSE
                    MINWRK = MAX( N+LWQP3,        2*N+N**2+LWCON,  2*N+LWQRF,         2*N+LWQP3, 2*N+N**2+N+LWLQF,  2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWSVDJ, 2*N+N**2+N+LWSVDJV, 2*N+N**2+N+LWUNMQR,2*N+N**2+N+LWUNMLQ, N+N**2+LWSVDJ,      N+LWUNMQRM )
                 END IF
                 MINIWRK = MINIWRK + N
                 IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
             ELSE
                 IF ( ERREST ) THEN
                    MINWRK = MAX( N+LWQP3, N+LWCON, 2*N+LWQRF,  2*N+N**2+LWSVDJV, 2*N+N**2+N+LWUNMQR, N+LWUNMQRM )
                 ELSE
                    MINWRK = MAX( N+LWQP3, 2*N+LWQRF,  2*N+N**2+LWSVDJV, 2*N+N**2+N+LWUNMQR, N+LWUNMQRM )
                 END IF
                 IF ( ROWPIV .OR. L2TRAN ) MINIWRK = MINIWRK + M
             END IF
             IF ( LQUERY ) THEN
                 CALL CUNMQR( 'L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR )
                 LWRK_CUNMQRM = INT( CDUMMY(1) )
                 CALL CUNMQR( 'L', 'N', N, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR )
                 LWRK_CUNMQR = INT( CDUMMY(1) )
                 IF ( .NOT. JRACC ) THEN
                     CALL CGEQP3( N,N, A, LDA, IWORK, CDUMMY,CDUMMY, -1, RDUMMY, IERR )
                     LWRK_CGEQP3N = INT( CDUMMY(1) )
                     CALL CGESVJ( 'L', 'U', 'N', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                     LWRK_CGESVJ = INT( CDUMMY(1) )
                     CALL CGESVJ( 'U', 'U', 'N', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                     LWRK_CGESVJU = INT( CDUMMY(1) )
                     CALL CGESVJ( 'L', 'U', 'V', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                     LWRK_CGESVJV = INT( CDUMMY(1) )
                     CALL CUNMLQ( 'L', 'C', N, N, N, A, LDA, CDUMMY, V, LDV, CDUMMY, -1, IERR )
                     LWRK_CUNMLQ = INT( CDUMMY(1) )
                     IF ( ERREST ) THEN
                       OPTWRK = MAX( N+LWRK_CGEQP3, N+LWCON,  2*N+N**2+LWCON, 2*N+LWRK_CGEQRF, 2*N+LWRK_CGEQP3N, 2*N+N**2+N+LWRK_CGELQF, 2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWRK_CGESVJ, 2*N+N**2+N+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR, 2*N+N**2+N+LWRK_CUNMLQ, N+N**2+LWRK_CGESVJU, N+LWRK_CUNMQRM )
                     ELSE
                       OPTWRK = MAX( N+LWRK_CGEQP3,   2*N+N**2+LWCON, 2*N+LWRK_CGEQRF, 2*N+LWRK_CGEQP3N, 2*N+N**2+N+LWRK_CGELQF, 2*N+N**2+N+N**2+LWCON, 2*N+N**2+N+LWRK_CGESVJ, 2*N+N**2+N+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR, 2*N+N**2+N+LWRK_CUNMLQ, N+N**2+LWRK_CGESVJU, N+LWRK_CUNMQRM )
                     END IF
                 ELSE
                     CALL CGESVJ( 'L', 'U', 'V', N, N, U, LDU, SVA, N, V, LDV, CDUMMY, -1, RDUMMY, -1, IERR )
                     LWRK_CGESVJV = INT( CDUMMY(1) )
                     CALL CUNMQR( 'L', 'N', N, N, N, CDUMMY, N, CDUMMY, V, LDV, CDUMMY, -1, IERR )
                     LWRK_CUNMQR = INT( CDUMMY(1) )
                     CALL CUNMQR( 'L', 'N', M, N, N, A, LDA, CDUMMY, U, LDU, CDUMMY, -1, IERR )
                     LWRK_CUNMQRM = INT( CDUMMY(1) )
                     IF ( ERREST ) THEN
                        OPTWRK = MAX( N+LWRK_CGEQP3, N+LWCON,    2*N+LWRK_CGEQRF, 2*N+N**2, 2*N+N**2+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR,N+LWRK_CUNMQRM )
                     ELSE
                        OPTWRK = MAX( N+LWRK_CGEQP3, 2*N+LWRK_CGEQRF,   2*N+N**2, 2*N+N**2+LWRK_CGESVJV, 2*N+N**2+N+LWRK_CUNMQR, N+LWRK_CUNMQRM )
                     END IF
                 END IF
             END IF
             IF ( L2TRAN .OR. ROWPIV ) THEN
                 MINRWRK = MAX( 7, 2*M,  LRWQP3, LRWSVDJ, LRWCON )
             ELSE
                 MINRWRK = MAX( 7, LRWQP3, LRWSVDJ, LRWCON )
             END IF
          END IF
          MINWRK = MAX( 2, MINWRK )
          OPTWRK = MAX( OPTWRK, MINWRK )
          IF ( LWORK  .LT. MINWRK  .AND. (.NOT.LQUERY) ) INFO = - 17
          IF ( LRWORK .LT. MINRWRK .AND. (.NOT.LQUERY) ) INFO = - 19
      END IF
*
      IF ( INFO .NE. 0 ) THEN
*       #:(
         CALL XERBLA( 'CGEJSV', - INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
          CWORK(1) = OPTWRK
          CWORK(2) = MINWRK
          RWORK(1) = MINRWRK
          IWORK(1) = MAX( 4, MINIWRK )
          RETURN
      END IF
*
*     Quick return for void matrix (Y3K safe)
* #:)
      IF ( ( M .EQ. 0 ) .OR. ( N .EQ. 0 ) ) THEN
         IWORK(1:4) = 0
         RWORK(1:7) = 0
         RETURN
      ENDIF
*
*     Determine whether the matrix U should be M x N or M x M
*
      IF ( LSVEC ) THEN
         N1 = N
         IF ( LSAME( JOBU, 'F' ) ) N1 = M
      END IF
*
*     Set numerical parameters
*
*!    NOTE: Make sure SLAMCH() does not fail on the target architecture.
*
      EPSLN = SLAMCH('Epsilon')
      SFMIN = SLAMCH('SafeMinimum')
      SMALL = SFMIN / EPSLN
      BIG   = SLAMCH('O')
*     BIG   = ONE / SFMIN
*
*     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
*
*(!)  If necessary, scale SVA() to protect the largest norm from
*     overflow. It is possible that this scaling pushes the smallest
*     column norm left from the underflow threshold (extreme case).
*
      SCALEM  = ONE / SQRT(REAL(M)*REAL(N))
      NOSCAL  = .TRUE.
      GOSCAL  = .TRUE.
      DO 1874 p = 1, N
         AAPP = ZERO
         AAQQ = ONE
         CALL CLASSQ( M, A(1,p), 1, AAPP, AAQQ )
         IF ( AAPP .GT. BIG ) THEN
            INFO = - 9
            CALL XERBLA( 'CGEJSV', -INFO )
            RETURN
         END IF
         AAQQ = SQRT(AAQQ)
         IF ( ( AAPP .LT. (BIG / AAQQ) ) .AND. NOSCAL  ) THEN
            SVA(p)  = AAPP * AAQQ
         ELSE
            NOSCAL  = .FALSE.
            SVA(p)  = AAPP * ( AAQQ * SCALEM )
            IF ( GOSCAL ) THEN
               GOSCAL = .FALSE.
               CALL SSCAL( p-1, SCALEM, SVA, 1 )
            END IF
         END IF
 1874 CONTINUE
*
      IF ( NOSCAL ) SCALEM = ONE
*
      AAPP = ZERO
      AAQQ = BIG
      DO 4781 p = 1, N
         AAPP = MAX( AAPP, SVA(p) )
         IF ( SVA(p) .NE. ZERO ) AAQQ = MIN( AAQQ, SVA(p) )
 4781 CONTINUE
*
*     Quick return for zero M x N matrix
* #:)
      IF ( AAPP .EQ. ZERO ) THEN
         IF ( LSVEC ) CALL CLASET( 'G', M, N1, CZERO, CONE, U, LDU )
         IF ( RSVEC ) CALL CLASET( 'G', N, N,  CZERO, CONE, V, LDV )
         RWORK(1) = ONE
         RWORK(2) = ONE
         IF ( ERREST ) RWORK(3) = ONE
         IF ( LSVEC .AND. RSVEC ) THEN
            RWORK(4) = ONE
            RWORK(5) = ONE
         END IF
         IF ( L2TRAN ) THEN
            RWORK(6) = ZERO
            RWORK(7) = ZERO
         END IF
         IWORK(1) = 0
         IWORK(2) = 0
         IWORK(3) = 0
         IWORK(4) = -1
         RETURN
      END IF
*
*     Issue warning if denormalized column norms detected. Override the
*     high relative accuracy request. Issue licence to kill nonzero columns
*     (set them to zero) whose norm is less than sigma_max / BIG (roughly).
* #:(
      WARNING = 0
      IF ( AAQQ .LE. SFMIN ) THEN
         L2RANK = .TRUE.
         L2KILL = .TRUE.
         WARNING = 1
      END IF
*
*     Quick return for one-column matrix
* #:)
      IF ( N .EQ. 1 ) THEN
*
         IF ( LSVEC ) THEN
            CALL CLASCL( 'G',0,0,SVA(1),SCALEM, M,1,A(1,1),LDA,IERR )
            CALL CLACPY( 'A', M, 1, A, LDA, U, LDU )
*           computing all M left singular vectors of the M x 1 matrix
            IF ( N1 .NE. N  ) THEN
              CALL CGEQRF( M, N, U,LDU, CWORK, CWORK(N+1),LWORK-N,IERR )
              CALL CUNGQR( M,N1,1, U,LDU,CWORK,CWORK(N+1),LWORK-N,IERR )
              CALL CCOPY( M, A(1,1), 1, U(1,1), 1 )
            END IF
         END IF
         IF ( RSVEC ) THEN
             V(1,1) = CONE
         END IF
         IF ( SVA(1) .LT. (BIG*SCALEM) ) THEN
            SVA(1)  = SVA(1) / SCALEM
            SCALEM  = ONE
         END IF
         RWORK(1) = ONE / SCALEM
         RWORK(2) = ONE
         IF ( SVA(1) .NE. ZERO ) THEN
            IWORK(1) = 1
            IF ( ( SVA(1) / SCALEM) .GE. SFMIN ) THEN
               IWORK(2) = 1
            ELSE
               IWORK(2) = 0
            END IF
         ELSE
            IWORK(1) = 0
            IWORK(2) = 0
         END IF
         IWORK(3) = 0
         IWORK(4) = -1
         IF ( ERREST ) RWORK(3) = ONE
         IF ( LSVEC .AND. RSVEC ) THEN
            RWORK(4) = ONE
            RWORK(5) = ONE
         END IF
         IF ( L2TRAN ) THEN
            RWORK(6) = ZERO
            RWORK(7) = ZERO
         END IF
         RETURN
*
      END IF
*
      TRANSP = .FALSE.
*
      AATMAX = -ONE
      AATMIN =  BIG
      IF ( ROWPIV .OR. L2TRAN ) THEN
*
*     Compute the row norms, needed to determine row pivoting sequence
*     (in the case of heavily row weighted A, row pivoting is strongly
*     advised) and to collect information needed to compare the
*     structures of A * A^* and A^* * A (in the case L2TRAN.EQ..TRUE.).
*
         IF ( L2TRAN ) THEN
            DO 1950 p = 1, M
               XSC   = ZERO
               TEMP1 = ONE
               CALL CLASSQ( N, A(p,1), LDA, XSC, TEMP1 )
*              CLASSQ gets both the ell_2 and the ell_infinity norm
*              in one pass through the vector
               RWORK(M+p)  = XSC * SCALEM
               RWORK(p)    = XSC * (SCALEM*SQRT(TEMP1))
               AATMAX = MAX( AATMAX, RWORK(p) )
               IF (RWORK(p) .NE. ZERO)  AATMIN = MIN(AATMIN,RWORK(p))
 1950       CONTINUE
         ELSE
            DO 1904 p = 1, M
               RWORK(M+p) = SCALEM*ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
               AATMAX = MAX( AATMAX, RWORK(M+p) )
               AATMIN = MIN( AATMIN, RWORK(M+p) )
 1904       CONTINUE
         END IF
*
      END IF
*
*     For square matrix A try to determine whether A^*  would be better
*     input for the preconditioned Jacobi SVD, with faster convergence.
*     The decision is based on an O(N) function of the vector of column
*     and row norms of A, based on the Shannon entropy. This should give
*     the right choice in most cases when the difference actually matters.
*     It may fail and pick the slower converging side.
*
      ENTRA  = ZERO
      ENTRAT = ZERO
      IF ( L2TRAN ) THEN
*
         XSC   = ZERO
         TEMP1 = ONE
         CALL SLASSQ( N, SVA, 1, XSC, TEMP1 )
         TEMP1 = ONE / TEMP1
*
         ENTRA = ZERO
         DO 1113 p = 1, N
            BIG1  = ( ( SVA(p) / XSC )**2 ) * TEMP1
            IF ( BIG1 .NE. ZERO ) ENTRA = ENTRA + BIG1 * ALOG(BIG1)
 1113    CONTINUE
         ENTRA = - ENTRA / ALOG(REAL(N))
*
*        Now, SVA().^2/Trace(A^* * A) is a point in the probability simplex.
*        It is derived from the diagonal of  A^* * A.  Do the same with the
*        diagonal of A * A^*, compute the entropy of the corresponding
*        probability distribution. Note that A * A^* and A^* * A have the
*        same trace.
*
         ENTRAT = ZERO
         DO 1114 p = 1, M
            BIG1 = ( ( RWORK(p) / XSC )**2 ) * TEMP1
            IF ( BIG1 .NE. ZERO ) ENTRAT = ENTRAT + BIG1 * ALOG(BIG1)
 1114    CONTINUE
         ENTRAT = - ENTRAT / ALOG(REAL(M))
*
*        Analyze the entropies and decide A or A^*. Smaller entropy
*        usually means better input for the algorithm.
*
         TRANSP = ( ENTRAT .LT. ENTRA )
*
*        If A^* is better than A, take the adjoint of A. This is allowed
*        only for square matrices, M=N.
         IF ( TRANSP ) THEN
*           In an optimal implementation, this trivial transpose
*           should be replaced with faster transpose.
            DO 1115 p = 1, N - 1
               A(p,p) = CONJG(A(p,p))
               DO 1116 q = p + 1, N
                   CTEMP = CONJG(A(q,p))
                  A(q,p) = CONJG(A(p,q))
                  A(p,q) = CTEMP
 1116          CONTINUE
 1115       CONTINUE
            A(N,N) = CONJG(A(N,N))
            DO 1117 p = 1, N
               RWORK(M+p) = SVA(p)
               SVA(p) = RWORK(p)
*              previously computed row 2-norms are now column 2-norms
*              of the transposed matrix
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
*
            ROWPIV = .TRUE.
         END IF
*
      END IF
*     END IF L2TRAN
*
*     Scale the matrix so that its maximal singular value remains less
*     than SQRT(BIG) -- the matrix is scaled so that its maximal column
*     has Euclidean norm equal to SQRT(BIG/N). The only reason to keep
*     SQRT(BIG) instead of BIG is the fact that CGEJSV uses LAPACK and
*     BLAS routines that, in some implementations, are not capable of
*     working in the full interval [SFMIN,BIG] and that they may provoke
*     overflows in the intermediate results. If the singular values spread
*     from SFMIN to BIG, then CGESVJ will compute them. So, in that case,
*     one should use CGESVJ instead of CGEJSV.
      BIG1   = SQRT( BIG )
      TEMP1  = SQRT( BIG / REAL(N) )
*     >> for future updates: allow bigger range, i.e. the largest column
*     will be allowed up to BIG/N and CGESVJ will do the rest. However, for
*     this all other (LAPACK) components must allow such a range.
*     TEMP1  = BIG/REAL(N)
*     TEMP1  = BIG * EPSLN  this should 'almost' work with current LAPACK components
      CALL SLASCL( 'G', 0, 0, AAPP, TEMP1, N, 1, SVA, N, IERR )
      IF ( AAQQ .GT. (AAPP * SFMIN) ) THEN
          AAQQ = ( AAQQ / AAPP ) * TEMP1
      ELSE
          AAQQ = ( AAQQ * TEMP1 ) / AAPP
      END IF
      TEMP1 = TEMP1 * SCALEM
      CALL CLASCL( 'G', 0, 0, AAPP, TEMP1, M, N, A, LDA, IERR )
*
*     To undo scaling at the end of this procedure, multiply the
*     computed singular values with USCAL2 / USCAL1.
*
      USCAL1 = TEMP1
      USCAL2 = AAPP
*
      IF ( L2KILL ) THEN
*        L2KILL enforces computation of nonzero singular values in
*        the restricted range of condition number of the initial A,
*        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN).
         XSC = SQRT( SFMIN )
      ELSE
         XSC = SMALL
*
*        Now, if the condition number of A is too big,
*        sigma_max(A) / sigma_min(A) .GT. SQRT(BIG/N) * EPSLN / SFMIN,
*        as a precaution measure, the full SVD is computed using CGESVJ
*        with accumulated Jacobi rotations. This provides numerically
*        more robust computation, at the cost of slightly increased run
*        time. Depending on the concrete implementation of BLAS and LAPACK
*        (i.e. how they behave in presence of extreme ill-conditioning) the
*        implementor may decide to remove this switch.
         IF ( ( AAQQ.LT.SQRT(SFMIN) ) .AND. LSVEC .AND. RSVEC ) THEN
            JRACC = .TRUE.
         END IF
*
      END IF
      IF ( AAQQ .LT. XSC ) THEN
         DO 700 p = 1, N
            IF ( SVA(p) .LT. XSC ) THEN
               CALL CLASET( 'A', M, 1, CZERO, CZERO, A(1,p), LDA )
               SVA(p) = ZERO
            END IF
 700     CONTINUE
      END IF
*
*     Preconditioning using QR factorization with pivoting
*
      IF ( ROWPIV ) THEN
*        Optional row permutation (Bjoerck row pivoting):
*        A result by Cox and Higham shows that the Bjoerck's
*        row pivoting combined with standard column pivoting
*        has similar effect as Powell-Reid complete pivoting.
*        The ell-infinity norms of A are made nonincreasing.
         IF ( ( LSVEC .AND. RSVEC ) .AND. .NOT.( JRACC ) ) THEN
              IWOFF = 2*N
         ELSE
              IWOFF = N
         END IF
         DO 1952 p = 1, M - 1
            q = ISAMAX( M-p+1, RWORK(M+p), 1 ) + p - 1
            IWORK(IWOFF+p) = q
            IF ( p .NE. q ) THEN
               TEMP1      = RWORK(M+p)
               RWORK(M+p) = RWORK(M+q)
               RWORK(M+q) = TEMP1
            END IF
 1952    CONTINUE
         CALL CLASWP( N, A, LDA, 1, M-1, IWORK(IWOFF+1), 1 )
      END IF
*
*     End of the preparation phase (scaling, optional sorting and
*     transposing, optional flushing of small columns).
*
*     Preconditioning
*
*     If the full SVD is needed, the right singular vectors are computed
*     from a matrix equation, and for that we need theoretical analysis
*     of the Businger-Golub pivoting. So we use CGEQP3 as the first RR QRF.
*     In all other cases the first RR QRF can be chosen by other criteria
*     (eg speed by replacing global with restricted window pivoting, such
*     as in xGEQPX from TOMS # 782). Good results will be obtained using
*     xGEQPX with properly (!) chosen numerical parameters.
*     Any improvement of CGEQP3 improves overall performance of CGEJSV.
*
*     A * P1 = Q1 * [ R1^* 0]^*:
      DO 1963 p = 1, N
*        .. all columns are free columns
         IWORK(p) = 0
 1963 CONTINUE
      CALL CGEQP3( M, N, A, LDA, IWORK, CWORK, CWORK(N+1), LWORK-N, RWORK, IERR )
*
*     The upper triangular matrix R1 from the first QRF is inspected for
*     rank deficiency and possibilities for deflation, or possible
*     ill-conditioning. Depending on the user specified flag L2RANK,
*     the procedure explores possibilities to reduce the numerical
*     rank by inspecting the computed upper triangular factor. If
*     L2RANK or L2ABER are up, then CGEJSV will compute the SVD of
*     A + dA, where ||dA|| <= f(M,N)*EPSLN.
*
      NR = 1
      IF ( L2ABER ) THEN
*        Standard absolute error bound suffices. All sigma_i with
*        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
*        aggressive enforcement of lower numerical rank by introducing a
*        backward error of the order of N*EPSLN*||A||.
         TEMP1 = SQRT(REAL(N))*EPSLN
         DO 3001 p = 2, N
            IF ( ABS(A(p,p)) .GE. (TEMP1*ABS(A(1,1))) ) THEN
               NR = NR + 1
            ELSE
               GO TO 3002
            END IF
 3001    CONTINUE
 3002    CONTINUE
      ELSE IF ( L2RANK ) THEN
*        .. similarly as above, only slightly more gentle (less aggressive).
*        Sudden drop on the diagonal of R1 is used as the criterion for
*        close-to-rank-deficient.
         TEMP1 = SQRT(SFMIN)
         DO 3401 p = 2, N
            IF ( ( ABS(A(p,p)) .LT. (EPSLN*ABS(A(p-1,p-1))) ) .OR. ( ABS(A(p,p)) .LT. SMALL ) .OR. ( L2KILL .AND. (ABS(A(p,p)) .LT. TEMP1) ) ) GO TO 3402
            NR = NR + 1
 3401    CONTINUE
 3402    CONTINUE
*
      ELSE
*        The goal is high relative accuracy. However, if the matrix
*        has high scaled condition number the relative accuracy is in
*        general not feasible. Later on, a condition number estimator
*        will be deployed to estimate the scaled condition number.
*        Here we just remove the underflowed part of the triangular
*        factor. This prevents the situation in which the code is
*        working hard to get the accuracy not warranted by the data.
         TEMP1  = SQRT(SFMIN)
         DO 3301 p = 2, N
            IF ( ( ABS(A(p,p)) .LT. SMALL ) .OR. ( L2KILL .AND. (ABS(A(p,p)) .LT. TEMP1) ) ) GO TO 3302
            NR = NR + 1
 3301    CONTINUE
 3302    CONTINUE
*
      END IF
*
      ALMORT = .FALSE.
      IF ( NR .EQ. N ) THEN
         MAXPRJ = ONE
         DO 3051 p = 2, N
            TEMP1  = ABS(A(p,p)) / SVA(IWORK(p))
            MAXPRJ = MIN( MAXPRJ, TEMP1 )
 3051    CONTINUE
         IF ( MAXPRJ**2 .GE. ONE - REAL(N)*EPSLN ) ALMORT = .TRUE.
      END IF
*
*
      SCONDA = - ONE
      CONDR1 = - ONE
      CONDR2 = - ONE
*
      IF ( ERREST ) THEN
         IF ( N .EQ. NR ) THEN
            IF ( RSVEC ) THEN
*              .. V is available as workspace
               CALL CLACPY( 'U', N, N, A, LDA, V, LDV )
               DO 3053 p = 1, N
                  TEMP1 = SVA(IWORK(p))
                  CALL CSSCAL( p, ONE/TEMP1, V(1,p), 1 )
 3053          CONTINUE
               IF ( LSVEC )THEN
                   CALL CPOCON( 'U', N, V, LDV, ONE, TEMP1, CWORK(N+1), RWORK, IERR )
               ELSE
                   CALL CPOCON( 'U', N, V, LDV, ONE, TEMP1, CWORK, RWORK, IERR )
               END IF
*
            ELSE IF ( LSVEC ) THEN
*              .. U is available as workspace
               CALL CLACPY( 'U', N, N, A, LDA, U, LDU )
               DO 3054 p = 1, N
                  TEMP1 = SVA(IWORK(p))
                  CALL CSSCAL( p, ONE/TEMP1, U(1,p), 1 )
 3054          CONTINUE
               CALL CPOCON( 'U', N, U, LDU, ONE, TEMP1, CWORK(N+1), RWORK, IERR )
            ELSE
               CALL CLACPY( 'U', N, N, A, LDA, CWORK, N )
*[]            CALL CLACPY( 'U', N, N, A, LDA, CWORK(N+1), N )
*              Change: here index shifted by N to the left, CWORK(1:N)
*              not needed for SIGMA only computation
               DO 3052 p = 1, N
                  TEMP1 = SVA(IWORK(p))
*[]               CALL CSSCAL( p, ONE/TEMP1, CWORK(N+(p-1)*N+1), 1 )
                  CALL CSSCAL( p, ONE/TEMP1, CWORK((p-1)*N+1), 1 )
 3052          CONTINUE
*           .. the columns of R are scaled to have unit Euclidean lengths.
*[]               CALL CPOCON( 'U', N, CWORK(N+1), N, ONE, TEMP1,
*[]     $              CWORK(N+N*N+1), RWORK, IERR )
               CALL CPOCON( 'U', N, CWORK, N, ONE, TEMP1, CWORK(N*N+1), RWORK, IERR )
*
            END IF
            IF ( TEMP1 .NE. ZERO ) THEN
               SCONDA = ONE / SQRT(TEMP1)
            ELSE
               SCONDA = - ONE
            END IF
*           SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1).
*           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
         ELSE
            SCONDA = - ONE
         END IF
      END IF
*
      L2PERT = L2PERT .AND. ( ABS( A(1,1)/A(NR,NR) ) .GT. SQRT(BIG1) )
*     If there is no violent scaling, artificial perturbation is not needed.
*
*     Phase 3:
*
      IF ( .NOT. ( RSVEC .OR. LSVEC ) ) THEN
*
*         Singular Values only
*
*         .. transpose A(1:NR,1:N)
         DO 1946 p = 1, MIN( N-1, NR )
            CALL CCOPY( N-p, A(p,p+1), LDA, A(p+1,p), 1 )
            CALL CLACGV( N-p+1, A(p,p), 1 )
 1946    CONTINUE
         IF ( NR .EQ. N ) A(N,N) = CONJG(A(N,N))
*
*        The following two DO-loops introduce small relative perturbation
*        into the strict upper triangle of the lower triangular matrix.
*        Small entries below the main diagonal are also changed.
*        This modification is useful if the computing environment does not
*        provide/allow FLUSH TO ZERO underflow, for it prevents many
*        annoying denormalized numbers in case of strongly scaled matrices.
*        The perturbation is structured so that it does not introduce any
*        new perturbation of the singular values, and it does not destroy
*        the job done by the preconditioner.
*        The licence for this perturbation is in the variable L2PERT, which
*        should be .FALSE. if FLUSH TO ZERO underflow is active.
*
         IF ( .NOT. ALMORT ) THEN
*
            IF ( L2PERT ) THEN
*              XSC = SQRT(SMALL)
               XSC = EPSLN / REAL(N)
               DO 4947 q = 1, NR
                  CTEMP = CMPLX(XSC*ABS(A(q,q)),ZERO)
                  DO 4949 p = 1, N
                     IF ( ( (p.GT.q) .AND. (ABS(A(p,q)).LE.TEMP1) ) .OR. ( p .LT. q ) )
*     $                     A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) )
     $                     A(p,q) = CTEMP
 4949             CONTINUE
 4947          CONTINUE
            ELSE
               CALL CLASET( 'U', NR-1,NR-1, CZERO,CZERO, A(1,2),LDA )
            END IF
*
*            .. second preconditioning using the QR factorization
*
            CALL CGEQRF( N,NR, A,LDA, CWORK, CWORK(N+1),LWORK-N, IERR )
*
*           .. and transpose upper to lower triangular
            DO 1948 p = 1, NR - 1
               CALL CCOPY( NR-p, A(p,p+1), LDA, A(p+1,p), 1 )
               CALL CLACGV( NR-p+1, A(p,p), 1 )
 1948       CONTINUE
*
         END IF
*
*           Row-cyclic Jacobi SVD algorithm with column pivoting
*
*           .. again some perturbation (a "background noise") is added
*           to drown denormals
            IF ( L2PERT ) THEN
*              XSC = SQRT(SMALL)
               XSC = EPSLN / REAL(N)
               DO 1947 q = 1, NR
                  CTEMP = CMPLX(XSC*ABS(A(q,q)),ZERO)
                  DO 1949 p = 1, NR
                     IF ( ( (p.GT.q) .AND. (ABS(A(p,q)).LE.TEMP1) ) .OR. ( p .LT. q ) )
*     $                   A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) )
     $                   A(p,q) = CTEMP
 1949             CONTINUE
 1947          CONTINUE
            ELSE
               CALL CLASET( 'U', NR-1, NR-1, CZERO, CZERO, A(1,2), LDA )
            END IF
*
*           .. and one-sided Jacobi rotations are started on a lower
*           triangular matrix (plus perturbation which is ignored in
*           the part which destroys triangular form (confusing?!))
*
            CALL CGESVJ( 'L', 'N', 'N', NR, NR, A, LDA, SVA, N, V, LDV, CWORK, LWORK, RWORK, LRWORK, INFO )
*
            SCALEM  = RWORK(1)
            NUMRANK = NINT(RWORK(2))
*
*
      ELSE IF ( ( RSVEC .AND. ( .NOT. LSVEC ) .AND. ( .NOT. JRACC ) )  .OR. ( JRACC .AND. ( .NOT. LSVEC ) .AND. ( NR .NE. N ) ) ) THEN
*
*        -> Singular Values and Right Singular Vectors <-
*
         IF ( ALMORT ) THEN
*
*           .. in this case NR equals N
            DO 1998 p = 1, NR
               CALL CCOPY( N-p+1, A(p,p), LDA, V(p,p), 1 )
               CALL CLACGV( N-p+1, V(p,p), 1 )
 1998       CONTINUE
            CALL CLASET( 'U', NR-1,NR-1, CZERO, CZERO, V(1,2), LDV )
*
            CALL CGESVJ( 'L','U','N', N, NR, V, LDV, SVA, NR, A, LDA, CWORK, LWORK, RWORK, LRWORK, INFO )
            SCALEM  = RWORK(1)
            NUMRANK = NINT(RWORK(2))

         ELSE
*
*        .. two more QR factorizations ( one QRF is not enough, two require
*        accumulated product of Jacobi rotations, three are perfect )
*
            CALL CLASET( 'L', NR-1,NR-1, CZERO, CZERO, A(2,1), LDA )
            CALL CGELQF( NR,N, A, LDA, CWORK, CWORK(N+1), LWORK-N, IERR)
            CALL CLACPY( 'L', NR, NR, A, LDA, V, LDV )
            CALL CLASET( 'U', NR-1,NR-1, CZERO, CZERO, V(1,2), LDV )
            CALL CGEQRF( NR, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR )
            DO 8998 p = 1, NR
               CALL CCOPY( NR-p+1, V(p,p), LDV, V(p,p), 1 )
               CALL CLACGV( NR-p+1, V(p,p), 1 )
 8998       CONTINUE
            CALL CLASET('U', NR-1, NR-1, CZERO, CZERO, V(1,2), LDV)
*
            CALL CGESVJ( 'L', 'U','N', NR, NR, V,LDV, SVA, NR, U, LDU, CWORK(N+1), LWORK-N, RWORK, LRWORK, INFO )
            SCALEM  = RWORK(1)
            NUMRANK = NINT(RWORK(2))
            IF ( NR .LT. N ) THEN
               CALL CLASET( 'A',N-NR, NR, CZERO,CZERO, V(NR+1,1),  LDV )
               CALL CLASET( 'A',NR, N-NR, CZERO,CZERO, V(1,NR+1),  LDV )
               CALL CLASET( 'A',N-NR,N-NR,CZERO,CONE, V(NR+1,NR+1),LDV )
            END IF
*
         CALL CUNMLQ( 'L', 'C', N, N, NR, A, LDA, CWORK, V, LDV, CWORK(N+1), LWORK-N, IERR )
*
         END IF
*         .. permute the rows of V
*         DO 8991 p = 1, N
*            CALL CCOPY( N, V(p,1), LDV, A(IWORK(p),1), LDA )
* 8991    CONTINUE
*         CALL CLACPY( 'All', N, N, A, LDA, V, LDV )
         CALL CLAPMR( .FALSE., N, N, V, LDV, IWORK )
*
          IF ( TRANSP ) THEN
            CALL CLACPY( 'A', N, N, V, LDV, U, LDU )
          END IF
*
      ELSE IF ( JRACC .AND. (.NOT. LSVEC) .AND. ( NR.EQ. N ) ) THEN
*
         CALL CLASET( 'L', N-1,N-1, CZERO, CZERO, A(2,1), LDA )
*
         CALL CGESVJ( 'U','N','V', N, N, A, LDA, SVA, N, V, LDV, CWORK, LWORK, RWORK, LRWORK, INFO )
          SCALEM  = RWORK(1)
          NUMRANK = NINT(RWORK(2))
          CALL CLAPMR( .FALSE., N, N, V, LDV, IWORK )
*
      ELSE IF ( LSVEC .AND. ( .NOT. RSVEC ) ) THEN
*
*        .. Singular Values and Left Singular Vectors                 ..
*
*        .. second preconditioning step to avoid need to accumulate
*        Jacobi rotations in the Jacobi iterations.
         DO 1965 p = 1, NR
            CALL CCOPY( N-p+1, A(p,p), LDA, U(p,p), 1 )
            CALL CLACGV( N-p+1, U(p,p), 1 )
 1965    CONTINUE
         CALL CLASET( 'U', NR-1, NR-1, CZERO, CZERO, U(1,2), LDU )
*
         CALL CGEQRF( N, NR, U, LDU, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR )
*
         DO 1967 p = 1, NR - 1
            CALL CCOPY( NR-p, U(p,p+1), LDU, U(p+1,p), 1 )
            CALL CLACGV( N-p+1, U(p,p), 1 )
 1967    CONTINUE
         CALL CLASET( 'U', NR-1, NR-1, CZERO, CZERO, U(1,2), LDU )
*
         CALL CGESVJ( 'L', 'U', 'N', NR,NR, U, LDU, SVA, NR, A, LDA, CWORK(N+1), LWORK-N, RWORK, LRWORK, INFO )
         SCALEM  = RWORK(1)
         NUMRANK = NINT(RWORK(2))
*
         IF ( NR .LT. M ) THEN
            CALL CLASET( 'A',  M-NR, NR,CZERO, CZERO, U(NR+1,1), LDU )
            IF ( NR .LT. N1 ) THEN
               CALL CLASET( 'A',NR, N1-NR, CZERO, CZERO, U(1,NR+1),LDU )
               CALL CLASET( 'A',M-NR,N1-NR,CZERO,CONE,U(NR+1,NR+1),LDU )
            END IF
         END IF
*
         CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR )
*
         IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
*
         DO 1974 p = 1, N1
            XSC = ONE / SCNRM2( M, U(1,p), 1 )
            CALL CSSCAL( M, XSC, U(1,p), 1 )
 1974    CONTINUE
*
         IF ( TRANSP ) THEN
            CALL CLACPY( 'A', N, N, U, LDU, V, LDV )
         END IF
*
      ELSE
*
*        .. Full SVD ..
*
         IF ( .NOT. JRACC ) THEN
*
         IF ( .NOT. ALMORT ) THEN
*
*           Second Preconditioning Step (QRF [with pivoting])
*           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
*           equivalent to an LQF CALL. Since in many libraries the QRF
*           seems to be better optimized than the LQF, we do explicit
*           transpose and use the QRF. This is subject to changes in an
*           optimized implementation of CGEJSV.
*
            DO 1968 p = 1, NR
               CALL CCOPY( N-p+1, A(p,p), LDA, V(p,p), 1 )
               CALL CLACGV( N-p+1, V(p,p), 1 )
 1968       CONTINUE
*
*           .. the following two loops perturb small entries to avoid
*           denormals in the second QR factorization, where they are
*           as good as zeros. This is done to avoid painfully slow
*           computation with denormals. The relative size of the perturbation
*           is a parameter that can be changed by the implementer.
*           This perturbation device will be obsolete on machines with
*           properly implemented arithmetic.
*           To switch it off, set L2PERT=.FALSE. To remove it from  the
*           code, remove the action under L2PERT=.TRUE., leave the ELSE part.
*           The following two loops should be blocked and fused with the
*           transposed copy above.
*
            IF ( L2PERT ) THEN
               XSC = SQRT(SMALL)
               DO 2969 q = 1, NR
                  CTEMP = CMPLX(XSC*ABS( V(q,q) ),ZERO)
                  DO 2968 p = 1, N
                     IF ( ( p .GT. q ) .AND. ( ABS(V(p,q)) .LE. TEMP1 ) .OR. ( p .LT. q ) )
*     $                   V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) )
     $                   V(p,q) = CTEMP
                     IF ( p .LT. q ) V(p,q) = - V(p,q)
 2968             CONTINUE
 2969          CONTINUE
            ELSE
               CALL CLASET( 'U', NR-1, NR-1, CZERO, CZERO, V(1,2), LDV )
            END IF
*
*           Estimate the row scaled condition number of R1
*           (If R1 is rectangular, N > NR, then the condition number
*           of the leading NR x NR submatrix is estimated.)
*
            CALL CLACPY( 'L', NR, NR, V, LDV, CWORK(2*N+1), NR )
            DO 3950 p = 1, NR
               TEMP1 = SCNRM2(NR-p+1,CWORK(2*N+(p-1)*NR+p),1)
               CALL CSSCAL(NR-p+1,ONE/TEMP1,CWORK(2*N+(p-1)*NR+p),1)
 3950       CONTINUE
            CALL CPOCON('L',NR,CWORK(2*N+1),NR,ONE,TEMP1, CWORK(2*N+NR*NR+1),RWORK,IERR)
            CONDR1 = ONE / SQRT(TEMP1)
*           .. here need a second opinion on the condition number
*           .. then assume worst case scenario
*           R1 is OK for inverse <=> CONDR1 .LT. REAL(N)
*           more conservative    <=> CONDR1 .LT. SQRT(REAL(N))
*
            COND_OK = SQRT(SQRT(REAL(NR)))
*[TP]       COND_OK is a tuning parameter.
*
            IF ( CONDR1 .LT. COND_OK ) THEN
*              .. the second QRF without pivoting. Note: in an optimized
*              implementation, this QRF should be implemented as the QRF
*              of a lower triangular matrix.
*              R1^* = Q2 * R2
               CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR )
*
               IF ( L2PERT ) THEN
                  XSC = SQRT(SMALL)/EPSLN
                  DO 3959 p = 2, NR
                     DO 3958 q = 1, p - 1
                        CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), ZERO)
                        IF ( ABS(V(q,p)) .LE. TEMP1 )
*     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) )
     $                     V(q,p) = CTEMP
 3958                CONTINUE
 3959             CONTINUE
               END IF
*
               IF ( NR .NE. N ) CALL CLACPY( 'A', N, NR, V, LDV, CWORK(2*N+1), N )
*              .. save ...
*
*           .. this transposed copy should be better than naive
               DO 1969 p = 1, NR - 1
                  CALL CCOPY( NR-p, V(p,p+1), LDV, V(p+1,p), 1 )
                  CALL CLACGV(NR-p+1, V(p,p), 1 )
 1969          CONTINUE
               V(NR,NR)=CONJG(V(NR,NR))
*
               CONDR2 = CONDR1
*
            ELSE
*
*              .. ill-conditioned case: second QRF with pivoting
*              Note that windowed pivoting would be equally good
*              numerically, and more run-time efficient. So, in
*              an optimal implementation, the next call to CGEQP3
*              should be replaced with eg. CALL CGEQPX (ACM TOMS #782)
*              with properly (carefully) chosen parameters.
*
*              R1^* * P2 = Q2 * R2
               DO 3003 p = 1, NR
                  IWORK(N+p) = 0
 3003          CONTINUE
               CALL CGEQP3( N, NR, V, LDV, IWORK(N+1), CWORK(N+1), CWORK(2*N+1), LWORK-2*N, RWORK, IERR )
**               CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1),
**     $              LWORK-2*N, IERR )
               IF ( L2PERT ) THEN
                  XSC = SQRT(SMALL)
                  DO 3969 p = 2, NR
                     DO 3968 q = 1, p - 1
                        CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), ZERO)
                        IF ( ABS(V(q,p)) .LE. TEMP1 )
*     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) )
     $                     V(q,p) = CTEMP
 3968                CONTINUE
 3969             CONTINUE
               END IF
*
               CALL CLACPY( 'A', N, NR, V, LDV, CWORK(2*N+1), N )
*
               IF ( L2PERT ) THEN
                  XSC = SQRT(SMALL)
                  DO 8970 p = 2, NR
                     DO 8971 q = 1, p - 1
                        CTEMP=CMPLX(XSC*MIN(ABS(V(p,p)),ABS(V(q,q))), ZERO)
*                        V(p,q) = - TEMP1*( V(q,p) / ABS(V(q,p)) )
                        V(p,q) = - CTEMP
 8971                CONTINUE
 8970             CONTINUE
               ELSE
                  CALL CLASET( 'L',NR-1,NR-1,CZERO,CZERO,V(2,1),LDV )
               END IF
*              Now, compute R2 = L3 * Q3, the LQ factorization.
               CALL CGELQF( NR, NR, V, LDV, CWORK(2*N+N*NR+1), CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR )
*              .. and estimate the condition number
               CALL CLACPY( 'L',NR,NR,V,LDV,CWORK(2*N+N*NR+NR+1),NR )
               DO 4950 p = 1, NR
                  TEMP1 = SCNRM2( p, CWORK(2*N+N*NR+NR+p), NR )
                  CALL CSSCAL( p, ONE/TEMP1, CWORK(2*N+N*NR+NR+p), NR )
 4950          CONTINUE
               CALL CPOCON( 'L',NR,CWORK(2*N+N*NR+NR+1),NR,ONE,TEMP1, CWORK(2*N+N*NR+NR+NR*NR+1),RWORK,IERR )
               CONDR2 = ONE / SQRT(TEMP1)
*
*
               IF ( CONDR2 .GE. COND_OK ) THEN
*                 .. save the Householder vectors used for Q3
*                 (this overwrites the copy of R2, as it will not be
*                 needed in this branch, but it does not overwrite the
*                 Huseholder vectors of Q2.).
                  CALL CLACPY( 'U', NR, NR, V, LDV, CWORK(2*N+1), N )
*                 .. and the rest of the information on Q3 is in
*                 WORK(2*N+N*NR+1:2*N+N*NR+N)
               END IF
*
            END IF
*
            IF ( L2PERT ) THEN
               XSC = SQRT(SMALL)
               DO 4968 q = 2, NR
                  CTEMP = XSC * V(q,q)
                  DO 4969 p = 1, q - 1
*                     V(p,q) = - TEMP1*( V(p,q) / ABS(V(p,q)) )
                     V(p,q) = - CTEMP
 4969             CONTINUE
 4968          CONTINUE
            ELSE
               CALL CLASET( 'U', NR-1,NR-1, CZERO,CZERO, V(1,2), LDV )
            END IF
*
*        Second preconditioning finished; continue with Jacobi SVD
*        The input matrix is lower triangular.
*
*        Recover the right singular vectors as solution of a well
*        conditioned triangular matrix equation.
*
            IF ( CONDR1 .LT. COND_OK ) THEN
*
               CALL CGESVJ( 'L','U','N',NR,NR,V,LDV,SVA,NR,U, LDU, CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,RWORK, LRWORK, INFO )
               SCALEM  = RWORK(1)
               NUMRANK = NINT(RWORK(2))
               DO 3970 p = 1, NR
                  CALL CCOPY(  NR, V(1,p), 1, U(1,p), 1 )
                  CALL CSSCAL( NR, SVA(p),    V(1,p), 1 )
 3970          CONTINUE

*        .. pick the right matrix equation and solve it
*
               IF ( NR .EQ. N ) THEN
* :))             .. best case, R1 is inverted. The solution of this matrix
*                 equation is Q2*V2 = the product of the Jacobi rotations
*                 used in CGESVJ, premultiplied with the orthogonal matrix
*                 from the second QR factorization.
                  CALL CTRSM('L','U','N','N', NR,NR,CONE, A,LDA, V,LDV)
               ELSE
*                 .. R1 is well conditioned, but non-square. Adjoint of R2
*                 is inverted to get the product of the Jacobi rotations
*                 used in CGESVJ. The Q-factor from the second QR
*                 factorization is then built in explicitly.
                  CALL CTRSM('L','U','C','N',NR,NR,CONE,CWORK(2*N+1), N,V,LDV)
                  IF ( NR .LT. N ) THEN
                  CALL CLASET('A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV)
                  CALL CLASET('A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV)
                  CALL CLASET('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV)
                  END IF
                  CALL CUNMQR('L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR)
               END IF
*
            ELSE IF ( CONDR2 .LT. COND_OK ) THEN
*
*              The matrix R2 is inverted. The solution of the matrix equation
*              is Q3^* * V3 = the product of the Jacobi rotations (applied to
*              the lower triangular L3 from the LQ factorization of
*              R2=L3*Q3), pre-multiplied with the transposed Q3.
               CALL CGESVJ( 'L', 'U', 'N', NR, NR, V, LDV, SVA, NR, U, LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, RWORK, LRWORK, INFO )
               SCALEM  = RWORK(1)
               NUMRANK = NINT(RWORK(2))
               DO 3870 p = 1, NR
                  CALL CCOPY( NR, V(1,p), 1, U(1,p), 1 )
                  CALL CSSCAL( NR, SVA(p),    U(1,p), 1 )
 3870          CONTINUE
               CALL CTRSM('L','U','N','N',NR,NR,CONE,CWORK(2*N+1),N, U,LDU)
*              .. apply the permutation from the second QR factorization
               DO 873 q = 1, NR
                  DO 872 p = 1, NR
                     CWORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q)
 872              CONTINUE
                  DO 874 p = 1, NR
                     U(p,q) = CWORK(2*N+N*NR+NR+p)
 874              CONTINUE
 873           CONTINUE
               IF ( NR .LT. N ) THEN
                  CALL CLASET( 'A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV )
                  CALL CLASET( 'A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV )
                  CALL CLASET('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV)
               END IF
               CALL CUNMQR( 'L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
            ELSE
*              Last line of defense.
* #:(          This is a rather pathological case: no scaled condition
*              improvement after two pivoted QR factorizations. Other
*              possibility is that the rank revealing QR factorization
*              or the condition estimator has failed, or the COND_OK
*              is set very close to ONE (which is unnecessary). Normally,
*              this branch should never be executed, but in rare cases of
*              failure of the RRQR or condition estimator, the last line of
*              defense ensures that CGEJSV completes the task.
*              Compute the full SVD of L3 using CGESVJ with explicit
*              accumulation of Jacobi rotations.
               CALL CGESVJ( 'L', 'U', 'V', NR, NR, V, LDV, SVA, NR, U, LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, RWORK, LRWORK, INFO )
               SCALEM  = RWORK(1)
               NUMRANK = NINT(RWORK(2))
               IF ( NR .LT. N ) THEN
                  CALL CLASET( 'A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV )
                  CALL CLASET( 'A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV )
                  CALL CLASET('A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV)
               END IF
               CALL CUNMQR( 'L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
*
               CALL CUNMLQ( 'L', 'C', NR, NR, NR, CWORK(2*N+1), N, CWORK(2*N+N*NR+1), U, LDU, CWORK(2*N+N*NR+NR+1), LWORK-2*N-N*NR-NR, IERR )
               DO 773 q = 1, NR
                  DO 772 p = 1, NR
                     CWORK(2*N+N*NR+NR+IWORK(N+p)) = U(p,q)
 772              CONTINUE
                  DO 774 p = 1, NR
                     U(p,q) = CWORK(2*N+N*NR+NR+p)
 774              CONTINUE
 773           CONTINUE
*
            END IF
*
*           Permute the rows of V using the (column) permutation from the
*           first QRF. Also, scale the columns to make them unit in
*           Euclidean norm. This applies to all cases.
*
            TEMP1 = SQRT(REAL(N)) * EPSLN
            DO 1972 q = 1, N
               DO 972 p = 1, N
                  CWORK(2*N+N*NR+NR+IWORK(p)) = V(p,q)
  972          CONTINUE
               DO 973 p = 1, N
                  V(p,q) = CWORK(2*N+N*NR+NR+p)
  973          CONTINUE
               XSC = ONE / SCNRM2( N, V(1,q), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL CSSCAL( N, XSC, V(1,q), 1 )
 1972       CONTINUE
*           At this moment, V contains the right singular vectors of A.
*           Next, assemble the left singular vector matrix U (M x N).
            IF ( NR .LT. M ) THEN
               CALL CLASET('A', M-NR, NR, CZERO, CZERO, U(NR+1,1), LDU)
               IF ( NR .LT. N1 ) THEN
                  CALL CLASET('A',NR,N1-NR,CZERO,CZERO,U(1,NR+1),LDU)
                  CALL CLASET('A',M-NR,N1-NR,CZERO,CONE, U(NR+1,NR+1),LDU)
               END IF
            END IF
*
*           The Q matrix from the first QRF is built into the left singular
*           matrix U. This applies to all cases.
*
            CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR )

*           The columns of U are normalized. The cost is O(M*N) flops.
            TEMP1 = SQRT(REAL(M)) * EPSLN
            DO 1973 p = 1, NR
               XSC = ONE / SCNRM2( M, U(1,p), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL CSSCAL( M, XSC, U(1,p), 1 )
 1973       CONTINUE
*
*           If the initial QRF is computed with row pivoting, the left
*           singular vectors must be adjusted.
*
            IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
*
         ELSE
*
*        .. the initial matrix A has almost orthogonal columns and
*        the second QRF is not needed
*
            CALL CLACPY( 'U', N, N, A, LDA, CWORK(N+1), N )
            IF ( L2PERT ) THEN
               XSC = SQRT(SMALL)
               DO 5970 p = 2, N
                  CTEMP = XSC * CWORK( N + (p-1)*N + p )
                  DO 5971 q = 1, p - 1
*                     CWORK(N+(q-1)*N+p)=-TEMP1 * ( CWORK(N+(p-1)*N+q) /
*     $                                        ABS(CWORK(N+(p-1)*N+q)) )
                     CWORK(N+(q-1)*N+p)=-CTEMP
 5971             CONTINUE
 5970          CONTINUE
            ELSE
               CALL CLASET( 'L',N-1,N-1,CZERO,CZERO,CWORK(N+2),N )
            END IF
*
            CALL CGESVJ( 'U', 'U', 'N', N, N, CWORK(N+1), N, SVA, N, U, LDU, CWORK(N+N*N+1), LWORK-N-N*N, RWORK, LRWORK, INFO )
*
            SCALEM  = RWORK(1)
            NUMRANK = NINT(RWORK(2))
            DO 6970 p = 1, N
               CALL CCOPY( N, CWORK(N+(p-1)*N+1), 1, U(1,p), 1 )
               CALL CSSCAL( N, SVA(p), CWORK(N+(p-1)*N+1), 1 )
 6970       CONTINUE
*
            CALL CTRSM( 'L', 'U', 'N', 'N', N, N, CONE, A, LDA, CWORK(N+1), N )
            DO 6972 p = 1, N
               CALL CCOPY( N, CWORK(N+p), N, V(IWORK(p),1), LDV )
 6972       CONTINUE
            TEMP1 = SQRT(REAL(N))*EPSLN
            DO 6971 p = 1, N
               XSC = ONE / SCNRM2( N, V(1,p), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL CSSCAL( N, XSC, V(1,p), 1 )
 6971       CONTINUE
*
*           Assemble the left singular vector matrix U (M x N).
*
            IF ( N .LT. M ) THEN
               CALL CLASET( 'A',  M-N, N, CZERO, CZERO, U(N+1,1), LDU )
               IF ( N .LT. N1 ) THEN
                  CALL CLASET('A',N,  N1-N, CZERO, CZERO,  U(1,N+1),LDU)
                  CALL CLASET( 'A',M-N,N1-N, CZERO, CONE,U(N+1,N+1),LDU)
               END IF
            END IF
            CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR )
            TEMP1 = SQRT(REAL(M))*EPSLN
            DO 6973 p = 1, N1
               XSC = ONE / SCNRM2( M, U(1,p), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL CSSCAL( M, XSC, U(1,p), 1 )
 6973       CONTINUE
*
            IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
*
         END IF
*
*        end of the  >> almost orthogonal case <<  in the full SVD
*
         ELSE
*
*        This branch deploys a preconditioned Jacobi SVD with explicitly
*        accumulated rotations. It is included as optional, mainly for
*        experimental purposes. It does perform well, and can also be used.
*        In this implementation, this branch will be automatically activated
*        if the  condition number sigma_max(A) / sigma_min(A) is predicted
*        to be greater than the overflow threshold. This is because the
*        a posteriori computation of the singular vectors assumes robust
*        implementation of BLAS and some LAPACK procedures, capable of working
*        in presence of extreme values, e.g. when the singular values spread from
*        the underflow to the overflow threshold.
*
         DO 7968 p = 1, NR
            CALL CCOPY( N-p+1, A(p,p), LDA, V(p,p), 1 )
            CALL CLACGV( N-p+1, V(p,p), 1 )
 7968    CONTINUE
*
         IF ( L2PERT ) THEN
            XSC = SQRT(SMALL/EPSLN)
            DO 5969 q = 1, NR
               CTEMP = CMPLX(XSC*ABS( V(q,q) ),ZERO)
               DO 5968 p = 1, N
                  IF ( ( p .GT. q ) .AND. ( ABS(V(p,q)) .LE. TEMP1 ) .OR. ( p .LT. q ) )
*     $                V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) )
     $                V(p,q) = CTEMP
                  IF ( p .LT. q ) V(p,q) = - V(p,q)
 5968          CONTINUE
 5969       CONTINUE
         ELSE
            CALL CLASET( 'U', NR-1, NR-1, CZERO, CZERO, V(1,2), LDV )
         END IF
          CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), LWORK-2*N, IERR )
         CALL CLACPY( 'L', N, NR, V, LDV, CWORK(2*N+1), N )
*
         DO 7969 p = 1, NR
            CALL CCOPY( NR-p+1, V(p,p), LDV, U(p,p), 1 )
            CALL CLACGV( NR-p+1, U(p,p), 1 )
 7969    CONTINUE

         IF ( L2PERT ) THEN
            XSC = SQRT(SMALL/EPSLN)
            DO 9970 q = 2, NR
               DO 9971 p = 1, q - 1
                  CTEMP = CMPLX(XSC * MIN(ABS(U(p,p)),ABS(U(q,q))), ZERO)
*                  U(p,q) = - TEMP1 * ( U(q,p) / ABS(U(q,p)) )
                  U(p,q) = - CTEMP
 9971          CONTINUE
 9970       CONTINUE
         ELSE
            CALL CLASET('U', NR-1, NR-1, CZERO, CZERO, U(1,2), LDU )
         END IF
          CALL CGESVJ( 'L', 'U', 'V', NR, NR, U, LDU, SVA, N, V, LDV, CWORK(2*N+N*NR+1), LWORK-2*N-N*NR, RWORK, LRWORK, INFO )
         SCALEM  = RWORK(1)
         NUMRANK = NINT(RWORK(2))

         IF ( NR .LT. N ) THEN
            CALL CLASET( 'A',N-NR,NR,CZERO,CZERO,V(NR+1,1),LDV )
            CALL CLASET( 'A',NR,N-NR,CZERO,CZERO,V(1,NR+1),LDV )
            CALL CLASET( 'A',N-NR,N-NR,CZERO,CONE,V(NR+1,NR+1),LDV )
         END IF
          CALL CUNMQR( 'L','N',N,N,NR,CWORK(2*N+1),N,CWORK(N+1), V,LDV,CWORK(2*N+N*NR+NR+1),LWORK-2*N-N*NR-NR,IERR )
*
*           Permute the rows of V using the (column) permutation from the
*           first QRF. Also, scale the columns to make them unit in
*           Euclidean norm. This applies to all cases.
*
            TEMP1 = SQRT(REAL(N)) * EPSLN
            DO 7972 q = 1, N
               DO 8972 p = 1, N
                  CWORK(2*N+N*NR+NR+IWORK(p)) = V(p,q)
 8972          CONTINUE
               DO 8973 p = 1, N
                  V(p,q) = CWORK(2*N+N*NR+NR+p)
 8973          CONTINUE
               XSC = ONE / SCNRM2( N, V(1,q), 1 )
               IF ( (XSC .LT. (ONE-TEMP1)) .OR. (XSC .GT. (ONE+TEMP1)) ) CALL CSSCAL( N, XSC, V(1,q), 1 )
 7972       CONTINUE
*
*           At this moment, V contains the right singular vectors of A.
*           Next, assemble the left singular vector matrix U (M x N).
*
         IF ( NR .LT. M ) THEN
            CALL CLASET( 'A',  M-NR, NR, CZERO, CZERO, U(NR+1,1), LDU )
            IF ( NR .LT. N1 ) THEN
               CALL CLASET('A',NR,  N1-NR, CZERO, CZERO,  U(1,NR+1),LDU)
               CALL CLASET('A',M-NR,N1-NR, CZERO, CONE,U(NR+1,NR+1),LDU)
            END IF
         END IF
*
         CALL CUNMQR( 'L', 'N', M, N1, N, A, LDA, CWORK, U, LDU, CWORK(N+1), LWORK-N, IERR )
*
            IF ( ROWPIV ) CALL CLASWP( N1, U, LDU, 1, M-1, IWORK(IWOFF+1), -1 )
*
*
         END IF
         IF ( TRANSP ) THEN
*           .. swap U and V because the procedure worked on A^*
            DO 6974 p = 1, N
               CALL CSWAP( N, U(1,p), 1, V(1,p), 1 )
 6974       CONTINUE
         END IF
*
      END IF
*     end of the full SVD
*
*     Undo scaling, if necessary (and possible)
*
      IF ( USCAL2 .LE. (BIG/SVA(1))*USCAL1 ) THEN
         CALL SLASCL( 'G', 0, 0, USCAL1, USCAL2, NR, 1, SVA, N, IERR )
         USCAL1 = ONE
         USCAL2 = ONE
      END IF
*
      IF ( NR .LT. N ) THEN
         DO 3004 p = NR+1, N
            SVA(p) = ZERO
 3004    CONTINUE
      END IF
*
      RWORK(1) = USCAL2 * SCALEM
      RWORK(2) = USCAL1
      IF ( ERREST ) RWORK(3) = SCONDA
      IF ( LSVEC .AND. RSVEC ) THEN
         RWORK(4) = CONDR1
         RWORK(5) = CONDR2
      END IF
      IF ( L2TRAN ) THEN
         RWORK(6) = ENTRA
         RWORK(7) = ENTRAT
      END IF
*
      IWORK(1) = NR
      IWORK(2) = NUMRANK
      IWORK(3) = WARNING
      IF ( TRANSP ) THEN
          IWORK(4) =  1
      ELSE
          IWORK(4) = -1
      END IF

*
      RETURN
*     ..
*     .. END OF CGEJSV
*     ..
      END
*
