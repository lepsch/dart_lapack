      SUBROUTINE SGESVDQ( JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA,
     $                    S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK,
     $                    WORK, LWORK, RWORK, LRWORK, INFO )
*     .. Scalar Arguments ..
      IMPLICIT    NONE
      CHARACTER   JOBA, JOBP, JOBR, JOBU, JOBV
      INTEGER     M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LWORK, LRWORK,
     $            INFO
*     ..
*     .. Array Arguments ..
      REAL        A( LDA, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      REAL        S( * ), RWORK( * )
      INTEGER     IWORK( * )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL        ZERO,         ONE
      PARAMETER ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER     IERR, IWOFF, NR, N1, OPTRATIO, p, q
      INTEGER     LWCON, LWQP3, LWRK_SGELQF, LWRK_SGESVD, LWRK_SGESVD2,
     $            LWRK_SGEQP3,  LWRK_SGEQRF, LWRK_SORMLQ, LWRK_SORMQR,
     $            LWRK_SORMQR2, LWLQF, LWQRF, LWSVD, LWSVD2, LWORQ,
     $            LWORQ2, LWUNLQ, MINWRK, MINWRK2, OPTWRK, OPTWRK2,
     $            IMINWRK, RMINWRK
      LOGICAL     ACCLA,  ACCLM, ACCLH, ASCALED, CONDA, DNTWU,  DNTWV,
     $            LQUERY, LSVC0, LSVEC, ROWPRM,  RSVEC, RTRANS, WNTUA,
     $            WNTUF,  WNTUR, WNTUS, WNTVA,   WNTVR
      REAL        BIG, EPSLN, RTMP, SCONDA, SFMIN
*     ..
*     .. Local Arrays
      REAL        RDUMMY(1)
*     ..
*     .. External Subroutines (BLAS, LAPACK)
      EXTERNAL    SGELQF, SGEQP3, SGEQRF, SGESVD, SLACPY, SLAPMT,
     $            SLASCL, SLASET, SLASWP, SSCAL,  SPOCON, SORMLQ,
     $            SORMQR, XERBLA
*     ..
*     .. External Functions (BLAS, LAPACK)
      LOGICAL    LSAME
      INTEGER    ISAMAX
      REAL        SLANGE, SNRM2, SLAMCH
      EXTERNAL    SLANGE, LSAME, ISAMAX, SNRM2, SLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC   ABS, MAX, MIN, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      WNTUS  = LSAME( JOBU, 'S' ) .OR. LSAME( JOBU, 'U' )
      WNTUR  = LSAME( JOBU, 'R' )
      WNTUA  = LSAME( JOBU, 'A' )
      WNTUF  = LSAME( JOBU, 'F' )
      LSVC0  = WNTUS .OR. WNTUR .OR. WNTUA
      LSVEC  = LSVC0 .OR. WNTUF
      DNTWU  = LSAME( JOBU, 'N' )
*
      WNTVR  = LSAME( JOBV, 'R' )
      WNTVA  = LSAME( JOBV, 'A' ) .OR. LSAME( JOBV, 'V' )
      RSVEC  = WNTVR .OR. WNTVA
      DNTWV  = LSAME( JOBV, 'N' )
*
      ACCLA  = LSAME( JOBA, 'A' )
      ACCLM  = LSAME( JOBA, 'M' )
      CONDA  = LSAME( JOBA, 'E' )
      ACCLH  = LSAME( JOBA, 'H' ) .OR. CONDA
*
      ROWPRM = LSAME( JOBP, 'P' )
      RTRANS = LSAME( JOBR, 'T' )
*
      IF ( ROWPRM ) THEN
         IF ( CONDA ) THEN
            IMINWRK = MAX( 1, N + M - 1 + N )
         ELSE
            IMINWRK = MAX( 1, N + M - 1 )
         END IF
         RMINWRK = MAX( 2, M )
      ELSE
         IF ( CONDA ) THEN
            IMINWRK = MAX( 1, N + N )
         ELSE
            IMINWRK = MAX( 1, N )
         END IF
         RMINWRK = 2
      END IF
      LQUERY = (LIWORK .EQ. -1 .OR. LWORK .EQ. -1 .OR. LRWORK .EQ. -1)
      INFO  = 0
      IF ( .NOT. ( ACCLA .OR. ACCLM .OR. ACCLH ) ) THEN
         INFO = -1
      ELSE IF ( .NOT.( ROWPRM .OR. LSAME( JOBP, 'N' ) ) ) THEN
          INFO = -2
      ELSE IF ( .NOT.( RTRANS .OR. LSAME( JOBR, 'N' ) ) ) THEN
          INFO = -3
      ELSE IF ( .NOT.( LSVEC .OR. DNTWU ) ) THEN
         INFO = -4
      ELSE IF ( WNTUR .AND. WNTVA ) THEN
         INFO = -5
      ELSE IF ( .NOT.( RSVEC .OR. DNTWV )) THEN
         INFO = -5
      ELSE IF ( M.LT.0 ) THEN
         INFO = -6
      ELSE IF ( ( N.LT.0 ) .OR. ( N.GT.M ) ) THEN
         INFO = -7
      ELSE IF ( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF ( LDU.LT.1 .OR. ( LSVC0 .AND. LDU.LT.M ) .OR.
     $       ( WNTUF .AND. LDU.LT.N ) ) THEN
         INFO = -12
      ELSE IF ( LDV.LT.1 .OR. ( RSVEC .AND. LDV.LT.N ) .OR.
     $          ( CONDA .AND. LDV.LT.N ) ) THEN
         INFO = -14
      ELSE IF ( LIWORK .LT. IMINWRK .AND. .NOT. LQUERY ) THEN
         INFO = -17
      END IF
*
*
      IF ( INFO .EQ. 0 ) THEN
*        .. compute the minimal and the optimal workspace lengths
*        [[The expressions for computing the minimal and the optimal
*        values of LWORK are written with a lot of redundancy and
*        can be simplified. However, this detailed form is easier for
*        maintenance and modifications of the code.]]
*
*        .. minimal workspace length for SGEQP3 of an M x N matrix
         LWQP3 = 3 * N + 1
*        .. minimal workspace length for SORMQR to build left singular vectors
         IF ( WNTUS .OR. WNTUR ) THEN
             LWORQ  = MAX( N  , 1 )
         ELSE IF ( WNTUA ) THEN
             LWORQ = MAX( M , 1 )
         END IF
*        .. minimal workspace length for SPOCON of an N x N matrix
         LWCON = 3 * N
*        .. SGESVD of an N x N matrix
         LWSVD = MAX( 5 * N, 1 )
         IF ( LQUERY ) THEN
             CALL SGEQP3( M, N, A, LDA, IWORK, RDUMMY, RDUMMY, -1,
     $           IERR )
             LWRK_SGEQP3 = INT( RDUMMY(1) )
             IF ( WNTUS .OR. WNTUR ) THEN
                 CALL SORMQR( 'L', 'N', M, N, N, A, LDA, RDUMMY, U,
     $                LDU, RDUMMY, -1, IERR )
                 LWRK_SORMQR = INT( RDUMMY(1) )
             ELSE IF ( WNTUA ) THEN
                 CALL SORMQR( 'L', 'N', M, M, N, A, LDA, RDUMMY, U,
     $                LDU, RDUMMY, -1, IERR )
                 LWRK_SORMQR = INT( RDUMMY(1) )
             ELSE
                 LWRK_SORMQR = 0
             END IF
         END IF
         MINWRK = 2
         OPTWRK = 2
         IF ( .NOT. (LSVEC .OR. RSVEC )) THEN
*            .. minimal and optimal sizes of the workspace if
*            only the singular values are requested
             IF ( CONDA ) THEN
                MINWRK = MAX( N+LWQP3, LWCON, LWSVD )
             ELSE
                MINWRK = MAX( N+LWQP3, LWSVD )
             END IF
             IF ( LQUERY ) THEN
                 CALL SGESVD( 'N', 'N', N, N, A, LDA, S, U, LDU,
     $                V, LDV, RDUMMY, -1, IERR )
                 LWRK_SGESVD = INT( RDUMMY(1) )
                 IF ( CONDA ) THEN
                    OPTWRK = MAX( N+LWRK_SGEQP3, N+LWCON, LWRK_SGESVD )
                 ELSE
                    OPTWRK = MAX( N+LWRK_SGEQP3, LWRK_SGESVD )
                 END IF
             END IF
         ELSE IF ( LSVEC .AND. (.NOT.RSVEC) ) THEN
*            .. minimal and optimal sizes of the workspace if the
*            singular values and the left singular vectors are requested
             IF ( CONDA ) THEN
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD, LWORQ )
             ELSE
                 MINWRK = N + MAX( LWQP3, LWSVD, LWORQ )
             END IF
             IF ( LQUERY ) THEN
                IF ( RTRANS ) THEN
                   CALL SGESVD( 'N', 'O', N, N, A, LDA, S, U, LDU,
     $                  V, LDV, RDUMMY, -1, IERR )
                ELSE
                   CALL SGESVD( 'O', 'N', N, N, A, LDA, S, U, LDU,
     $                  V, LDV, RDUMMY, -1, IERR )
                END IF
                LWRK_SGESVD = INT( RDUMMY(1) )
                IF ( CONDA ) THEN
                    OPTWRK = N + MAX( LWRK_SGEQP3, LWCON, LWRK_SGESVD,
     $                               LWRK_SORMQR )
                ELSE
                    OPTWRK = N + MAX( LWRK_SGEQP3, LWRK_SGESVD,
     $                               LWRK_SORMQR )
                END IF
             END IF
         ELSE IF ( RSVEC .AND. (.NOT.LSVEC) ) THEN
*            .. minimal and optimal sizes of the workspace if the
*            singular values and the right singular vectors are requested
             IF ( CONDA ) THEN
                 MINWRK = N + MAX( LWQP3, LWCON, LWSVD )
             ELSE
                 MINWRK = N + MAX( LWQP3, LWSVD )
             END IF
             IF ( LQUERY ) THEN
                 IF ( RTRANS ) THEN
                     CALL SGESVD( 'O', 'N', N, N, A, LDA, S, U, LDU,
     $                    V, LDV, RDUMMY, -1, IERR )
                 ELSE
                     CALL SGESVD( 'N', 'O', N, N, A, LDA, S, U, LDU,
     $                    V, LDV, RDUMMY, -1, IERR )
                 END IF
                 LWRK_SGESVD = INT( RDUMMY(1) )
                 IF ( CONDA ) THEN
                     OPTWRK = N + MAX( LWRK_SGEQP3, LWCON, LWRK_SGESVD )
                 ELSE
                     OPTWRK = N + MAX( LWRK_SGEQP3, LWRK_SGESVD )
                 END IF
             END IF
         ELSE
*            .. minimal and optimal sizes of the workspace if the
*            full SVD is requested
             IF ( RTRANS ) THEN
                 MINWRK = MAX( LWQP3, LWSVD, LWORQ )
                 IF ( CONDA ) MINWRK = MAX( MINWRK, LWCON )
                 MINWRK = MINWRK + N
                 IF ( WNTVA ) THEN
*                   .. minimal workspace length for N x N/2 SGEQRF
                    LWQRF  = MAX( N/2, 1 )
*                   .. minimal workspace length for N/2 x N/2 SGESVD
                    LWSVD2 = MAX( 5 * (N/2), 1 )
                    LWORQ2 = MAX( N, 1 )
                    MINWRK2 = MAX( LWQP3, N/2+LWQRF, N/2+LWSVD2,
     $                        N/2+LWORQ2, LWORQ )
                    IF ( CONDA ) MINWRK2 = MAX( MINWRK2, LWCON )
                    MINWRK2 = N + MINWRK2
                    MINWRK = MAX( MINWRK, MINWRK2 )
                 END IF
             ELSE
                 MINWRK = MAX( LWQP3, LWSVD, LWORQ )
                 IF ( CONDA ) MINWRK = MAX( MINWRK, LWCON )
                 MINWRK = MINWRK + N
                 IF ( WNTVA ) THEN
*                   .. minimal workspace length for N/2 x N SGELQF
                    LWLQF  = MAX( N/2, 1 )
                    LWSVD2 = MAX( 5 * (N/2), 1 )
                    LWUNLQ = MAX( N , 1 )
                    MINWRK2 = MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2,
     $                        N/2+LWUNLQ, LWORQ )
                    IF ( CONDA ) MINWRK2 = MAX( MINWRK2, LWCON )
                    MINWRK2 = N + MINWRK2
                    MINWRK = MAX( MINWRK, MINWRK2 )
                 END IF
             END IF
             IF ( LQUERY ) THEN
                IF ( RTRANS ) THEN
                   CALL SGESVD( 'O', 'A', N, N, A, LDA, S, U, LDU,
     $                  V, LDV, RDUMMY, -1, IERR )
                   LWRK_SGESVD = INT( RDUMMY(1) )
                   OPTWRK = MAX(LWRK_SGEQP3,LWRK_SGESVD,LWRK_SORMQR)
                   IF ( CONDA ) OPTWRK = MAX( OPTWRK, LWCON )
                   OPTWRK = N + OPTWRK
                   IF ( WNTVA ) THEN
                       CALL SGEQRF(N,N/2,U,LDU,RDUMMY,RDUMMY,-1,IERR)
                       LWRK_SGEQRF = INT( RDUMMY(1) )
                       CALL SGESVD( 'S', 'O', N/2,N/2, V,LDV, S, U,LDU,
     $                      V, LDV, RDUMMY, -1, IERR )
                       LWRK_SGESVD2 = INT( RDUMMY(1) )
                       CALL SORMQR( 'R', 'C', N, N, N/2, U, LDU, RDUMMY,
     $                      V, LDV, RDUMMY, -1, IERR )
                       LWRK_SORMQR2 = INT( RDUMMY(1) )
                       OPTWRK2 = MAX( LWRK_SGEQP3, N/2+LWRK_SGEQRF,
     $                           N/2+LWRK_SGESVD2, N/2+LWRK_SORMQR2 )
                       IF ( CONDA ) OPTWRK2 = MAX( OPTWRK2, LWCON )
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   END IF
                ELSE
                   CALL SGESVD( 'S', 'O', N, N, A, LDA, S, U, LDU,
     $                  V, LDV, RDUMMY, -1, IERR )
                   LWRK_SGESVD = INT( RDUMMY(1) )
                   OPTWRK = MAX(LWRK_SGEQP3,LWRK_SGESVD,LWRK_SORMQR)
                   IF ( CONDA ) OPTWRK = MAX( OPTWRK, LWCON )
                   OPTWRK = N + OPTWRK
                   IF ( WNTVA ) THEN
                      CALL SGELQF(N/2,N,U,LDU,RDUMMY,RDUMMY,-1,IERR)
                      LWRK_SGELQF = INT( RDUMMY(1) )
                      CALL SGESVD( 'S','O', N/2,N/2, V, LDV, S, U, LDU,
     $                     V, LDV, RDUMMY, -1, IERR )
                      LWRK_SGESVD2 = INT( RDUMMY(1) )
                      CALL SORMLQ( 'R', 'N', N, N, N/2, U, LDU, RDUMMY,
     $                     V, LDV, RDUMMY,-1,IERR )
                      LWRK_SORMLQ = INT( RDUMMY(1) )
                      OPTWRK2 = MAX( LWRK_SGEQP3, N/2+LWRK_SGELQF,
     $                           N/2+LWRK_SGESVD2, N/2+LWRK_SORMLQ )
                       IF ( CONDA ) OPTWRK2 = MAX( OPTWRK2, LWCON )
                       OPTWRK2 = N + OPTWRK2
                       OPTWRK = MAX( OPTWRK, OPTWRK2 )
                   END IF
                END IF
             END IF
         END IF
*
         MINWRK = MAX( 2, MINWRK )
         OPTWRK = MAX( 2, OPTWRK )
         IF ( LWORK .LT. MINWRK .AND. (.NOT.LQUERY) ) INFO = -19
*
      END IF
*
      IF (INFO .EQ. 0 .AND. LRWORK .LT. RMINWRK .AND. .NOT. LQUERY) THEN
         INFO = -21
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGESVDQ', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
*
*     Return optimal workspace
*
          IWORK(1) = IMINWRK
          WORK(1) = OPTWRK
          WORK(2) = MINWRK
          RWORK(1) = RMINWRK
          RETURN
      END IF
*
*     Quick return if the matrix is void.
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) THEN
*     .. all output is void.
         RETURN
      END IF
*
      BIG = SLAMCH('O')
      ASCALED = .FALSE.
      IWOFF = 1
      IF ( ROWPRM ) THEN
            IWOFF = M
*           .. reordering the rows in decreasing sequence in the
*           ell-infinity norm - this enhances numerical robustness in
*           the case of differently scaled rows.
            DO 1904 p = 1, M
*               RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
*               [[SLANGE will return NaN if an entry of the p-th row is Nan]]
                RWORK(p) = SLANGE( 'M', 1, N, A(p,1), LDA, RDUMMY )
*               .. check for NaN's and Inf's
                IF ( ( RWORK(p) .NE. RWORK(p) ) .OR.
     $               ( (RWORK(p)*ZERO) .NE. ZERO ) ) THEN
                    INFO = -8
                    CALL XERBLA( 'SGESVDQ', -INFO )
                    RETURN
                END IF
 1904       CONTINUE
            DO 1952 p = 1, M - 1
            q = ISAMAX( M-p+1, RWORK(p), 1 ) + p - 1
            IWORK(N+p) = q
            IF ( p .NE. q ) THEN
               RTMP     = RWORK(p)
               RWORK(p) = RWORK(q)
               RWORK(q) = RTMP
            END IF
 1952       CONTINUE
*
            IF ( RWORK(1) .EQ. ZERO ) THEN
*              Quick return: A is the M x N zero matrix.
               NUMRANK = 0
               CALL SLASET( 'G', N, 1, ZERO, ZERO, S, N )
               IF ( WNTUS ) CALL SLASET('G', M, N, ZERO, ONE, U, LDU)
               IF ( WNTUA ) CALL SLASET('G', M, M, ZERO, ONE, U, LDU)
               IF ( WNTVA ) CALL SLASET('G', N, N, ZERO, ONE, V, LDV)
               IF ( WNTUF ) THEN
                   CALL SLASET( 'G', N, 1, ZERO, ZERO, WORK, N )
                   CALL SLASET( 'G', M, N, ZERO,  ONE, U, LDU )
               END IF
               DO 5001 p = 1, N
                   IWORK(p) = p
 5001          CONTINUE
               IF ( ROWPRM ) THEN
                   DO 5002 p = N + 1, N + M - 1
                       IWORK(p) = p - N
 5002              CONTINUE
               END IF
               IF ( CONDA ) RWORK(1) = -1
               RWORK(2) = -1
               RETURN
            END IF
*
            IF ( RWORK(1) .GT. BIG / SQRT(REAL(M)) ) THEN
*               .. to prevent overflow in the QR factorization, scale the
*               matrix by 1/sqrt(M) if too large entry detected
                CALL SLASCL('G',0,0,SQRT(REAL(M)),ONE, M,N, A,LDA, IERR)
                ASCALED = .TRUE.
            END IF
            CALL SLASWP( N, A, LDA, 1, M-1, IWORK(N+1), 1 )
      END IF
*
*    .. At this stage, preemptive scaling is done only to avoid column
*    norms overflows during the QR factorization. The SVD procedure should
*    have its own scaling to save the singular values from overflows and
*    underflows. That depends on the SVD procedure.
*
      IF ( .NOT.ROWPRM ) THEN
          RTMP = SLANGE( 'M', M, N, A, LDA, RDUMMY )
          IF ( ( RTMP .NE. RTMP ) .OR.
     $         ( (RTMP*ZERO) .NE. ZERO ) ) THEN
               INFO = -8
               CALL XERBLA( 'SGESVDQ', -INFO )
               RETURN
          END IF
          IF ( RTMP .GT. BIG / SQRT(REAL(M)) ) THEN
*             .. to prevent overflow in the QR factorization, scale the
*             matrix by 1/sqrt(M) if too large entry detected
              CALL SLASCL('G',0,0, SQRT(REAL(M)),ONE, M,N, A,LDA, IERR)
              ASCALED = .TRUE.
          END IF
      END IF
*
*     .. QR factorization with column pivoting
*
*     A * P = Q * [ R ]
*                 [ 0 ]
*
      DO 1963 p = 1, N
*        .. all columns are free columns
         IWORK(p) = 0
 1963 CONTINUE
      CALL SGEQP3( M, N, A, LDA, IWORK, WORK, WORK(N+1), LWORK-N,
     $      IERR )
*
*    If the user requested accuracy level allows truncation in the
*    computed upper triangular factor, the matrix R is examined and,
*    if possible, replaced with its leading upper trapezoidal part.
*
      EPSLN = SLAMCH('E')
      SFMIN = SLAMCH('S')
*     SMALL = SFMIN / EPSLN
      NR = N
*
      IF ( ACCLA ) THEN
*
*        Standard absolute error bound suffices. All sigma_i with
*        sigma_i < N*EPS*||A||_F are flushed to zero. This is an
*        aggressive enforcement of lower numerical rank by introducing a
*        backward error of the order of N*EPS*||A||_F.
         NR = 1
         RTMP = SQRT(REAL(N))*EPSLN
         DO 3001 p = 2, N
            IF ( ABS(A(p,p)) .LT. (RTMP*ABS(A(1,1))) ) GO TO 3002
               NR = NR + 1
 3001    CONTINUE
 3002    CONTINUE
*
      ELSEIF ( ACCLM ) THEN
*        .. similarly as above, only slightly more gentle (less aggressive).
*        Sudden drop on the diagonal of R is used as the criterion for being
*        close-to-rank-deficient. The threshold is set to EPSLN=SLAMCH('E').
*        [[This can be made more flexible by replacing this hard-coded value
*        with a user specified threshold.]] Also, the values that underflow
*        will be truncated.
         NR = 1
         DO 3401 p = 2, N
            IF ( ( ABS(A(p,p)) .LT. (EPSLN*ABS(A(p-1,p-1))) ) .OR.
     $           ( ABS(A(p,p)) .LT. SFMIN ) ) GO TO 3402
            NR = NR + 1
 3401    CONTINUE
 3402    CONTINUE
*
      ELSE
*        .. RRQR not authorized to determine numerical rank except in the
*        obvious case of zero pivots.
*        .. inspect R for exact zeros on the diagonal;
*        R(i,i)=0 => R(i:N,i:N)=0.
         NR = 1
         DO 3501 p = 2, N
            IF ( ABS(A(p,p)) .EQ. ZERO ) GO TO 3502
            NR = NR + 1
 3501    CONTINUE
 3502    CONTINUE
*
         IF ( CONDA ) THEN
*           Estimate the scaled condition number of A. Use the fact that it is
*           the same as the scaled condition number of R.
*              .. V is used as workspace
               CALL SLACPY( 'U', N, N, A, LDA, V, LDV )
*              Only the leading NR x NR submatrix of the triangular factor
*              is considered. Only if NR=N will this give a reliable error
*              bound. However, even for NR < N, this can be used on an
*              expert level and obtain useful information in the sense of
*              perturbation theory.
               DO 3053 p = 1, NR
                  RTMP = SNRM2( p, V(1,p), 1 )
                  CALL SSCAL( p, ONE/RTMP, V(1,p), 1 )
 3053          CONTINUE
               IF ( .NOT. ( LSVEC .OR. RSVEC ) ) THEN
                   CALL SPOCON( 'U', NR, V, LDV, ONE, RTMP,
     $                  WORK, IWORK(N+IWOFF), IERR )
               ELSE
                   CALL SPOCON( 'U', NR, V, LDV, ONE, RTMP,
     $                  WORK(N+1), IWORK(N+IWOFF), IERR )
               END IF
               SCONDA = ONE / SQRT(RTMP)
*           For NR=N, SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1),
*           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
*           See the reference [1] for more details.
         END IF
*
      ENDIF
*
      IF ( WNTUR ) THEN
          N1 = NR
      ELSE IF ( WNTUS .OR. WNTUF) THEN
          N1 = N
      ELSE IF ( WNTUA ) THEN
          N1 = M
      END IF
*
      IF ( .NOT. ( RSVEC .OR. LSVEC ) ) THEN
*.......................................................................
*        .. only the singular values are requested
*.......................................................................
         IF ( RTRANS ) THEN
*
*         .. compute the singular values of R**T = [A](1:NR,1:N)**T
*           .. set the lower triangle of [A] to [A](1:NR,1:N)**T and
*           the upper triangle of [A] to zero.
            DO 1146 p = 1, MIN( N, NR )
               DO 1147 q = p + 1, N
                  A(q,p) = A(p,q)
                  IF ( q .LE. NR ) A(p,q) = ZERO
 1147          CONTINUE
 1146       CONTINUE
*
            CALL SGESVD( 'N', 'N', N, NR, A, LDA, S, U, LDU,
     $           V, LDV, WORK, LWORK, INFO )
*
         ELSE
*
*           .. compute the singular values of R = [A](1:NR,1:N)
*
            IF ( NR .GT. 1 )
     $          CALL SLASET( 'L', NR-1,NR-1, ZERO,ZERO, A(2,1), LDA )
            CALL SGESVD( 'N', 'N', NR, N, A, LDA, S, U, LDU,
     $           V, LDV, WORK, LWORK, INFO )
*
         END IF
*
      ELSE IF ( LSVEC .AND. ( .NOT. RSVEC) ) THEN
*.......................................................................
*       .. the singular values and the left singular vectors requested
*.......................................................................""""""""
         IF ( RTRANS ) THEN
*            .. apply SGESVD to R**T
*            .. copy R**T into [U] and overwrite [U] with the right singular
*            vectors of R
            DO 1192 p = 1, NR
               DO 1193 q = p, N
                  U(q,p) = A(p,q)
 1193          CONTINUE
 1192       CONTINUE
            IF ( NR .GT. 1 )
     $          CALL SLASET( 'U', NR-1,NR-1, ZERO,ZERO, U(1,2), LDU )
*           .. the left singular vectors not computed, the NR right singular
*           vectors overwrite [U](1:NR,1:NR) as transposed. These
*           will be pre-multiplied by Q to build the left singular vectors of A.
               CALL SGESVD( 'N', 'O', N, NR, U, LDU, S, U, LDU,
     $              U, LDU, WORK(N+1), LWORK-N, INFO )
*
               DO 1119 p = 1, NR
                   DO 1120 q = p + 1, NR
                      RTMP   = U(q,p)
                      U(q,p) = U(p,q)
                      U(p,q) = RTMP
 1120              CONTINUE
 1119          CONTINUE
*
         ELSE
*            .. apply SGESVD to R
*            .. copy R into [U] and overwrite [U] with the left singular vectors
             CALL SLACPY( 'U', NR, N, A, LDA, U, LDU )
             IF ( NR .GT. 1 )
     $         CALL SLASET( 'L', NR-1, NR-1, ZERO, ZERO, U(2,1), LDU )
*            .. the right singular vectors not computed, the NR left singular
*            vectors overwrite [U](1:NR,1:NR)
                CALL SGESVD( 'O', 'N', NR, N, U, LDU, S, U, LDU,
     $               V, LDV, WORK(N+1), LWORK-N, INFO )
*               .. now [U](1:NR,1:NR) contains the NR left singular vectors of
*               R. These will be pre-multiplied by Q to build the left singular
*               vectors of A.
         END IF
*
*           .. assemble the left singular vector matrix U of dimensions
*              (M x NR) or (M x N) or (M x M).
         IF ( ( NR .LT. M ) .AND. ( .NOT.WNTUF ) ) THEN
             CALL SLASET('A', M-NR, NR, ZERO, ZERO, U(NR+1,1), LDU)
             IF ( NR .LT. N1 ) THEN
                CALL SLASET( 'A',NR,N1-NR,ZERO,ZERO,U(1,NR+1), LDU )
                CALL SLASET( 'A',M-NR,N1-NR,ZERO,ONE,
     $               U(NR+1,NR+1), LDU )
             END IF
         END IF
*
*           The Q matrix from the first QRF is built into the left singular
*           vectors matrix U.
*
         IF ( .NOT.WNTUF )
     $       CALL SORMQR( 'L', 'N', M, N1, N, A, LDA, WORK, U,
     $            LDU, WORK(N+1), LWORK-N, IERR )
         IF ( ROWPRM .AND. .NOT.WNTUF )
     $          CALL SLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 )
*
      ELSE IF ( RSVEC .AND. ( .NOT. LSVEC ) ) THEN
*.......................................................................
*       .. the singular values and the right singular vectors requested
*.......................................................................
          IF ( RTRANS ) THEN
*            .. apply SGESVD to R**T
*            .. copy R**T into V and overwrite V with the left singular vectors
            DO 1165 p = 1, NR
               DO 1166 q = p, N
                  V(q,p) = (A(p,q))
 1166          CONTINUE
 1165       CONTINUE
            IF ( NR .GT. 1 )
     $          CALL SLASET( 'U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV )
*           .. the left singular vectors of R**T overwrite V, the right singular
*           vectors not computed
            IF ( WNTVR .OR. ( NR .EQ. N ) ) THEN
               CALL SGESVD( 'O', 'N', N, NR, V, LDV, S, U, LDU,
     $              U, LDU, WORK(N+1), LWORK-N, INFO )
*
               DO 1121 p = 1, NR
                   DO 1122 q = p + 1, NR
                      RTMP   = V(q,p)
                      V(q,p) = V(p,q)
                      V(p,q) = RTMP
 1122              CONTINUE
 1121          CONTINUE
*
               IF ( NR .LT. N ) THEN
                   DO 1103 p = 1, NR
                      DO 1104 q = NR + 1, N
                          V(p,q) = V(q,p)
 1104                 CONTINUE
 1103              CONTINUE
               END IF
               CALL SLAPMT( .FALSE., NR, N, V, LDV, IWORK )
            ELSE
*               .. need all N right singular vectors and NR < N
*               [!] This is simple implementation that augments [V](1:N,1:NR)
*               by padding a zero block. In the case NR << N, a more efficient
*               way is to first use the QR factorization. For more details
*               how to implement this, see the " FULL SVD " branch.
                CALL SLASET('G', N, N-NR, ZERO, ZERO, V(1,NR+1), LDV)
                CALL SGESVD( 'O', 'N', N, N, V, LDV, S, U, LDU,
     $               U, LDU, WORK(N+1), LWORK-N, INFO )
*
                DO 1123 p = 1, N
                   DO 1124 q = p + 1, N
                      RTMP   = V(q,p)
                      V(q,p) = V(p,q)
                      V(p,q) = RTMP
 1124              CONTINUE
 1123           CONTINUE
                CALL SLAPMT( .FALSE., N, N, V, LDV, IWORK )
            END IF
*
          ELSE
*            .. aply SGESVD to R
*            .. copy R into V and overwrite V with the right singular vectors
             CALL SLACPY( 'U', NR, N, A, LDA, V, LDV )
             IF ( NR .GT. 1 )
     $         CALL SLASET( 'L', NR-1, NR-1, ZERO, ZERO, V(2,1), LDV )
*            .. the right singular vectors overwrite V, the NR left singular
*            vectors stored in U(1:NR,1:NR)
             IF ( WNTVR .OR. ( NR .EQ. N ) ) THEN
                CALL SGESVD( 'N', 'O', NR, N, V, LDV, S, U, LDU,
     $               V, LDV, WORK(N+1), LWORK-N, INFO )
                CALL SLAPMT( .FALSE., NR, N, V, LDV, IWORK )
*               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
             ELSE
*               .. need all N right singular vectors and NR < N
*               [!] This is simple implementation that augments [V](1:NR,1:N)
*               by padding a zero block. In the case NR << N, a more efficient
*               way is to first use the LQ factorization. For more details
*               how to implement this, see the " FULL SVD " branch.
                 CALL SLASET('G', N-NR, N, ZERO,ZERO, V(NR+1,1), LDV)
                 CALL SGESVD( 'N', 'O', N, N, V, LDV, S, U, LDU,
     $                V, LDV, WORK(N+1), LWORK-N, INFO )
                 CALL SLAPMT( .FALSE., N, N, V, LDV, IWORK )
             END IF
*            .. now [V] contains the transposed matrix of the right singular
*            vectors of A.
          END IF
*
      ELSE
*.......................................................................
*       .. FULL SVD requested
*.......................................................................
         IF ( RTRANS ) THEN
*
*            .. apply SGESVD to R**T [[this option is left for R&D&T]]
*
            IF ( WNTVR .OR. ( NR .EQ. N ) ) THEN
*            .. copy R**T into [V] and overwrite [V] with the left singular
*            vectors of R**T
            DO 1168 p = 1, NR
               DO 1169 q = p, N
                  V(q,p) = A(p,q)
 1169          CONTINUE
 1168       CONTINUE
            IF ( NR .GT. 1 )
     $          CALL SLASET( 'U', NR-1,NR-1, ZERO,ZERO, V(1,2), LDV )
*
*           .. the left singular vectors of R**T overwrite [V], the NR right
*           singular vectors of R**T stored in [U](1:NR,1:NR) as transposed
               CALL SGESVD( 'O', 'A', N, NR, V, LDV, S, V, LDV,
     $              U, LDU, WORK(N+1), LWORK-N, INFO )
*              .. assemble V
               DO 1115 p = 1, NR
                  DO 1116 q = p + 1, NR
                     RTMP   = V(q,p)
                     V(q,p) = V(p,q)
                     V(p,q) = RTMP
 1116             CONTINUE
 1115          CONTINUE
               IF ( NR .LT. N ) THEN
                   DO 1101 p = 1, NR
                      DO 1102 q = NR+1, N
                         V(p,q) = V(q,p)
 1102                 CONTINUE
 1101              CONTINUE
               END IF
               CALL SLAPMT( .FALSE., NR, N, V, LDV, IWORK )
*
                DO 1117 p = 1, NR
                   DO 1118 q = p + 1, NR
                      RTMP   = U(q,p)
                      U(q,p) = U(p,q)
                      U(p,q) = RTMP
 1118              CONTINUE
 1117           CONTINUE
*
                IF ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) THEN
                  CALL SLASET('A', M-NR,NR, ZERO,ZERO, U(NR+1,1), LDU)
                  IF ( NR .LT. N1 ) THEN
                     CALL SLASET('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU)
                     CALL SLASET( 'A',M-NR,N1-NR,ZERO,ONE,
     $                    U(NR+1,NR+1), LDU )
                  END IF
               END IF
*
            ELSE
*               .. need all N right singular vectors and NR < N
*            .. copy R**T into [V] and overwrite [V] with the left singular
*            vectors of R**T
*               [[The optimal ratio N/NR for using QRF instead of padding
*                 with zeros. Here hard coded to 2; it must be at least
*                 two due to work space constraints.]]
*               OPTRATIO = ILAENV(6, 'SGESVD', 'S' // 'O', NR,N,0,0)
*               OPTRATIO = MAX( OPTRATIO, 2 )
                OPTRATIO = 2
                IF ( OPTRATIO*NR .GT. N ) THEN
                   DO 1198 p = 1, NR
                      DO 1199 q = p, N
                         V(q,p) = A(p,q)
 1199                 CONTINUE
 1198              CONTINUE
                   IF ( NR .GT. 1 )
     $             CALL SLASET('U',NR-1,NR-1, ZERO,ZERO, V(1,2),LDV)
*
                   CALL SLASET('A',N,N-NR,ZERO,ZERO,V(1,NR+1),LDV)
                   CALL SGESVD( 'O', 'A', N, N, V, LDV, S, V, LDV,
     $                  U, LDU, WORK(N+1), LWORK-N, INFO )
*
                   DO 1113 p = 1, N
                      DO 1114 q = p + 1, N
                         RTMP   = V(q,p)
                         V(q,p) = V(p,q)
                         V(p,q) = RTMP
 1114                 CONTINUE
 1113              CONTINUE
                   CALL SLAPMT( .FALSE., N, N, V, LDV, IWORK )
*              .. assemble the left singular vector matrix U of dimensions
*              (M x N1), i.e. (M x N) or (M x M).
*
                   DO 1111 p = 1, N
                      DO 1112 q = p + 1, N
                         RTMP   = U(q,p)
                         U(q,p) = U(p,q)
                         U(p,q) = RTMP
 1112                 CONTINUE
 1111              CONTINUE
*
                   IF ( ( N .LT. M ) .AND. .NOT.(WNTUF)) THEN
                      CALL SLASET('A',M-N,N,ZERO,ZERO,U(N+1,1),LDU)
                      IF ( N .LT. N1 ) THEN
                        CALL SLASET('A',N,N1-N,ZERO,ZERO,U(1,N+1),LDU)
                        CALL SLASET('A',M-N,N1-N,ZERO,ONE,
     $                       U(N+1,N+1), LDU )
                      END IF
                   END IF
                ELSE
*                  .. copy R**T into [U] and overwrite [U] with the right
*                  singular vectors of R
                   DO 1196 p = 1, NR
                      DO 1197 q = p, N
                         U(q,NR+p) = A(p,q)
 1197                 CONTINUE
 1196              CONTINUE
                   IF ( NR .GT. 1 )
     $             CALL SLASET('U',NR-1,NR-1,ZERO,ZERO,U(1,NR+2),LDU)
                   CALL SGEQRF( N, NR, U(1,NR+1), LDU, WORK(N+1),
     $                  WORK(N+NR+1), LWORK-N-NR, IERR )
                   DO 1143 p = 1, NR
                       DO 1144 q = 1, N
                           V(q,p) = U(p,NR+q)
 1144                  CONTINUE
 1143              CONTINUE
                  CALL SLASET('U',NR-1,NR-1,ZERO,ZERO,V(1,2),LDV)
                  CALL SGESVD( 'S', 'O', NR, NR, V, LDV, S, U, LDU,
     $                 V,LDV, WORK(N+NR+1),LWORK-N-NR, INFO )
                  CALL SLASET('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV)
                  CALL SLASET('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV)
                  CALL SLASET('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV)
                  CALL SORMQR('R','C', N, N, NR, U(1,NR+1), LDU,
     $                 WORK(N+1),V,LDV,WORK(N+NR+1),LWORK-N-NR,IERR)
                  CALL SLAPMT( .FALSE., N, N, V, LDV, IWORK )
*                 .. assemble the left singular vector matrix U of dimensions
*                 (M x NR) or (M x N) or (M x M).
                  IF ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) THEN
                     CALL SLASET('A',M-NR,NR,ZERO,ZERO,U(NR+1,1),LDU)
                     IF ( NR .LT. N1 ) THEN
                     CALL SLASET('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU)
                     CALL SLASET( 'A',M-NR,N1-NR,ZERO,ONE,
     $                    U(NR+1,NR+1),LDU)
                     END IF
                  END IF
                END IF
            END IF
*
         ELSE
*
*            .. apply SGESVD to R [[this is the recommended option]]
*
             IF ( WNTVR .OR. ( NR .EQ. N ) ) THEN
*                .. copy R into [V] and overwrite V with the right singular vectors
                 CALL SLACPY( 'U', NR, N, A, LDA, V, LDV )
                IF ( NR .GT. 1 )
     $          CALL SLASET( 'L', NR-1,NR-1, ZERO,ZERO, V(2,1), LDV )
*               .. the right singular vectors of R overwrite [V], the NR left
*               singular vectors of R stored in [U](1:NR,1:NR)
                CALL SGESVD( 'S', 'O', NR, N, V, LDV, S, U, LDU,
     $               V, LDV, WORK(N+1), LWORK-N, INFO )
                CALL SLAPMT( .FALSE., NR, N, V, LDV, IWORK )
*               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
*               .. assemble the left singular vector matrix U of dimensions
*              (M x NR) or (M x N) or (M x M).
               IF ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) THEN
                  CALL SLASET('A', M-NR,NR, ZERO,ZERO, U(NR+1,1), LDU)
                  IF ( NR .LT. N1 ) THEN
                     CALL SLASET('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU)
                     CALL SLASET( 'A',M-NR,N1-NR,ZERO,ONE,
     $                    U(NR+1,NR+1), LDU )
                  END IF
               END IF
*
             ELSE
*              .. need all N right singular vectors and NR < N
*              .. the requested number of the left singular vectors
*               is then N1 (N or M)
*               [[The optimal ratio N/NR for using LQ instead of padding
*                 with zeros. Here hard coded to 2; it must be at least
*                 two due to work space constraints.]]
*               OPTRATIO = ILAENV(6, 'SGESVD', 'S' // 'O', NR,N,0,0)
*               OPTRATIO = MAX( OPTRATIO, 2 )
               OPTRATIO = 2
               IF ( OPTRATIO * NR .GT. N ) THEN
                  CALL SLACPY( 'U', NR, N, A, LDA, V, LDV )
                  IF ( NR .GT. 1 )
     $            CALL SLASET('L', NR-1,NR-1, ZERO,ZERO, V(2,1),LDV)
*              .. the right singular vectors of R overwrite [V], the NR left
*                 singular vectors of R stored in [U](1:NR,1:NR)
                  CALL SLASET('A', N-NR,N, ZERO,ZERO, V(NR+1,1),LDV)
                  CALL SGESVD( 'S', 'O', N, N, V, LDV, S, U, LDU,
     $                 V, LDV, WORK(N+1), LWORK-N, INFO )
                  CALL SLAPMT( .FALSE., N, N, V, LDV, IWORK )
*                 .. now [V] contains the transposed matrix of the right
*                 singular vectors of A. The leading N left singular vectors
*                 are in [U](1:N,1:N)
*                 .. assemble the left singular vector matrix U of dimensions
*                 (M x N1), i.e. (M x N) or (M x M).
                  IF ( ( N .LT. M ) .AND. .NOT.(WNTUF)) THEN
                      CALL SLASET('A',M-N,N,ZERO,ZERO,U(N+1,1),LDU)
                      IF ( N .LT. N1 ) THEN
                        CALL SLASET('A',N,N1-N,ZERO,ZERO,U(1,N+1),LDU)
                        CALL SLASET( 'A',M-N,N1-N,ZERO,ONE,
     $                       U(N+1,N+1), LDU )
                      END IF
                  END IF
               ELSE
                  CALL SLACPY( 'U', NR, N, A, LDA, U(NR+1,1), LDU )
                  IF ( NR .GT. 1 )
     $            CALL SLASET('L',NR-1,NR-1,ZERO,ZERO,U(NR+2,1),LDU)
                  CALL SGELQF( NR, N, U(NR+1,1), LDU, WORK(N+1),
     $                 WORK(N+NR+1), LWORK-N-NR, IERR )
                  CALL SLACPY('L',NR,NR,U(NR+1,1),LDU,V,LDV)
                  IF ( NR .GT. 1 )
     $            CALL SLASET('U',NR-1,NR-1,ZERO,ZERO,V(1,2),LDV)
                  CALL SGESVD( 'S', 'O', NR, NR, V, LDV, S, U, LDU,
     $                 V, LDV, WORK(N+NR+1), LWORK-N-NR, INFO )
                  CALL SLASET('A',N-NR,NR,ZERO,ZERO,V(NR+1,1),LDV)
                  CALL SLASET('A',NR,N-NR,ZERO,ZERO,V(1,NR+1),LDV)
                  CALL SLASET('A',N-NR,N-NR,ZERO,ONE,V(NR+1,NR+1),LDV)
                  CALL SORMLQ('R','N',N,N,NR,U(NR+1,1),LDU,WORK(N+1),
     $                 V, LDV, WORK(N+NR+1),LWORK-N-NR,IERR)
                  CALL SLAPMT( .FALSE., N, N, V, LDV, IWORK )
*               .. assemble the left singular vector matrix U of dimensions
*              (M x NR) or (M x N) or (M x M).
                  IF ( ( NR .LT. M ) .AND. .NOT.(WNTUF)) THEN
                     CALL SLASET('A',M-NR,NR,ZERO,ZERO,U(NR+1,1),LDU)
                     IF ( NR .LT. N1 ) THEN
                     CALL SLASET('A',NR,N1-NR,ZERO,ZERO,U(1,NR+1),LDU)
                     CALL SLASET( 'A',M-NR,N1-NR,ZERO,ONE,
     $                    U(NR+1,NR+1), LDU )
                     END IF
                  END IF
               END IF
             END IF
*        .. end of the "R**T or R" branch
         END IF
*
*           The Q matrix from the first QRF is built into the left singular
*           vectors matrix U.
*
         IF ( .NOT. WNTUF )
     $       CALL SORMQR( 'L', 'N', M, N1, N, A, LDA, WORK, U,
     $            LDU, WORK(N+1), LWORK-N, IERR )
         IF ( ROWPRM .AND. .NOT.WNTUF )
     $          CALL SLASWP( N1, U, LDU, 1, M-1, IWORK(N+1), -1 )
*
*     ... end of the "full SVD" branch
      END IF
*
*     Check whether some singular values are returned as zeros, e.g.
*     due to underflow, and update the numerical rank.
      p = NR
      DO 4001 q = p, 1, -1
          IF ( S(q) .GT. ZERO ) GO TO 4002
          NR = NR - 1
 4001 CONTINUE
 4002 CONTINUE
*
*     .. if numerical rank deficiency is detected, the truncated
*     singular values are set to zero.
      IF ( NR .LT. N ) CALL SLASET( 'G', N-NR,1, ZERO,ZERO, S(NR+1), N )
*     .. undo scaling; this may cause overflow in the largest singular
*     values.
      IF ( ASCALED )
     $   CALL SLASCL( 'G',0,0, ONE,SQRT(REAL(M)), NR,1, S, N, IERR )
      IF ( CONDA ) RWORK(1) = SCONDA
      RWORK(2) = p - NR
*     .. p-NR is the number of singular values that are computed as
*     exact zeros in SGESVD() applied to the (possibly truncated)
*     full row rank triangular (trapezoidal) factor of A.
      NUMRANK = NR
*
      RETURN
*
*     End of SGESVDQ
*
      END