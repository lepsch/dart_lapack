      SUBROUTINE CTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,
     $                   LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK,
     $                   IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, JOB
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      REAL               DIF( * ), S( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      INTEGER            IDIFJB
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, IDIFJB = 3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SOMCON, WANTBH, WANTDF, WANTS
      INTEGER            I, IERR, IFST, ILST, K, KS, LWMIN, N1, N2
      REAL               BIGNUM, COND, EPS, LNRM, RNRM, SCALE, SMLNUM
      COMPLEX            YHAX, YHBX
*     ..
*     .. Local Arrays ..
      COMPLEX            DUMMY( 1 ), DUMMY1( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SCNRM2, SLAMCH, SLAPY2, SROUNDUP_LWORK
      COMPLEX            CDOTC
      EXTERNAL           LSAME, SCNRM2, SLAMCH, SLAPY2, SROUNDUP_LWORK,
     $                   CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMV, CLACPY, CTGEXC, CTGSYL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTDF = LSAME( JOB, 'V' ) .OR. WANTBH
*
      SOMCON = LSAME( HOWMNY, 'S' )
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( .NOT.WANTS .AND. .NOT.WANTDF ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( WANTS .AND. LDVL.LT.N ) THEN
         INFO = -10
      ELSE IF( WANTS .AND. LDVR.LT.N ) THEN
         INFO = -12
      ELSE
*
*        Set M to the number of eigenpairs for which condition numbers
*        are required, and test MM.
*
         IF( SOMCON ) THEN
            M = 0
            DO 10 K = 1, N
               IF( SELECT( K ) )
     $            M = M + 1
   10       CONTINUE
         ELSE
            M = N
         END IF
*
         IF( N.EQ.0 ) THEN
            LWMIN = 1
         ELSE IF( LSAME( JOB, 'V' ) .OR. LSAME( JOB, 'B' ) ) THEN
            LWMIN = 2*N*N
         ELSE
            LWMIN = N
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
*
         IF( MM.LT.M ) THEN
            INFO = -15
         ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTGSNA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      KS = 0
      DO 20 K = 1, N
*
*        Determine whether condition numbers are required for the k-th
*        eigenpair.
*
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( K ) )
     $         GO TO 20
         END IF
*
         KS = KS + 1
*
         IF( WANTS ) THEN
*
*           Compute the reciprocal condition number of the k-th
*           eigenvalue.
*
            RNRM = SCNRM2( N, VR( 1, KS ), 1 )
            LNRM = SCNRM2( N, VL( 1, KS ), 1 )
            CALL CGEMV( 'N', N, N, CMPLX( ONE, ZERO ), A, LDA,
     $                  VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 )
            YHAX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            CALL CGEMV( 'N', N, N, CMPLX( ONE, ZERO ), B, LDB,
     $                  VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 )
            YHBX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            COND = SLAPY2( ABS( YHAX ), ABS( YHBX ) )
            IF( COND.EQ.ZERO ) THEN
               S( KS ) = -ONE
            ELSE
               S( KS ) = COND / ( RNRM*LNRM )
            END IF
         END IF
*
         IF( WANTDF ) THEN
            IF( N.EQ.1 ) THEN
               DIF( KS ) = SLAPY2( ABS( A( 1, 1 ) ), ABS( B( 1, 1 ) ) )
            ELSE
*
*              Estimate the reciprocal condition number of the k-th
*              eigenvectors.
*
*              Copy the matrix (A, B) to the array WORK and move the
*              (k,k)th pair to the (1,1) position.
*
               CALL CLACPY( 'Full', N, N, A, LDA, WORK, N )
               CALL CLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
               IFST = K
               ILST = 1
*
               CALL CTGEXC( .FALSE., .FALSE., N, WORK, N, WORK( N*N+1 ),
     $                      N, DUMMY, 1, DUMMY1, 1, IFST, ILST, IERR )
*
               IF( IERR.GT.0 ) THEN
*
*                 Ill-conditioned problem - swap rejected.
*
                  DIF( KS ) = ZERO
               ELSE
*
*                 Reordering successful, solve generalized Sylvester
*                 equation for R and L,
*                            A22 * R - L * A11 = A12
*                            B22 * R - L * B11 = B12,
*                 and compute estimate of Difl[(A11,B11), (A22, B22)].
*
                  N1 = 1
                  N2 = N - N1
                  I = N*N + 1
                  CALL CTGSYL( 'N', IDIFJB, N2, N1, WORK( N*N1+N1+1 ),
     $                         N, WORK, N, WORK( N1+1 ), N,
     $                         WORK( N*N1+N1+I ), N, WORK( I ), N,
     $                         WORK( N1+I ), N, SCALE, DIF( KS ), DUMMY,
     $                         1, IWORK, IERR )
               END IF
            END IF
         END IF
*
   20 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RETURN
*
*     End of CTGSNA
*
      END