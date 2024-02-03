      SUBROUTINE ZTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             HOWMNY, JOB;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N;
*     ..
*     .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      double             DIF( * ), S( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      int                IDIFJB;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, IDIFJB = 3 )
*     ..
*     .. Local Scalars ..
      bool               LQUERY, SOMCON, WANTBH, WANTDF, WANTS;
      int                I, IERR, IFST, ILST, K, KS, LWMIN, N1, N2;
      double             BIGNUM, COND, EPS, LNRM, RNRM, SCALE, SMLNUM;
      COMPLEX*16         YHAX, YHBX
*     ..
*     .. Local Arrays ..
      COMPLEX*16         DUMMY( 1 ), DUMMY1( 1 )
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLAPY2, DZNRM2;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, DLAMCH, DLAPY2, DZNRM2, ZDOTC
*     ..
*     .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZLACPY, ZTGEXC, ZTGSYL
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, MAX
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
               IF( SELECT( K ) ) M = M + 1
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
         WORK( 1 ) = LWMIN
*
         IF( MM.LT.M ) THEN
            INFO = -15
         ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSNA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      KS = 0
      DO 20 K = 1, N
*
*        Determine whether condition numbers are required for the k-th
*        eigenpair.
*
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( K ) ) GO TO 20
         END IF
*
         KS = KS + 1
*
         IF( WANTS ) THEN
*
*           Compute the reciprocal condition number of the k-th
*           eigenvalue.
*
            RNRM = DZNRM2( N, VR( 1, KS ), 1 )
            LNRM = DZNRM2( N, VL( 1, KS ), 1 )
            CALL ZGEMV( 'N', N, N, DCMPLX( ONE, ZERO ), A, LDA, VR( 1, KS ), 1, DCMPLX( ZERO, ZERO ), WORK, 1 )
            YHAX = ZDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            CALL ZGEMV( 'N', N, N, DCMPLX( ONE, ZERO ), B, LDB, VR( 1, KS ), 1, DCMPLX( ZERO, ZERO ), WORK, 1 )
            YHBX = ZDOTC( N, WORK, 1, VL( 1, KS ), 1 )
            COND = DLAPY2( ABS( YHAX ), ABS( YHBX ) )
            IF( COND.EQ.ZERO ) THEN
               S( KS ) = -ONE
            ELSE
               S( KS ) = COND / ( RNRM*LNRM )
            END IF
         END IF
*
         IF( WANTDF ) THEN
            IF( N.EQ.1 ) THEN
               DIF( KS ) = DLAPY2( ABS( A( 1, 1 ) ), ABS( B( 1, 1 ) ) )
            ELSE
*
*              Estimate the reciprocal condition number of the k-th
*              eigenvectors.
*
*              Copy the matrix (A, B) to the array WORK and move the
*              (k,k)th pair to the (1,1) position.
*
               CALL ZLACPY( 'Full', N, N, A, LDA, WORK, N )
               CALL ZLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
               IFST = K
               ILST = 1
*
               CALL ZTGEXC( .FALSE., .FALSE., N, WORK, N, WORK( N*N+1 ), N, DUMMY, 1, DUMMY1, 1, IFST, ILST, IERR )
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
                  CALL ZTGSYL( 'N', IDIFJB, N2, N1, WORK( N*N1+N1+1 ), N, WORK, N, WORK( N1+1 ), N, WORK( N*N1+N1+I ), N, WORK( I ), N, WORK( N1+I ), N, SCALE, DIF( KS ), DUMMY, 1, IWORK, IERR )
               END IF
            END IF
         END IF
*
   20 CONTINUE
      WORK( 1 ) = LWMIN
      RETURN
*
*     End of ZTGSNA
*
      END
