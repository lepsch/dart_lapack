      SUBROUTINE DSPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK, SWORK, ITER, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                INFO, ITER, LDA, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               SWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( N, * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      LOGICAL            DOITREF
      PARAMETER          ( DOITREF = .TRUE. )
*
      int                ITERMAX
      PARAMETER          ( ITERMAX = 30 )
*
      DOUBLE PRECISION   BWDMAX
      PARAMETER          ( BWDMAX = 1.0E+00 )
*
      DOUBLE PRECISION   NEGONE, ONE
      PARAMETER          ( NEGONE = -1.0D+0, ONE = 1.0D+0 )
*
*     .. Local Scalars ..
      int                I, IITER, PTSA, PTSX
      DOUBLE PRECISION   ANRM, CTE, EPS, RNRM, XNRM
*
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSYMM, DLACPY, DLAT2S, DLAG2S, SLAG2D, SPOTRF, SPOTRS, DPOTRF, DPOTRS, XERBLA
*     ..
*     .. External Functions ..
      int                IDAMAX
      DOUBLE PRECISION   DLAMCH, DLANSY
      LOGICAL            LSAME
      EXTERNAL           IDAMAX, DLAMCH, DLANSY, LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      ITER = 0
*
*     Test the input parameters.
*
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSPOSV', -INFO )
         RETURN
      END IF
*
*     Quick return if (N.EQ.0).
*
      IF( N.EQ.0 ) RETURN
*
*     Skip single precision iterative refinement if a priori slower
*     than double precision factorization.
*
      IF( .NOT.DOITREF ) THEN
         ITER = -1
         GO TO 40
      END IF
*
*     Compute some constants.
*
      ANRM = DLANSY( 'I', UPLO, N, A, LDA, WORK )
      EPS = DLAMCH( 'Epsilon' )
      CTE = ANRM*EPS*SQRT( DBLE( N ) )*BWDMAX
*
*     Set the indices PTSA, PTSX for referencing SA and SX in SWORK.
*
      PTSA = 1
      PTSX = PTSA + N*N
*
*     Convert B from double precision to single precision and store the
*     result in SX.
*
      CALL DLAG2S( N, NRHS, B, LDB, SWORK( PTSX ), N, INFO )
*
      IF( INFO.NE.0 ) THEN
         ITER = -2
         GO TO 40
      END IF
*
*     Convert A from double precision to single precision and store the
*     result in SA.
*
      CALL DLAT2S( UPLO, N, A, LDA, SWORK( PTSA ), N, INFO )
*
      IF( INFO.NE.0 ) THEN
         ITER = -2
         GO TO 40
      END IF
*
*     Compute the Cholesky factorization of SA.
*
      CALL SPOTRF( UPLO, N, SWORK( PTSA ), N, INFO )
*
      IF( INFO.NE.0 ) THEN
         ITER = -3
         GO TO 40
      END IF
*
*     Solve the system SA*SX = SB.
*
      CALL SPOTRS( UPLO, N, NRHS, SWORK( PTSA ), N, SWORK( PTSX ), N, INFO )
*
*     Convert SX back to double precision
*
      CALL SLAG2D( N, NRHS, SWORK( PTSX ), N, X, LDX, INFO )
*
*     Compute R = B - AX (R is WORK).
*
      CALL DLACPY( 'All', N, NRHS, B, LDB, WORK, N )
*
      CALL DSYMM( 'Left', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N )
*
*     Check whether the NRHS normwise backward errors satisfy the
*     stopping criterion. If yes, set ITER=0 and return.
*
      DO I = 1, NRHS
         XNRM = ABS( X( IDAMAX( N, X( 1, I ), 1 ), I ) )
         RNRM = ABS( WORK( IDAMAX( N, WORK( 1, I ), 1 ), I ) )
         IF( RNRM.GT.XNRM*CTE ) GO TO 10
      END DO
*
*     If we are here, the NRHS normwise backward errors satisfy the
*     stopping criterion. We are good to exit.
*
      ITER = 0
      RETURN
*
   10 CONTINUE
*
      DO 30 IITER = 1, ITERMAX
*
*        Convert R (in WORK) from double precision to single precision
*        and store the result in SX.
*
         CALL DLAG2S( N, NRHS, WORK, N, SWORK( PTSX ), N, INFO )
*
         IF( INFO.NE.0 ) THEN
            ITER = -2
            GO TO 40
         END IF
*
*        Solve the system SA*SX = SR.
*
         CALL SPOTRS( UPLO, N, NRHS, SWORK( PTSA ), N, SWORK( PTSX ), N, INFO )
*
*        Convert SX back to double precision and update the current
*        iterate.
*
         CALL SLAG2D( N, NRHS, SWORK( PTSX ), N, WORK, N, INFO )
*
         DO I = 1, NRHS
            CALL DAXPY( N, ONE, WORK( 1, I ), 1, X( 1, I ), 1 )
         END DO
*
*        Compute R = B - AX (R is WORK).
*
         CALL DLACPY( 'All', N, NRHS, B, LDB, WORK, N )
*
         CALL DSYMM( 'L', UPLO, N, NRHS, NEGONE, A, LDA, X, LDX, ONE, WORK, N )
*
*        Check whether the NRHS normwise backward errors satisfy the
*        stopping criterion. If yes, set ITER=IITER>0 and return.
*
         DO I = 1, NRHS
            XNRM = ABS( X( IDAMAX( N, X( 1, I ), 1 ), I ) )
            RNRM = ABS( WORK( IDAMAX( N, WORK( 1, I ), 1 ), I ) )
            IF( RNRM.GT.XNRM*CTE ) GO TO 20
         END DO
*
*        If we are here, the NRHS normwise backward errors satisfy the
*        stopping criterion, we are good to exit.
*
         ITER = IITER
*
         RETURN
*
   20    CONTINUE
*
   30 CONTINUE
*
*     If we are at this place of the code, this is because we have
*     performed ITER=ITERMAX iterations and never satisfied the
*     stopping criterion, set up the ITER flag accordingly and follow
*     up on double precision routine.
*
      ITER = -ITERMAX - 1
*
   40 CONTINUE
*
*     Single-precision iterative refinement failed to converge to a
*     satisfactory solution, so we resort to double precision.
*
      CALL DPOTRF( UPLO, N, A, LDA, INFO )
*
      IF( INFO.NE.0 ) RETURN
*
      CALL DLACPY( 'All', N, NRHS, B, LDB, X, LDX )
      CALL DPOTRS( UPLO, N, NRHS, A, LDA, X, LDX, INFO )
*
      RETURN
*
*     End of DSPOSV
*
      END
