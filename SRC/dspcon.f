      SUBROUTINE DSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             AP( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IP, KASE;
      double             AINVNM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DSPTRS, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( ANORM.LT.ZERO ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('DSPCON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N.EQ.0 ) {
         RCOND = ONE
         RETURN
      } else if ( ANORM.LE.ZERO ) {
         RETURN
      }

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         IP = N*( N+1 ) / 2
         DO 10 I = N, 1, -1
            IF( IPIV( I ).GT.0 .AND. AP( IP ).EQ.ZERO ) RETURN
            IP = IP - I
   10    CONTINUE
      } else {

         // Lower triangular storage: examine D from top to bottom.

         IP = 1
         DO 20 I = 1, N
            IF( IPIV( I ).GT.0 .AND. AP( IP ).EQ.ZERO ) RETURN
            IP = IP + N - I + 1
   20    CONTINUE
      }

      // Estimate the 1-norm of the inverse.

      KASE = 0
   30 CONTINUE
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {

         // Multiply by inv(L*D*L**T) or inv(U*D*U**T).

         dsptrs(UPLO, N, 1, AP, IPIV, WORK, N, INFO );
         GO TO 30
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

      RETURN

      // End of DSPCON

      }
