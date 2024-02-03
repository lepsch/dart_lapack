      SUBROUTINE CHPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            AP( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IP, KASE;
      REAL               AINVNM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPTRS, CLACN2, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( ANORM.LT.ZERO ) {
         INFO = -5
      }
      if ( INFO != 0 ) {
         xerbla('CHPCON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N == 0 ) {
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
            IF( IPIV( I ).GT.0 && AP( IP ) == ZERO ) RETURN
            IP = IP - I
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         IP = 1
         for (I = 1; I <= N; I++) { // 20
            IF( IPIV( I ).GT.0 && AP( IP ) == ZERO ) RETURN
            IP = IP + N - I + 1
         } // 20
      }

      // Estimate the 1-norm of the inverse.

      KASE = 0
      } // 30
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {

         // Multiply by inv(L*D*L**H) or inv(U*D*U**H).

         chptrs(UPLO, N, 1, AP, IPIV, WORK, N, INFO );
         GO TO 30
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      RETURN

      // End of CHPCON

      }
