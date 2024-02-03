      SUBROUTINE DGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             D( * ), DL( * ), DU( * ), DU2( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               ONENRM;
      int                I, KASE, KASE1;
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
      // EXTERNAL DGTTRS, DLACN2, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      ONENRM = NORM == '1' || LSAME( NORM, 'O' )
      if ( .NOT.ONENRM && .NOT.LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( ANORM < ZERO ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('DGTCON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N == 0 ) {
         RCOND = ONE
         RETURN
      } else if ( ANORM == ZERO ) {
         RETURN
      }

      // Check that D(1:N) is non-zero.

      for (I = 1; I <= N; I++) { // 10
         IF( D( I ) == ZERO ) RETURN
      } // 10

      AINVNM = ZERO
      if ( ONENRM ) {
         KASE1 = 1
      } else {
         KASE1 = 2
      }
      KASE = 0
      } // 20
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(U)*inv(L).

            dgttrs('No transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO );
         } else {

            // Multiply by inv(L**T)*inv(U**T).

            dgttrs('Transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO );
         }
         GO TO 20
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      RETURN

      // End of DGTCON

      }
