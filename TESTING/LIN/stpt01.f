      SUBROUTINE STPT01( UPLO, DIAG, N, AP, AINVP, RCOND, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                N;
      REAL               RCOND, RESID
      // ..
      // .. Array Arguments ..
      REAL               AINVP( * ), AP( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UNITD;
      int                J, JC;
      REAL               AINVNM, ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANTP
      // EXTERNAL LSAME, SLAMCH, SLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RCOND = ONE
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANTP( '1', UPLO, DIAG, N, AP, WORK )
      AINVNM = SLANTP( '1', UPLO, DIAG, N, AINVP, WORK )
      if ( ANORM.LE.ZERO || AINVNM.LE.ZERO ) {
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      }
      RCOND = ( ONE / ANORM ) / AINVNM

      // Compute A * AINV, overwriting AINV.

      UNITD = LSAME( DIAG, 'U' )
      if ( LSAME( UPLO, 'U' ) ) {
         JC = 1
         for (J = 1; J <= N; J++) { // 10
            if (UNITD) AINVP( JC+J-1 ) = ONE;

            // Form the j-th column of A*AINV

            stpmv('Upper', 'No transpose', DIAG, J, AP, AINVP( JC ), 1 );

            // Subtract 1 from the diagonal

            AINVP( JC+J-1 ) = AINVP( JC+J-1 ) - ONE
            JC = JC + J
         } // 10
      } else {
         JC = 1
         for (J = 1; J <= N; J++) { // 20
            if (UNITD) AINVP( JC ) = ONE;

            // Form the j-th column of A*AINV

            stpmv('Lower', 'No transpose', DIAG, N-J+1, AP( JC ), AINVP( JC ), 1 );

            // Subtract 1 from the diagonal

            AINVP( JC ) = AINVP( JC ) - ONE
            JC = JC + N - J + 1
         } // 20
      }

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = SLANTP( '1', UPLO, 'Non-unit', N, AINVP, WORK )

      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS

      RETURN

      // End of STPT01

      }
