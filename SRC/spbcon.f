      SUBROUTINE SPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AB( LDAB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      REAL               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ISAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLATBS, SRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KD.LT.0 ) {
         INFO = -3
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -5
      } else if ( ANORM.LT.ZERO ) {
         INFO = -6
      }
      if ( INFO != 0 ) {
         xerbla('SPBCON', -INFO );
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

      SMLNUM = SLAMCH( 'Safe minimum' )

      // Estimate the 1-norm of the inverse.

      KASE = 0
      NORMIN = 'N'
      } // 10
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( UPPER ) {

            // Multiply by inv(U**T).

            slatbs('Upper', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y'

            // Multiply by inv(U).

            slatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         } else {

            // Multiply by inv(L).

            slatbs('Lower', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y'

            // Multiply by inv(L**T).

            slatbs('Lower', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         }

         // Multiply by 1/SCALE if doing so will not cause overflow.

         SCALE = SCALEL*SCALEU
         if ( SCALE != ONE ) {
            IX = ISAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 20
            srscl(N, SCALE, WORK, 1 );
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 20

      RETURN

      // End of SPBCON

      }
