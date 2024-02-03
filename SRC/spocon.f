      SUBROUTINE SPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      REAL               ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      REAL               AINVNM, SCALE, SCALEL, SCALEU, SMLNUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SLAMCH;
      // EXTERNAL LSAME, ISAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLATRS, SRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      } else if ( ANORM < ZERO ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('SPOCON', -INFO );
         RETURN;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         RETURN;
      } else if ( ANORM == ZERO ) {
         RETURN;
      }

      SMLNUM = SLAMCH( 'Safe minimum' );

      // Estimate the 1-norm of inv(A).

      KASE = 0;
      NORMIN = 'N';
      } // 10
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( UPPER ) {

            // Multiply by inv(U**T).

            slatrs('Upper', 'Transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y';

            // Multiply by inv(U).

            slatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         } else {

            // Multiply by inv(L).

            slatrs('Lower', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y';

            // Multiply by inv(L**T).

            slatrs('Lower', 'Transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         }

         // Multiply by 1/SCALE if doing so will not cause overflow.

         SCALE = SCALEL*SCALEU;
         if ( SCALE != ONE ) {
            IX = ISAMAX( N, WORK, 1 );
            IF( SCALE < ABS( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 20;
            srscl(N, SCALE, WORK, 1 );
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 20
      RETURN;

      // End of SPOCON

      }
