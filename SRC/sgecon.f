      SUBROUTINE SGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, LDA, N;
      REAL               ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1;
      REAL               AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      int                ISAMAX;
      REAL               SLAMCH;
      // EXTERNAL LSAME, ISAMAX, SLAMCH, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SLATRS, SRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      HUGEVAL = SLAMCH( 'Overflow' );

      // Test the input parameters.

      INFO = 0;
      ONENRM = NORM == '1' || LSAME( NORM, 'O' );
      if ( !ONENRM && !LSAME( NORM, 'I' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( ANORM < ZERO ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('SGECON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM == ZERO ) {
         return;
      } else if ( SISNAN( ANORM ) ) {
         RCOND = ANORM;
         INFO = -5;
         return;
      } else if ( ANORM > HUGEVAL ) {
         INFO = -5;
         return;
      }

      SMLNUM = SLAMCH( 'Safe minimum' );

      // Estimate the norm of inv(A).

      AINVNM = ZERO;
      NORMIN = 'N';
      if ( ONENRM ) {
         KASE1 = 1;
      } else {
         KASE1 = 2;
      }
      KASE = 0;
      } // 10
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(L).

            slatrs('Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO );

            // Multiply by inv(U).

            slatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO );
         } else {

            // Multiply by inv(U**T).

            slatrs('Upper', 'Transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO );

            // Multiply by inv(L**T).

            slatrs('Lower', 'Transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO );
         }

         // Divide X by 1/(SL*SU) if doing so will not cause overflow.

         SCALE = SL*SU;
         NORMIN = 'Y';
         if ( SCALE != ONE ) {
            IX = ISAMAX( N, WORK, 1 );
            if( SCALE < ABS( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 20;
            srscl(N, SCALE, WORK, 1 );
         }
         GO TO 10;
      }

      // Compute the estimate of the reciprocal condition number.

      if ( AINVNM != ZERO ) {
         RCOND = ( ONE / AINVNM ) / ANORM;
      } else {
         INFO = 1;
         return;
      }

      // Check for NaNs and Infs

      if( SISNAN( RCOND ) || RCOND > HUGEVAL ) INFO = 1;

      } // 20
      return;

      // End of SGECON

      }
