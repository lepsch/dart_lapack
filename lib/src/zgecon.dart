      void zgecon(NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, LDA, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL;
      Complex         ZDUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame, DISNAN;
      //- int                IZAMAX;
      //- double             DLAMCH;
      // EXTERNAL lsame, IZAMAX, DLAMCH, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDRSCL, ZLACN2, ZLATRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      HUGEVAL = DLAMCH( 'Overflow' );

      // Test the input parameters.

      INFO = 0;
      ONENRM = NORM == '1' || lsame( NORM, 'O' );
      if ( !ONENRM && !lsame( NORM, 'I' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( ANORM < ZERO ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('ZGECON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM == ZERO ) {
         return;
      } else if ( DISNAN( ANORM ) ) {
         RCOND = ANORM;
         INFO = -5;
         return;
      } else if ( ANORM > HUGEVAL ) {
         INFO = -5;
         return;
      }

      SMLNUM = DLAMCH( 'Safe minimum' );

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
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(L).

            zlatrs('Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO );

            // Multiply by inv(U).

            zlatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO );
         } else {

            // Multiply by inv(U**H).

            zlatrs('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO );

            // Multiply by inv(L**H).

            zlatrs('Lower', 'Conjugate transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO );
         }

         // Divide X by 1/(SL*SU) if doing so will not cause overflow.

         SCALE = SL*SU;
         NORMIN = 'Y';
         if ( SCALE != ONE ) {
            IX = IZAMAX( N, WORK, 1 );
            if( SCALE < CABS1( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 20;
            zdrscl(N, SCALE, WORK, 1 );
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

      if( DISNAN( RCOND ) || RCOND > HUGEVAL ) INFO = 1;

      } // 20
      return;
      }