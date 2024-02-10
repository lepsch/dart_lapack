      void cgecon(NORM, N, final Matrix<double> A, final int LDA, ANORM, RCOND, WORK, final Array<double> RWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM;
      int                INFO, LDA, N;
      double               ANORM, RCOND;
      double               RWORK( * );
      Complex            A( LDA, * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double               AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL;
      Complex            ZDUM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      //- int                ICAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ICAMAX, SLAMCH, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACN2, CLATRS, CSRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();

      HUGEVAL = SLAMCH( 'Overflow' );

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
         xerbla('CGECON', -INFO );
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
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(L).

            clatrs('Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO );

            // Multiply by inv(U).

            clatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO );
         } else {

            // Multiply by inv(U**H).

            clatrs('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, RWORK( N+1 ), INFO );

            // Multiply by inv(L**H).

            clatrs('Lower', 'Conjugate transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, RWORK, INFO );
         }

         // Divide X by 1/(SL*SU) if doing so will not cause overflow.

         SCALE = SL*SU;
         NORMIN = 'Y';
         if ( SCALE != ONE ) {
            IX = ICAMAX( N, WORK, 1 );
            if( SCALE < CABS1( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 20;
            csrscl(N, SCALE, WORK, 1 );
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
      }
