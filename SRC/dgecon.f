      SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, LDA, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               ONENRM;
      String             NORMIN;
      int                IX, KASE, KASE1;
      double             AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLATRS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      HUGEVAL = DLAMCH( 'Overflow' )

      // Test the input parameters.

      INFO = 0
      ONENRM = NORM == '1' || LSAME( NORM, 'O' )
      if ( !ONENRM && !LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      } else if ( ANORM < ZERO ) {
         INFO = -5
      }
      if ( INFO != 0 ) {
         xerbla('DGECON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N == 0 ) {
         RCOND = ONE
         RETURN
      } else if ( ANORM == ZERO ) {
         RETURN
      } else if ( DISNAN( ANORM ) ) {
         RCOND = ANORM
         INFO = -5
         RETURN
      } else if ( ANORM > HUGEVAL ) {
         INFO = -5
         RETURN
      }

      SMLNUM = DLAMCH( 'Safe minimum' )

      // Estimate the norm of inv(A).

      AINVNM = ZERO
      NORMIN = 'N'
      if ( ONENRM ) {
         KASE1 = 1
      } else {
         KASE1 = 2
      }
      KASE = 0
      } // 10
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(L).

            dlatrs('Lower', 'No transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO );

            // Multiply by inv(U).

            dlatrs('Upper', 'No transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO );
         } else {

            // Multiply by inv(U**T).

            dlatrs('Upper', 'Transpose', 'Non-unit', NORMIN, N, A, LDA, WORK, SU, WORK( 3*N+1 ), INFO );

            // Multiply by inv(L**T).

            dlatrs('Lower', 'Transpose', 'Unit', NORMIN, N, A, LDA, WORK, SL, WORK( 2*N+1 ), INFO );
         }

         // Divide X by 1/(SL*SU) if doing so will not cause overflow.

         SCALE = SL*SU
         NORMIN = 'Y'
         if ( SCALE != ONE ) {
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE < ABS( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 20
            drscl(N, SCALE, WORK, 1 );
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      if ( AINVNM != ZERO ) {
         RCOND = ( ONE / AINVNM ) / ANORM
      } else {
         INFO = 1
         RETURN
      }

      // Check for NaNs and Infs

      IF( DISNAN( RCOND ) || RCOND > HUGEVAL ) INFO = 1

      } // 20
      RETURN

      // End of DGECON

      }
