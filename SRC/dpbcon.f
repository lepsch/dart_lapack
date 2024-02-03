      SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AB( LDAB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      String             NORMIN;
      int                IX, KASE;
      double             AINVNM, SCALE, SCALEL, SCALEU, SMLNUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACN2, DLATBS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
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
      if ( INFO.NE.0 ) {
         xerbla('DPBCON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N.EQ.0 ) {
         RCOND = ONE
         RETURN
      } else if ( ANORM.EQ.ZERO ) {
         RETURN
      }

      SMLNUM = DLAMCH( 'Safe minimum' )

      // Estimate the 1-norm of the inverse.

      KASE = 0
      NORMIN = 'N'
      } // 10
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( UPPER ) {

            // Multiply by inv(U**T).

            dlatbs('Upper', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y'

            // Multiply by inv(U).

            dlatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         } else {

            // Multiply by inv(L).

            dlatbs('Lower', 'No transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEL, WORK( 2*N+1 ), INFO );
            NORMIN = 'Y'

            // Multiply by inv(L**T).

            dlatbs('Lower', 'Transpose', 'Non-unit', NORMIN, N, KD, AB, LDAB, WORK, SCALEU, WORK( 2*N+1 ), INFO );
         }

         // Multiply by 1/SCALE if doing so will not cause overflow.

         SCALE = SCALEL*SCALEU
         if ( SCALE.NE.ONE ) {
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO ) GO TO 20
            drscl(N, SCALE, WORK, 1 );
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM.NE.ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 20

      RETURN

      // End of DPBCON

      }
