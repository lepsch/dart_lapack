      SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, KL, KU, LDAB, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             AB( LDAB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LNOTI, ONENRM;
      String             NORMIN;
      int                IX, J, JP, KASE, KASE1, KD, LM;
      double             AINVNM, SCALE, SMLNUM, T;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DDOT, DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DDOT, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLACN2, DLATBS, DRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      ONENRM = NORM == '1' || LSAME( NORM, 'O' )
      if ( .NOT.ONENRM && .NOT.LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( KL < 0 ) {
         INFO = -3
      } else if ( KU < 0 ) {
         INFO = -4
      } else if ( LDAB < 2*KL+KU+1 ) {
         INFO = -6
      } else if ( ANORM < ZERO ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('DGBCON', -INFO );
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

      SMLNUM = DLAMCH( 'Safe minimum' )

      // Estimate the norm of inv(A).

      AINVNM = ZERO
      NORMIN = 'N'
      if ( ONENRM ) {
         KASE1 = 1
      } else {
         KASE1 = 2
      }
      KD = KL + KU + 1
      LNOTI = KL > 0
      KASE = 0
      } // 10
      dlacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(L).

            if ( LNOTI ) {
               for (J = 1; J <= N - 1; J++) { // 20
                  LM = MIN( KL, N-J )
                  JP = IPIV( J )
                  T = WORK( JP )
                  if ( JP != J ) {
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  }
                  daxpy(LM, -T, AB( KD+1, J ), 1, WORK( J+1 ), 1 );
               } // 20
            }

            // Multiply by inv(U).

            dlatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );
         } else {

            // Multiply by inv(U**T).

            dlatbs('Upper', 'Transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, WORK( 2*N+1 ), INFO );

            // Multiply by inv(L**T).

            if ( LNOTI ) {
               DO 30 J = N - 1, 1, -1
                  LM = MIN( KL, N-J )
                  WORK( J ) = WORK( J ) - DDOT( LM, AB( KD+1, J ), 1, WORK( J+1 ), 1 )
                  JP = IPIV( J )
                  if ( JP != J ) {
                     T = WORK( JP )
                     WORK( JP ) = WORK( J )
                     WORK( J ) = T
                  }
               } // 30
            }
         }

         // Divide X by 1/SCALE if doing so will not cause overflow.

         NORMIN = 'Y'
         if ( SCALE != ONE ) {
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE < ABS( WORK( IX ) )*SMLNUM || SCALE == ZERO ) GO TO 40
            drscl(N, SCALE, WORK, 1 );
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 40
      RETURN

      // End of DGBCON

      }
