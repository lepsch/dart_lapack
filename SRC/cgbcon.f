      SUBROUTINE CGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                INFO, KL, KU, LDAB, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               RWORK( * )
      COMPLEX            AB( LDAB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LNOTI, ONENRM;
      String             NORMIN;
      int                IX, J, JP, KASE, KASE1, KD, LM;
      REAL               AINVNM, SCALE, SMLNUM
      COMPLEX            T, ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH
      COMPLEX            CDOTC
      // EXTERNAL LSAME, ICAMAX, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CLACN2, CLATBS, CSRSCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      ONENRM = NORM == '1' .OR. LSAME( NORM, 'O' )
      if ( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KL.LT.0 ) {
         INFO = -3
      } else if ( KU.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.2*KL+KU+1 ) {
         INFO = -6
      } else if ( ANORM.LT.ZERO ) {
         INFO = -8
      }
      if ( INFO != 0 ) {
         xerbla('CGBCON', -INFO );
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

      // Estimate the norm of inv(A).

      AINVNM = ZERO
      NORMIN = 'N'
      if ( ONENRM ) {
         KASE1 = 1
      } else {
         KASE1 = 2
      }
      KD = KL + KU + 1
      LNOTI = KL.GT.0
      KASE = 0
      } // 10
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
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
                  caxpy(LM, -T, AB( KD+1, J ), 1, WORK( J+1 ), 1 );
               } // 20
            }

            // Multiply by inv(U).

            clatbs('Upper', 'No transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, RWORK, INFO );
         } else {

            // Multiply by inv(U**H).

            clatbs('Upper', 'Conjugate transpose', 'Non-unit', NORMIN, N, KL+KU, AB, LDAB, WORK, SCALE, RWORK, INFO );

            // Multiply by inv(L**H).

            if ( LNOTI ) {
               DO 30 J = N - 1, 1, -1
                  LM = MIN( KL, N-J )
                  WORK( J ) = WORK( J ) - CDOTC( LM, AB( KD+1, J ), 1, WORK( J+1 ), 1 )
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
            IX = ICAMAX( N, WORK, 1 )
            IF( SCALE.LT.CABS1( WORK( IX ) )*SMLNUM .OR. SCALE == ZERO ) GO TO 40
            csrscl(N, SCALE, WORK, 1 );
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      } // 40
      RETURN

      // End of CGBCON

      }
