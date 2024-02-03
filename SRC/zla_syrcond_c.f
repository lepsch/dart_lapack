      double           FUNCTION ZLA_SYRCOND_C( UPLO, N, A, LDA, AF, LDAF, IPIV, C, CAPPLY, INFO, WORK, RWORK );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      bool               CAPPLY;
      int                N, LDA, LDAF, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * )
      double             C( * ), RWORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                KASE;
      double             AINVNM, ANORM, TMP;
      int                I, J;
      bool               UP, UPPER;
      COMPLEX*16         ZDUM
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACN2, ZSYTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Statement Functions ..
      double           CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      ZLA_SYRCOND_C = 0.0D+0

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLA_SYRCOND_C', -INFO );
         RETURN
      }
      UP = false;
      IF ( LSAME( UPLO, 'U' ) ) UP = true;

      // Compute norm of op(A)*op2(C).

      ANORM = 0.0D+0
      if ( UP ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CAPPLY ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( J, I ) ) / C( J )
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( I, J ) ) / C( J )
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( J, I ) )
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( I, J ) )
               }
            }
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0D+0
            if ( CAPPLY ) {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( I, J ) ) / C( J )
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( J, I ) ) / C( J )
               }
            } else {
               for (J = 1; J <= I; J++) {
                  TMP = TMP + CABS1( A( I, J ) )
               }
               for (J = I+1; J <= N; J++) {
                  TMP = TMP + CABS1( A( J, I ) )
               }
            }
            RWORK( I ) = TMP
            ANORM = MAX( ANORM, TMP )
         }
      }

      // Quick return if possible.

      if ( N == 0 ) {
         ZLA_SYRCOND_C = 1.0D+0
         RETURN
      } else if ( ANORM == 0.0D+0 ) {
         RETURN
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0D+0

      KASE = 0
      } // 10
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE == 2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            }

            if ( UP ) {
               zsytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               zsytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               }
            }
         } else {

            // Multiply by inv(C**T).

            if ( CAPPLY ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               }
            }

            if ( UP ) {
               zsytrs('U', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               zsytrs('L', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * RWORK( I )
            }
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM .NE. 0.0D+0) ZLA_SYRCOND_C = 1.0D+0 / AINVNM;

      RETURN

      // End of ZLA_SYRCOND_C

      }
