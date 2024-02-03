      REAL FUNCTION SLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, CMODE, C, INFO, WORK, IWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, LDA, LDAF, INFO, CMODE;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )
*    ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J;
      REAL               AINVNM, TMP
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SGETRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      SLA_GERCOND = 0.0

      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      if ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT. LSAME(TRANS, 'C') ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -6
      }
      if ( INFO != 0 ) {
         xerbla('SLA_GERCOND', -INFO );
         RETURN
      }
      if ( N == 0 ) {
         SLA_GERCOND = 1.0
         RETURN
      }

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      if (NOTRANS) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0
            if ( CMODE == 1 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) * C( J ) )
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) )
               }
            } else {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( I, J ) / C( J ) )
               }
            }
            WORK( 2*N+I ) = TMP
         }
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0
            if ( CMODE == 1 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) * C( J ) )
               }
            } else if ( CMODE == 0 ) {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) )
               }
            } else {
               for (J = 1; J <= N; J++) {
                  TMP = TMP + ABS( A( J, I ) / C( J ) )
               }
            }
            WORK( 2*N+I ) = TMP
         }
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0

      KASE = 0
      } // 10
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == 2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK(I) = WORK(I) * WORK(2*N+I)
            }

            if (NOTRANS) {
               sgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               sgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I )
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               }
            }
         } else {

            // Multiply by inv(C**T).

            if ( CMODE == 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I )
               }
            } else if ( CMODE == -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               }
            }

            if (NOTRANS) {
               sgetrs('Transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            } else {
               sgetrs('No transpose', N, 1, AF, LDAF, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            }
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != 0.0) SLA_GERCOND = ( 1.0 / AINVNM );

      RETURN

      // End of SLA_GERCOND

      }
